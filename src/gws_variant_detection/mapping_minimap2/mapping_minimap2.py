# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os
import glob
import pandas as pd
import numpy as np
import plotly.graph_objects as go

from gws_core import (
    ConfigParams, ConfigSpecs, File, InputSpec, Folder, ShellProxy,
    InputSpecs, IntParam, OutputSpec, OutputSpecs, Task,
    TaskInputs, TaskOutputs, task_decorator, FloatParam,
    PlotlyResource, ResourceSet
)

from ..base_env.Alignment_env import AlignmentShellProxyHelper

@task_decorator(
    "ReadsMapping",
    human_name="NanoReadsMapping",
    short_description=(
        "Minimap2 mapping of ONT reads with primary-only filtering, optional low-tail removal "
        "(≤ q), exact-X downsampling by QNAME; returns BAMs, HTML plots, and a ResourceSet of Plotly distributions"
    )
)
class NanoMapping(Task):
    """
    Steps:
      1) Mapping with Minimap2 (preset: map-ont), excluding secondary and supplementary alignments.
      2) Optional: removal of the shortest reads (≤ low_tail_q fraction of the length distribution).
      3) Exact downsampling by QNAME (after filtering) if 'target_reads' > 0.
      4) Per-sample outputs:
           - <sample>.sorted.bam (+ index)
           - HTML plots (histogram + KDE) for Read Length and Read Coverage (final primary alignments)
           - Measurement CSVs: <sample>.final_lengths.csv (ReadLength), <sample>.final_coverage.csv (CoveragePct)
      5) Additionally, this Task builds a ResourceSet of Plotly figures
    """
    input_specs = InputSpecs({
        'input_folder': InputSpec(Folder, human_name="Folder containing FASTQ(.gz) files"),
        'input_file':   InputSpec(File,   human_name="Reference FASTA")
    })

    output_specs = OutputSpecs({
        'output_path': OutputSpec(
            Folder,
            human_name="Aligned reads",
            short_description="Final BAMs (<sample>.sorted.bam)"
        ),
        'read_length_plots': OutputSpec(
            Folder,
            human_name="Read-length plots (HTML)",
            short_description="Per-sample HTML histograms/KDE for FINAL BAMs"
        ),
        'distributions': OutputSpec(
            ResourceSet,
            human_name="Read length & coverage (Plotly)",
            short_description="One interactive Plotly figure per sample and per metric (Length / Coverage) with a stats sidebar"
        ),
    })

    config_specs: ConfigSpecs = {
        "num_threads":  IntParam(default_value=15, min_value=2, short_description="Number of threads"),
        "target_reads": IntParam(default_value=5000, min_value=1, short_description="Exact number of primary reads to keep per sample"),
        "low_tail_q":   FloatParam(default_value=0.05, min_value=0.0, max_value=0.9999,
                                   short_description="Fraction to drop from the low-length tail (e.g., 0.05 = drop ≤ 5th percentile; 0 = no filtering)")
    }

    python_file_path = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "_mapping_minimap2.py"
    )

    shell_proxy: ShellProxy = None

    # ---------- Helpers: stats text + histogram/KDE with the same annotation
    @staticmethod
    def _stats_text(values: np.ndarray) -> str:
        s = pd.Series(values, dtype="float64")
        extra = [0.01, 0.025, 0.05, 0.10, 0.90, 0.95, 0.975, 0.99]
        desc = s.describe(percentiles=[0.25, 0.50, 0.75] + extra).to_frame(name="value")
        desc.loc["missing"] = int(len(s.index) - s.notna().sum())

        def _fmt(v):
            try:
                return f"{float(v):.3f}"
            except Exception:
                return str(v)
        return "<br>".join([f"{idx}: {_fmt(val[0])}" for idx, val in desc.iterrows()])

    @staticmethod
    def _kde(values: np.ndarray) -> tuple[np.ndarray, np.ndarray] | None:
        try:
            from sklearn.neighbors import KernelDensity
            from sklearn.model_selection import GridSearchCV
        except Exception:
            return None
        if values.size < 10:
            return None
        x = values.reshape(-1, 1).astype(float)
        std = values.std(ddof=1)
        base_bw = max(1e-6, values.mean() * 0.01) if (std <= 0) else 1.06 * std * (values.shape[0] ** (-1/5))
        grid = np.unique(np.clip(base_bw * np.array([0.25, 0.5, 1.0, 2.0]), 1e-6, None))
        gs = GridSearchCV(KernelDensity(kernel="gaussian"), {"bandwidth": grid},
                          cv=min(5, max(2, values.shape[0] // 5)))
        gs.fit(x)
        kde = gs.best_estimator_
        xv = np.linspace(values.min(), values.max(), 2048)[:, None]
        yv = np.exp(kde.score_samples(xv))
        return xv.ravel(), yv

    @staticmethod
    def _hist_kde_with_stats(values: np.ndarray, title: str, xlab: str, bins: int = 500) -> go.Figure:
        fig = go.Figure()
        fig.add_trace(go.Histogram(x=values, nbinsx=bins, histnorm="probability density",
                                   opacity=0.65, name="Histogram"))
        kde = NanoMapping._kde(values)
        if kde is not None:
            xv, yv = kde
            fig.add_trace(go.Scatter(x=xv, y=yv, mode="lines", name="KDE"))

        stats_txt = NanoMapping._stats_text(values)

        fig.update_layout(
            title=title,
            xaxis_title=xlab,
            yaxis_title="Density",
            template="plotly_white",
            barmode="overlay",
            # >>> same push to the far right + generous right margin
            annotations=[dict(
                x=1.35, y=1, xref="paper", yref="paper",
                xanchor="left",
                text=stats_txt, showarrow=False, align="left", font=dict(size=10)
            )],
            margin=dict(r=560, l=60, t=60, b=60),
        )
        return fig

    @staticmethod
    def _make_length_fig(csv_path: str, sample: str) -> PlotlyResource:
        d = pd.read_csv(csv_path)
        vals = d["ReadLength"].dropna().to_numpy()
        fig = NanoMapping._hist_kde_with_stats(vals, f"{sample} — Read length", "Read length (bp)")
        pr = PlotlyResource(fig)
        pr.name = f"{sample} — Length"
        return pr

    @staticmethod
    def _make_coverage_fig(csv_path: str, sample: str) -> PlotlyResource:
        d = pd.read_csv(csv_path)
        vals = d["CoveragePct"].dropna().to_numpy()
        fig = NanoMapping._hist_kde_with_stats(vals, f"{sample} — Read coverage", "Coverage (%)", bins=200)
        pr = PlotlyResource(fig)
        pr.name = f"{sample} — Coverage"
        return pr

    # ------------------------- run --------------------------------
    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        input_folder: Folder = inputs['input_folder']
        input_file:   File   = inputs['input_file']
        num_threads   = int(params["num_threads"])
        target_reads  = int(params["target_reads"])
        low_tail_q    = float(params["low_tail_q"])

        shell_proxy: ShellProxy = AlignmentShellProxyHelper.create_proxy(self.message_dispatcher)
        result_path = os.path.join(shell_proxy.working_dir, 'Aligned_reads')
        plots_path  = os.path.join(shell_proxy.working_dir, 'read_length_plots')
        os.makedirs(result_path, exist_ok=True)
        os.makedirs(plots_path,  exist_ok=True)

        cmd = (
            f"python3 {self.python_file_path} "
            f"{input_folder.path} {input_file.path} {num_threads} {result_path} {target_reads} {plots_path} {low_tail_q}"
        )
        rc = shell_proxy.run(cmd, shell_mode=True)
        if rc != 0:
            raise RuntimeError("Mapping script failed – check logs.")
        self.shell_proxy = shell_proxy

        # ---- Build ResourceSet from generated CSVs ----
        rs = ResourceSet()

        len_csvs = sorted(glob.glob(os.path.join(plots_path, "*.final_lengths.csv")))
        cov_csvs = sorted(glob.glob(os.path.join(plots_path, "*.final_coverage.csv")))

        def sample_name(fp: str) -> str:
            base = os.path.basename(fp)
            return base.replace(".final_lengths.csv", "").replace(".final_coverage.csv", "")

        for fp in len_csvs:
            s = sample_name(fp)
            rs.add_resource(self._make_length_fig(fp, s), f"{s} — Length")

        for fp in cov_csvs:
            s = sample_name(fp)
            rs.add_resource(self._make_coverage_fig(fp, s), f"{s} — Coverage")

        return {
            'output_path': Folder(result_path),
            'read_length_plots': Folder(plots_path),  # contains HTML + CSV
            'distributions': rs,                      # ResourceSet of PlotlyResource
        }

    def run_after_task(self):
        if self.shell_proxy:
            self.shell_proxy.clean_working_dir()
