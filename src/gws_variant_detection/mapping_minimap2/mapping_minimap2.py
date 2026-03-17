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
           - HTML plots (histogram only) for Read Length and Read Coverage (final primary alignments)
           - Measurement CSVs: <sample>.final_lengths.csv (ReadLength), <sample>.final_coverage.csv (CoveragePct)
      5) Additionally, this Task builds a ResourceSet of Plotly figures (histogram + stats sidebar, no KDE).
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
            short_description="Per-sample HTML histograms for FINAL BAMs"
        ),
        'distributions': OutputSpec(
            ResourceSet,
            human_name="Read length & coverage (Plotly)",
            short_description="One interactive Plotly figure per sample and per metric (Length / Coverage) with a stats sidebar"
        ),
    })

    config_specs = ConfigSpecs({
        "num_threads":  IntParam(default_value=15, min_value=2, short_description="Number of threads"),
        "target_reads": IntParam(default_value=5000, min_value=1, short_description="Exact number of primary reads to keep per sample"),
        "low_tail_q":   FloatParam(default_value=0.05, min_value=0.0, max_value=0.9999,
                                   short_description="Fraction to drop from the low-length tail (e.g., 0.05 = drop ≤ 5th percentile; 0 = no filtering)")
    })

    python_file_path = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "_mapping_minimap2.py"
    )

    shell_proxy: ShellProxy = None

    # ---------- Helpers: stats text (same order as filtering) + histogram-only ----------
    @staticmethod
    def _stats_text(values: np.ndarray) -> str:
        """
        Order: count, mean, std, min, max, 1%, 2.5%, 5%, 10%, 25%, 50%, 75%, 90%, 95%, 97.5%, 99%, missing
        """
        desired_order = [
            "count", "mean", "std", "min", "max",
            "1%", "2.5%", "5%", "10%", "25%", "50%", "75%",
            "90%", "95%", "97.5%", "99%", "missing"
        ]
        s = pd.Series(values).dropna().astype(float)

        def _fmt(v):
            try:
                if isinstance(v, (int, np.integer)) or float(v).is_integer():
                    return f"{int(round(float(v)))}"
                return f"{float(v):.3f}"
            except Exception:
                return str(v)

        if s.empty:
            zero_map = {k: 0 for k in desired_order}
            return "<br>".join(f"{k}: {zero_map[k]}" for k in desired_order)

        desc = s.describe(percentiles=[0.25, 0.50, 0.75])
        qmap = {
            "1%":    s.quantile(0.01),
            "2.5%":  s.quantile(0.025),
            "5%":    s.quantile(0.05),
            "10%":   s.quantile(0.10),
            "90%":   s.quantile(0.90),
            "95%":   s.quantile(0.95),
            "97.5%": s.quantile(0.975),
            "99%":   s.quantile(0.99),
        }
        missing = int(len(values) - s.notna().sum())

        values_map = {
            "count": desc.get("count", 0),
            "mean":  desc.get("mean", 0.0),
            "std":   desc.get("std", 0.0),
            "min":   desc.get("min", 0.0),
            "max":   desc.get("max", 0.0),
            "1%":    qmap["1%"],
            "2.5%":  qmap["2.5%"],
            "5%":    qmap["5%"],
            "10%":   qmap["10%"],
            "25%":   desc.get("25%", 0.0),
            "50%":   desc.get("50%", 0.0),
            "75%":   desc.get("75%", 0.0),
            "90%":   qmap["90%"],
            "95%":   qmap["95%"],
            "97.5%": qmap["97.5%"],
            "99%":   qmap["99%"],
            "missing": missing,
        }
        return "<br>".join(f"{k}: {_fmt(values_map[k])}" for k in desired_order)

    @staticmethod
    def _bins(n: int) -> int:
        # same policy used elsewhere: 10% of N, min 50, max 5000
        return min(max(int(n * 0.10), 50), 5000)

    @staticmethod
    def _hist_with_stats(values: np.ndarray, title: str, xlab: str) -> go.Figure:
        s = pd.Series(values, dtype="float64")
        n = s.notna().sum()
        nbins = NanoMapping._bins(int(n)) if n else 50

        fig = go.Figure()
        fig.add_trace(go.Histogram(
            x=s, nbinsx=nbins, histnorm="probability density",
            opacity=0.75, name="Histogram"
        ))

        stats_txt = NanoMapping._stats_text(values)

        fig.update_layout(
            title=f"{title} — N={int(n)}, bins={nbins}",
            xaxis_title=xlab,
            yaxis_title="Density",
            template="plotly_white",
            legend=dict(x=1.02, y=1, xanchor="left"),
            # Push stats far right and give a large right margin (like filtering)
            annotations=[dict(
                x=1.42, y=1, xref="paper", yref="paper", xanchor="left",
                text=stats_txt, showarrow=False, align="left", font=dict(size=10)
            )],
            margin=dict(r=620, l=60, t=60, b=60),
        )
        return fig

    @staticmethod
    def _make_length_fig(csv_path: str, sample: str) -> PlotlyResource:
        d = pd.read_csv(csv_path)
        vals = d["ReadLength"].dropna().to_numpy()
        fig = NanoMapping._hist_with_stats(vals, f"{sample} — Read length", "Read length (bp)")
        pr = PlotlyResource(fig)
        pr.name = f"{sample} — Length"
        return pr

    @staticmethod
    def _make_coverage_fig(csv_path: str, sample: str) -> PlotlyResource:
        d = pd.read_csv(csv_path)
        vals = d["CoveragePct"].dropna().to_numpy()
        fig = NanoMapping._hist_with_stats(vals, f"{sample} — Read coverage", "Coverage (%)")
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

        # ---- Check output vs input file counts ----
        n_input = len([f for f in os.listdir(input_folder.path) if f.endswith(".fastq.gz")])
        n_output = len([f for f in os.listdir(result_path) if f.endswith(".sorted.bam")])

        if n_output == 0:
            raise RuntimeError(
                f"No output BAM files were produced (0/{n_input} samples). "
                "Please check: (1) the filtering_summary and the coverage distribution plots "
                "generated by the 'NanoReadsFiltering' task, (2) reparametrize the 'NanoReadsFiltering' task "
                "and restart the analysis. "
                "To know more please consult the documentation: https://constellab.community/bricks/gws_variant_detection/latest/doc/variant-detection/38d49814-d0d2-41c5-95f6-4f79ea083b6d"
            )
        elif n_output < n_input:
            self.log_warning_message(
                f"Only {n_output}/{n_input} samples produced a BAM file. "
                "Please check the coverage distribution plots and the filtering_summary table in the 'NanoReadsFiltering' task "
                "and adjust its filtering parameters if necessary. "
                "To know more please consult the documentation: https://constellab.community/bricks/gws_variant_detection/latest/doc/variant-detection/38d49814-d0d2-41c5-95f6-4f79ea083b6d"
            )

        return {
            'output_path': Folder(result_path),
            'read_length_plots': Folder(plots_path),  # contains HTML + CSV
            'distributions': rs,                      # ResourceSet of PlotlyResource (hist-only)
        }

    def run_after_task(self):
        if self.shell_proxy:
            self.shell_proxy.clean_working_dir()
