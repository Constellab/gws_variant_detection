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
    TaskInputs, TaskOutputs, task_decorator, Table, TableImporter,
    PlotlyResource, ResourceSet
)
from ..base_env.Filtering_env import FilteringShellProxyHelper


@task_decorator(
    "ReadsQualityFiltering",
    human_name="NanoReadsFiltering",
    short_description=(
        "Phred filtering → alignment (map-ont) → identity & coverage thresholds. "
        "Outputs staged HTML+CSV and three Plotly ResourceSets (RAW / Phred / ID-COV) with the same stats panel."
    )
)
class ReadsQualityFiltering(Task):
    """
    Steps:
    1) Phred filtering with `chopper -q <phred>` on each FASTQ.gz (quality trimming by score).
    2) Per-read alignment against the reference using minimap2/mappy (preset: map-ont).
    3) Keep reads that satisfy BOTH:
        - Identity (%) = (mlen / blen) × 100, where mlen is the number of matching bases
            on the aligned segment and blen is the alignment block length.
        - Coverage (%) = ((q_en − q_st) / read_length) × 100, where q_st/q_en are the
            read-aligned coordinates and read_length is the full read length (read-centric).
    4) Per-sample outputs:
        - Filtered FASTQ (after identity+coverage thresholds) and intermediate files.
        - HTML histograms (with stats panel) and CSV value dumps:
            • RAW (before Phred): length & coverage
            • After Phred:        length & coverage
            • After ID/COV:       length & coverage
    5) The Task also builds three Plotly ResourceSets (RAW / after Phred / after ID/COV),
        each containing interactive figures that include the same on-graph statistics panel.
    """

    # ---------- I/O ----------
    input_specs = InputSpecs({
        'input_folder': InputSpec(Folder, human_name="Folder containing FASTQ(.gz) files"),
        'reference_fasta': InputSpec(File, human_name="Combined reference sequence (FASTA)")
    })

    output_specs = OutputSpecs({
        'final_filtered_reads': OutputSpec(
            Folder,
            human_name="Filtered reads (identity + coverage)",
            short_description="Folder with all reads passing quality, identity, and coverage filters"
        ),
        'number_reads_filtered': OutputSpec(
            Table,
            human_name="Filtering summary table",
            short_description="Summary of reads kept per sample"
        ),

        # HTML folder outputs (as you requested)
        'plots_raw': OutputSpec(
            Folder,
            human_name="Plots — RAW input",
            short_description="Per-sample histograms before Phred filtering (HTML + CSV)"
        ),
        'plots_after_phred': OutputSpec(
            Folder,
            human_name="Plots — after Phred",
            short_description="Per-sample histograms after Phred filtering (HTML + CSV)"
        ),
        'plots_after_identity_coverage': OutputSpec(
            Folder,
            human_name="Plots — after identity/coverage",
            short_description="Per-sample histograms after identity/coverage (HTML + CSV)"
        ),

        # NEW: three separate ResourceSets (one per stage)
        'plots_raw_rs': OutputSpec(ResourceSet, human_name="RAW distributions (Plotly)"),
        'plots_phred_rs': OutputSpec(ResourceSet, human_name="After Phred distributions (Plotly)"),
        'plots_ic_rs': OutputSpec(ResourceSet, human_name="After identity/coverage distributions (Plotly)"),
    })

    # ---------- Config ----------
    config_specs: ConfigSpecs = {
        "num_threads": IntParam(default_value=15, min_value=2, short_description="Number of threads"),
        "qualified_quality_phred": IntParam(default_value=10, min_value=4, short_description="Phred quality threshold"),
        "min_identity": IntParam(default_value=60, min_value=10, short_description="Minimum identity (%)"),
        "min_coverage": IntParam(default_value=60, min_value=10, short_description="Minimum coverage (%)")
    }

    python_file_path = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "_filtering_chopper.py"
    )

    shell_proxy: ShellProxy = None

    # ---------- Helpers: Plotly with same stats panel ----------
    @staticmethod
    def _stats_text(series: pd.Series) -> str:
        s = series.dropna().astype(float)
        if s.empty:
            idx = ["count","mean","std","min","max","1%","2.5%","5%","10%","25%","50%","75%","90%","95%","97.5%","99%","missing"]
            vals = [0]*len(idx)
            return "<br>".join([f"{k}: {v}" for k, v in zip(idx, vals)])
        extra = [0.01, 0.025, 0.05, 0.10, 0.90, 0.95, 0.975, 0.99]
        desc = s.describe(percentiles=[0.25, 0.50, 0.75] + extra)
        missing = int(len(s.index) - s.notna().sum())
        def _fmt(v):
            try:
                return f"{float(v):.3f}"
            except Exception:
                return str(v)
        lines = [f"{idx}: {_fmt(desc[idx])}" for idx in
                 ["count","mean","std","min","25%","50%","75%","max"]]
        extra_pairs = [("1%", 0.01), ("2.5%", 0.025), ("5%", 0.05), ("10%", 0.10),
                       ("90%", 0.90), ("95%", 0.95), ("97.5%", 0.975), ("99%", 0.99)]
        for name, p in extra_pairs:
            val = s.quantile(p)
            lines.append(f"{name}: {_fmt(val)}")
        lines.append(f"missing: {missing}")
        return "<br>".join(lines)

    @staticmethod
    def _bins(n: int) -> int:
        return min(max(int(n * 0.10), 50), 5000)

    @staticmethod
    def _hist_with_stats(values: np.ndarray, title: str, xlab: str) -> go.Figure:
        s = pd.Series(values, dtype="float64")
        n = s.notna().sum()
        nbins = ReadsQualityFiltering._bins(int(n)) if n else 50

        fig = go.Figure()
        fig.add_trace(go.Histogram(
            x=s, nbinsx=nbins, histnorm="probability density",
            opacity=0.75, name="Histogram"
        ))

        stats_txt = ReadsQualityFiltering._stats_text(s)
        fig.update_layout(
            title=f"{title} — N={int(n)}, bins={nbins}",
            xaxis_title=xlab,
            yaxis_title="Density",
            template="plotly_white",
            legend=dict(x=1.02, y=1, xanchor="left"),
            annotations=[dict(
                x=1.42, y=1, xref="paper", yref="paper", xanchor="left",
                text=stats_txt, showarrow=False, align="left", font=dict(size=10)
            )],
            margin=dict(r=620, l=60, t=60, b=60)
        )
        return fig

    @staticmethod
    def _make_plotly_from_csv(csv_path: str, sample: str, metric: str, title_suffix: str, xlab: str) -> PlotlyResource:
        d = pd.read_csv(csv_path)
        col = "ReadLength" if metric == "Length" else "CoveragePct"
        vals = d[col].dropna().to_numpy()
        fig = ReadsQualityFiltering._hist_with_stats(vals, f"{sample} — {title_suffix}", xlab)
        pr = PlotlyResource(fig)
        # Resource names per RS: just sample + metric (no stage needed inside a per-stage RS)
        pr.name = f"{sample} — {metric}"
        return pr

    @staticmethod
    def _stage_files(plots_base: str, stage_key: str) -> list[tuple[str, str, str]]:
        """
        Return list of (sample, metric, csv_path) for a given stage folder.
        stage_key in {"RAW","PHRED","IDCOV"}.
        """
        if stage_key == "RAW":
            len_glob = os.path.join(plots_base, "*__raw_input_lengths.csv")
            cov_glob = os.path.join(plots_base, "*__raw_input_coverage.csv")
        elif stage_key == "PHRED":
            len_glob = os.path.join(plots_base, "*__after_phred_lengths.csv")
            cov_glob = os.path.join(plots_base, "*__after_phred_coverage.csv")
        else:  # IDCOV
            len_glob = os.path.join(plots_base, "*__after_identity_coverage_lengths.csv")
            cov_glob = os.path.join(plots_base, "*__after_identity_coverage_coverage.csv")

        out = []
        for fp in sorted(glob.glob(len_glob)):
            base = os.path.basename(fp)
            sample = base.replace("__raw_input_lengths.csv","")\
                         .replace("__after_phred_lengths.csv","")\
                         .replace("__after_identity_coverage_lengths.csv","")
            out.append((sample, "Length", fp))
        for fp in sorted(glob.glob(cov_glob)):
            base = os.path.basename(fp)
            sample = base.replace("__raw_input_coverage.csv","")\
                         .replace("__after_phred_coverage.csv","")\
                         .replace("__after_identity_coverage_coverage.csv","")
            out.append((sample, "Coverage", fp))
        return out

    # ---------- run ----------
    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        input_folder: Folder = inputs['input_folder']
        reference_fasta: File = inputs['reference_fasta']
        num_threads = params["num_threads"]
        qualified_quality_phred = params["qualified_quality_phred"]
        min_identity = params["min_identity"]
        min_coverage = params["min_coverage"]

        shell_proxy: ShellProxy = FilteringShellProxyHelper.create_proxy(self.message_dispatcher)

        quality_output_folder = os.path.join(shell_proxy.working_dir, "quality_output_folder")
        identity_coverage_output_folder = os.path.join(shell_proxy.working_dir, "identity_coverage_output_folder")
        final_output_folder = os.path.join(shell_proxy.working_dir, "final_output_folder")
        summary_output_folder = os.path.join(shell_proxy.working_dir, "summary_output_folder")

        plots_raw_folder = os.path.join(shell_proxy.working_dir, "plots_raw")
        plots_after_phred_folder = os.path.join(shell_proxy.working_dir, "plots_after_phred")
        plots_after_ic_folder = os.path.join(shell_proxy.working_dir, "plots_after_identity_coverage")

        for p in [
            quality_output_folder, identity_coverage_output_folder, final_output_folder,
            summary_output_folder, plots_raw_folder, plots_after_phred_folder, plots_after_ic_folder
        ]:
            os.makedirs(p, exist_ok=True)

        cmd = (
            f"python3 {self.python_file_path} "
            f"{input_folder.path} {quality_output_folder} {identity_coverage_output_folder} "
            f"{final_output_folder} {summary_output_folder} "
            f"{int(num_threads)} {int(qualified_quality_phred)} "
            f"{float(min_identity)} {float(min_coverage)} {reference_fasta.path} "
            f"{plots_raw_folder} {plots_after_phred_folder} {plots_after_ic_folder}"
        )
        shell_proxy.run(cmd, shell_mode=True)
        self.shell_proxy = shell_proxy

        # Import summary
        summary_table_path = os.path.join(summary_output_folder, "filtering_summary.csv")
        summary_table = TableImporter.call(
            File(summary_table_path),
            {"delimiter": ",", "header": 0, "file_format": "csv"}
        )

        # Build three separate ResourceSets (one per stage)
        rs_raw   = ResourceSet()
        rs_phred = ResourceSet()
        rs_ic    = ResourceSet()

        # RAW
        for sample, metric, csv_path in self._stage_files(plots_raw_folder, "RAW"):
            if metric == "Length":
                pr = self._make_plotly_from_csv(csv_path, sample, "Length",
                                                "Read length (RAW)", "Read length (bp)")
            else:
                pr = self._make_plotly_from_csv(csv_path, sample, "Coverage",
                                                "Read coverage (RAW)", "Coverage (%)")
            rs_raw.add_resource(pr, pr.name)

        # After Phred
        for sample, metric, csv_path in self._stage_files(plots_after_phred_folder, "PHRED"):
            if metric == "Length":
                pr = self._make_plotly_from_csv(csv_path, sample, "Length",
                                                "Read length after Phred", "Read length (bp)")
            else:
                pr = self._make_plotly_from_csv(csv_path, sample, "Coverage",
                                                "Coverage after Phred", "Coverage (%)")
            rs_phred.add_resource(pr, pr.name)

        # After identity/coverage
        for sample, metric, csv_path in self._stage_files(plots_after_ic_folder, "IDCOV"):
            if metric == "Length":
                pr = self._make_plotly_from_csv(csv_path, sample, "Length",
                                                "Read length after ID/COV", "Read length (bp)")
            else:
                pr = self._make_plotly_from_csv(csv_path, sample, "Coverage",
                                                "Coverage after ID/COV", "Coverage (%)")
            rs_ic.add_resource(pr, pr.name)

        return {
            'final_filtered_reads': Folder(final_output_folder),
            'number_reads_filtered': summary_table,

            'plots_raw': Folder(plots_raw_folder),
            'plots_after_phred': Folder(plots_after_phred_folder),
            'plots_after_identity_coverage': Folder(plots_after_ic_folder),

            'plots_raw_rs': rs_raw,
            'plots_phred_rs': rs_phred,
            'plots_ic_rs': rs_ic,
        }

    def run_after_task(self):
        if self.shell_proxy:
            self.shell_proxy.clean_working_dir()
