#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import subprocess
import gzip
import mappy as mp
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from pathlib import Path


# =========================
#     Utils: Plot/Stats
# =========================

def _descriptive_stats_from_values(values: list[float], name: str) -> pd.DataFrame:
    """Return a DataFrame of descriptive stats from a list of numeric values."""
    if not values:
        idx = ["count","mean","std","min","max","1%","2.5%","5%","10%","25%","50%","75%","90%","95%","97.5%","99%","missing"]
        return pd.DataFrame({"value":[0]*len(idx)}, index=idx)

    s = pd.Series(values, dtype="float64", name=name)
    extra = [0.01, 0.025, 0.05, 0.10, 0.90, 0.95, 0.975, 0.99]
    desc = s.describe(percentiles=[0.25, 0.50, 0.75] + extra).to_frame(name="value")
    desc.loc["missing"] = int(len(s.index) - s.notna().sum())
    return desc


def _export_values_csv(values: list[float], out_csv: Path, colname: str):
    """Save raw values to CSV so the Task can build PlotlyResource later."""
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({colname: values}).to_csv(out_csv, index=False)


def _plot_hist_html(values: list[float], out_html: Path, title: str, xlab: str):
    """
    Write a Plotly histogram + stats. Bins = max(10% of N, 50) capped at 5000.
    The stats panel is pushed far to the right to avoid any overlap with legends.
    """
    out_html.parent.mkdir(parents=True, exist_ok=True)

    if not values:
        out_html.write_text(f"<html><body><h3>{title}</h3><p>No data.</p></body></html>")
        return

    n = len(values)
    bins = int(max(n * 0.10, 50))
    bins = min(bins, 5000)

    s = pd.Series(values, dtype="float64", name=xlab)
    stats_df = _descriptive_stats_from_values(values, xlab)

    hist = go.Histogram(x=s, nbinsx=bins, histnorm="probability density",
                        opacity=0.75, name="Histogram")

    def _fmt(v):
        try:
            return f"{float(v):.3f}"
        except Exception:
            return str(v)
    stats_txt = "<br>".join([f"{idx}: {_fmt(val[0])}" for idx, val in stats_df.iterrows()])

    fig = go.Figure([hist])
    fig.update_layout(
        title=title + f" — N={n}, bins={bins}",
        xaxis_title=xlab,
        yaxis_title="Density",
        template="plotly_white",
        legend=dict(x=1.02, y=1, xanchor="left"),
        annotations=[dict(
            x=1.42, y=1, xref="paper", yref="paper",
            xanchor="left",
            text=stats_txt, showarrow=False, align="left", font=dict(size=10)
        )],
        margin=dict(r=620, l=60, t=60, b=60)
    )
    fig.write_html(str(out_html))


# =========================
#   Main Filtering Class
# =========================

class FastpFiltering:
    """
      1) Phred filter via 'chopper -q <phred>' on FASTQ.gz
      2) Align each read with minimap2/mappy (preset 'map-ont')
      3) Keep reads passing:
           - identity = mlen / blen * 100  (on aligned segment)
           - coverage = (q_en - q_st) / len(read) * 100   [read-centric]
      4) Write ALL passing reads (no randomization/no quota)
      5) Produce per-sample histograms for:
           - Read length: RAW, after Phred, after ID/COV
           - Coverage (%): RAW, after Phred, after ID/COV
         Each plot is exported as HTML **and** a CSV (values) to feed PlotlyResource.
    """

    def __init__(
        self,
        input_folder,
        quality_output_folder,
        identity_coverage_output_folder,
        final_output_folder,
        summary_output_folder,
        num_threads,
        qualified_quality_phred,
        min_identity,
        min_coverage,
        reference_genome,
        plots_raw_folder,            # per-sample RAW plots (length & coverage)
        plots_after_phred_folder,    # after chopper (length & coverage)
        plots_after_ic_folder,       # after identity/coverage (length & coverage)
    ):
        self.input_folder = input_folder
        self.quality_output_folder = quality_output_folder
        self.identity_coverage_output_folder = identity_coverage_output_folder
        self.final_output_folder = final_output_folder
        self.summary_output_folder = summary_output_folder
        self.num_threads = int(num_threads)
        self.qualified_quality_phred = int(qualified_quality_phred)
        self.min_identity = float(min_identity)
        self.min_coverage = float(min_coverage)
        self.reference_genome = reference_genome
        self.plots_raw_folder = plots_raw_folder
        self.plots_after_phred_folder = plots_after_phred_folder
        self.plots_after_ic_folder = plots_after_ic_folder

        self.aligner = mp.Aligner(reference_genome, preset="map-ont")
        if not self.aligner:
            raise RuntimeError(f"Failed to load reference: {reference_genome}")

        for p in [
            self.quality_output_folder,
            self.identity_coverage_output_folder,
            self.final_output_folder,
            self.summary_output_folder,
            self.plots_raw_folder,
            self.plots_after_phred_folder,
            self.plots_after_ic_folder,
        ]:
            os.makedirs(p, exist_ok=True)

    @staticmethod
    def _read_fastq_record(handle):
        name = handle.readline()
        if not name:
            return None
        seq = handle.readline()
        plus = handle.readline()
        qual = handle.readline()
        if not seq or not plus or not qual:
            return None
        return name.rstrip("\n"), seq.rstrip("\n"), plus.rstrip("\n"), qual.rstrip("\n")

    def _best_identity_and_coverage(self, read_sequence: str):
        """
        identity (%) = mlen / blen * 100
        coverage (%) = (q_en - q_st) / len(read) * 100
        Returns (identity, coverage) of the BEST passing hit; otherwise (None, None).
        """
        read_len = len(read_sequence)
        if read_len == 0:
            return None, None

        best_ident = None
        best_cov = None

        for hit in self.aligner.map(read_sequence):
            aligned_q_len = hit.q_en - hit.q_st
            blen = getattr(hit, "blen", 0)
            mlen = getattr(hit, "mlen", 0)
            if aligned_q_len <= 0 or blen <= 0:
                continue

            ident = (mlen / blen) * 100.0
            cov = (aligned_q_len / read_len) * 100.0

            if ident >= self.min_identity and cov >= self.min_coverage:
                if best_ident is None or ident > best_ident:
                    best_ident, best_cov = ident, cov

        return best_ident, best_cov

    def _best_coverage_any(self, read_sequence: str):
        """Best coverage (%) across all hits (ignores thresholds)."""
        read_len = len(read_sequence)
        if read_len == 0:
            return None

        best_cov = None
        for hit in self.aligner.map(read_sequence):
            aligned_q_len = hit.q_en - hit.q_st
            blen = getattr(hit, "blen", 0)
            if aligned_q_len <= 0 or blen <= 0:
                continue
            cov = (aligned_q_len / read_len) * 100.0
            if (best_cov is None) or (cov > best_cov):
                best_cov = cov
        return best_cov

    def _collect_raw_len_and_cov(self, gz_fastq_path: str):
        """One pass to get RAW lengths and coverage (best hit, no thresholds)."""
        lengths, coverages = [], []
        with gzip.open(gz_fastq_path, "rt") as fh:
            while True:
                rec = self._read_fastq_record(fh)
                if rec is None:
                    break
                _, seq, _, _ = rec
                lengths.append(len(seq))
                cov = self._best_coverage_any(seq)
                if cov is not None:
                    coverages.append(cov)
        return lengths, coverages

    def filter_reads(self):
        summary_rows = []
        fastq_files = sorted([f for f in os.listdir(self.input_folder) if f.endswith(".fastq.gz")])

        for fastq_file in fastq_files:
            sample = os.path.splitext(os.path.splitext(fastq_file)[0])[0]  # remove .fastq.gz
            input_path = os.path.join(self.input_folder, fastq_file)
            quality_output_path = os.path.join(self.quality_output_folder, fastq_file)
            ic_output_path = os.path.join(self.identity_coverage_output_folder, fastq_file)
            final_output_path = os.path.join(self.final_output_folder, fastq_file)

            # ---------- RAW ----------
            raw_lengths, raw_coverages = self._collect_raw_len_and_cov(input_path)

            _plot_hist_html(
                raw_lengths,
                Path(self.plots_raw_folder) / f"{sample}__raw_input_lengths.html",
                title=f"{sample} — Read length distribution (RAW, before chopper)",
                xlab="Read length (bp)"
            )
            _export_values_csv(raw_lengths,
                               Path(self.plots_raw_folder) / f"{sample}__raw_input_lengths.csv",
                               "ReadLength")

            _plot_hist_html(
                raw_coverages,
                Path(self.plots_raw_folder) / f"{sample}__raw_input_coverage.html",
                title=f"{sample} — Coverage distribution (RAW, before chopper)",
                xlab="Coverage (%)"
            )
            _export_values_csv(raw_coverages,
                               Path(self.plots_raw_folder) / f"{sample}__raw_input_coverage.csv",
                               "CoveragePct")

            # ---------- Phred filter (chopper) ----------
            quality_cmd = (
                f"gunzip -c {input_path} | "
                f"chopper -q {self.qualified_quality_phred} --threads {self.num_threads} | "
                f"gzip > {quality_output_path}"
            )
            result = subprocess.run(quality_cmd, shell=True)
            if result.returncode != 0:
                print(f"[ERROR] Phred filtering failed on {fastq_file}", file=sys.stderr)
                continue

            lengths_after_phred: list[int] = []
            cov_after_phred: list[float] = []
            lengths_after_ic: list[int] = []
            cov_after_ic: list[float] = []
            total_after_phred = 0
            total_after_ic = 0

            # ---------- Align + identity/coverage filters ----------
            with gzip.open(quality_output_path, "rt") as infile, \
                 gzip.open(ic_output_path, "wt") as ic_out, \
                 gzip.open(final_output_path, "wt") as final_out:

                while True:
                    rec = self._read_fastq_record(infile)
                    if rec is None:
                        break
                    name, seq, plus, qual = rec
                    read_len = len(seq)
                    total_after_phred += 1
                    lengths_after_phred.append(read_len)

                    cov_any = self._best_coverage_any(seq)
                    if cov_any is not None:
                        cov_after_phred.append(cov_any)

                    ident, cov = self._best_identity_and_coverage(seq)
                    if ident is not None:
                        for handle in (ic_out, final_out):
                            handle.write(f"{name}\n{seq}\n+\n{qual}\n")
                        total_after_ic += 1
                        lengths_after_ic.append(read_len)
                        cov_after_ic.append(cov)

            # ---------- Summary row ----------
            summary_rows.append({
                "File": fastq_file,
                "Total_Reads_RAW": len(raw_lengths),
                "Total_Reads_After_Phred": total_after_phred,
                "Total_Reads_Passing_Identity_Coverage": total_after_ic,
                "Min_Identity": self.min_identity,
                "Min_Coverage": self.min_coverage,
                "Phred_Threshold": self.qualified_quality_phred,
            })

            # ---------- Plots (HTML + CSV at each stage) ----------
            # After Phred
            _plot_hist_html(
                lengths_after_phred,
                Path(self.plots_after_phred_folder) / f"{sample}__after_phred_lengths.html",
                title=f"{sample} — Read length distribution after Phred (chopper)",
                xlab="Read length (bp)"
            )
            _export_values_csv(lengths_after_phred,
                               Path(self.plots_after_phred_folder) / f"{sample}__after_phred_lengths.csv",
                               "ReadLength")

            _plot_hist_html(
                cov_after_phred,
                Path(self.plots_after_phred_folder) / f"{sample}__after_phred_coverage.html",
                title=f"{sample} — Coverage distribution after Phred (best hit, no thresholds)",
                xlab="Coverage (%)"
            )
            _export_values_csv(cov_after_phred,
                               Path(self.plots_after_phred_folder) / f"{sample}__after_phred_coverage.csv",
                               "CoveragePct")

            # After identity/coverage (kept reads only)
            _plot_hist_html(
                lengths_after_ic,
                Path(self.plots_after_ic_folder) / f"{sample}__after_identity_coverage_lengths.html",
                title=f"{sample} — Read length distribution after identity/coverage filters",
                xlab="Read length (bp)"
            )
            _export_values_csv(lengths_after_ic,
                               Path(self.plots_after_ic_folder) / f"{sample}__after_identity_coverage_lengths.csv",
                               "ReadLength")

            _plot_hist_html(
                cov_after_ic,
                Path(self.plots_after_ic_folder) / f"{sample}__after_identity_coverage_coverage.html",
                title=f"{sample} — Coverage distribution after identity/coverage filters",
                xlab="Coverage (%)"
            )
            _export_values_csv(cov_after_ic,
                               Path(self.plots_after_ic_folder) / f"{sample}__after_identity_coverage_coverage.csv",
                               "CoveragePct")

        # ---------- Global summary CSV ----------
        summary_df = pd.DataFrame(summary_rows)
        summary_csv_path = os.path.join(self.summary_output_folder, "filtering_summary.csv")
        summary_df.to_csv(summary_csv_path, index=False)
        print(f"Summary saved: {summary_csv_path}")
        print(f"RAW plots (HTML+CSV): {self.plots_raw_folder}")
        print(f"After Phred plots (HTML+CSV): {self.plots_after_phred_folder}")
        print(f"After ID/COV plots (HTML+CSV): {self.plots_after_ic_folder}")


# =========================
#   CLI entry point
# =========================
if __name__ == "__main__":
    if len(sys.argv) != 14:
        print(
            "Usage: filtering_identity_coverage.py "
            "<input_folder> <quality_out> <identity_coverage_out> <final_out> <summary_out> "
            "<threads> <phred> <min_identity> <min_coverage> <reference.fa> "
            "<plots_raw_out> <plots_after_phred_out> <plots_after_ic_out>",
            file=sys.stderr
        )
        sys.exit(1)

    input_folder = sys.argv[1]
    quality_output_folder = sys.argv[2]
    identity_coverage_output_folder = sys.argv[3]
    final_output_folder = sys.argv[4]
    summary_output_folder = sys.argv[5]
    num_threads = int(sys.argv[6])
    qualified_quality_phred = int(sys.argv[7])
    min_identity = float(sys.argv[8])
    min_coverage = float(sys.argv[9])
    reference_genome = sys.argv[10]
    plots_raw_folder = sys.argv[11]
    plots_after_phred_folder = sys.argv[12]
    plots_after_ic_folder = sys.argv[13]

    filt = FastpFiltering(
        input_folder=input_folder,
        quality_output_folder=quality_output_folder,
        identity_coverage_output_folder=identity_coverage_output_folder,
        final_output_folder=final_output_folder,
        summary_output_folder=summary_output_folder,
        num_threads=num_threads,
        qualified_quality_phred=qualified_quality_phred,
        min_identity=min_identity,
        min_coverage=min_coverage,
        reference_genome=reference_genome,
        plots_raw_folder=plots_raw_folder,
        plots_after_phred_folder=plots_after_phred_folder,
        plots_after_ic_folder=plots_after_ic_folder,
    )
    filt.filter_reads()
