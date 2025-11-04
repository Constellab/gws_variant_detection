#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import subprocess
import random
import pysam

# --- REQUIRED PLOTTING IMPORTS ---
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from sklearn.neighbors import KernelDensity
from sklearn.model_selection import GridSearchCV, KFold
# ---------------------------------

def _ensure_dir(p):
    os.makedirs(p, exist_ok=True)

def _descriptive_stats(s: pd.Series) -> pd.DataFrame:
    """Return a describe() table with extra percentiles + 'missing' row."""
    extra = [0.01, 0.025, 0.05, 0.10, 0.90, 0.95, 0.975, 0.99]
    desc = s.describe(percentiles=[0.25, 0.50, 0.75] + extra).to_frame(name="value")
    desc.loc["missing"] = int(len(s.index) - s.notna().sum())
    return desc

def _fit_multimodal_kde_1d(s: pd.Series):
    """Cross-validated Gaussian KDE; returns (xs, ys, best_bandwidth)."""
    x = s.values.reshape(-1, 1)
    n = len(s)
    std = s.std(ddof=1)
    base_bw = max(1e-6, s.mean() * 0.01) if (std <= 0 or n <= 1) else 1.06 * std * (n ** (-1/5))
    bw_grid = np.unique(np.clip(base_bw * np.array([0.15, 0.25, 0.35, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0]), 1e-6, None))
    cv = KFold(n_splits=5 if n >= 10 else max(2, min(5, n)), shuffle=True, random_state=42)
    grid = GridSearchCV(KernelDensity(kernel="gaussian"), param_grid={"bandwidth": bw_grid}, cv=cv)
    grid.fit(x)
    best_kde = grid.best_estimator_
    bw = float(grid.best_params_["bandwidth"])
    xs = np.linspace(s.min(), s.max(), 2048).reshape(-1, 1)
    ys = np.exp(best_kde.score_samples(xs))
    return xs.ravel(), ys, bw

def _plot_series_html(s: pd.Series, out_html_path, title, xaxis_title, bins=500):
    """Histogram + KDE + stats annotation (HTML output)."""
    _ensure_dir(os.path.dirname(out_html_path))
    if s.empty:
        with open(out_html_path, "w") as fh:
            fh.write(f"<html><body><h3>{title}</h3><p>No data.</p></body></html>")
        return

    stats_df = _descriptive_stats(s)
    hist = go.Histogram(x=s, nbinsx=bins, histnorm="probability density",
                        opacity=0.55, name="Histogram")
    xs, ys, bw = _fit_multimodal_kde_1d(s)
    kde_line = go.Scatter(x=xs, y=ys, mode="lines", name=f"KDE (bw={bw:.4g})")

    def _fmt(v):
        try:
            return f"{float(v):.3f}"
        except Exception:
            return str(v)

    stats_txt = "<br>".join([f"{idx}: {_fmt(val[0])}" for idx, val in stats_df.iterrows()])

    fig = go.Figure([hist, kde_line])
    fig.update_layout(
        title=title,
        xaxis_title=xaxis_title,
        yaxis_title="Density",
        template="plotly_white",
        # >>> pushed far right with extra margin so it never overlaps the legend
        annotations=[dict(
            x=1.35, y=1, xref="paper", yref="paper",
            xanchor="left",
            text=stats_txt, showarrow=False, align="left", font=dict(size=10)
        )],
        margin=dict(r=560, l=60, t=60, b=60)
    )
    fig.write_html(out_html_path)

def _plot_lengths_html(lengths, out_html_path, title="Read length distribution", bins=500):
    s = pd.Series(lengths, name="ReadLength", dtype="float64") if lengths else pd.Series(dtype="float64")
    _plot_series_html(s, out_html_path, title, "Read length (bp)", bins=bins)

def _plot_coverage_html(coverages, out_html_path, title="Read coverage distribution", bins=500):
    s = pd.Series(coverages, name="CoveragePct", dtype="float64") if coverages else pd.Series(dtype="float64")
    _plot_series_html(s, out_html_path, title, "Coverage (%)", bins=bins)

_INTERNAL_SEED = 42

class ReadsMappingMinimap:
    """
    Steps:
      1) Minimap2 mapping (map-ont), exclude secondary/supplementary.
      2) Optional removal of the shortest fraction (≤ low_tail_q) on primaries.
      3) Exact downsampling by QNAME after filtering if target_reads > 0.
      4) Per sample outputs:
           - <sample>.sorted.bam (+ index)
           - HTML plots (histogram + KDE) for Read Length & Read Coverage (final primaries)
           - CSVs: <sample>.final_lengths.csv (ReadLength), <sample>.final_coverage.csv (CoveragePct)
    """
    def __init__(self, ref_genome, fastq_dir, threads, output_dir, target_reads, plots_dir,
                 low_tail_q=0.05, inclusive=True):
        self.ref_genome   = ref_genome
        self.fastq_dir    = fastq_dir
        self.threads      = int(threads)
        self.output_dir   = output_dir
        self.target_reads = max(0, int(target_reads))
        self.plots_dir    = plots_dir
        self.low_tail_q   = float(low_tail_q)
        self.inclusive    = bool(inclusive)

    def create_output_directory(self):
        os.makedirs(self.output_dir, exist_ok=True)
        _ensure_dir(self.plots_dir)

    def run_ok(self, cmd):
        subprocess.run(cmd, check=True)

    def index_reference(self):
        if not os.path.exists(self.ref_genome + ".mmi"):
            self.run_ok(["minimap2", "-d", f"{self.ref_genome}.mmi", self.ref_genome])
        if not os.path.exists(self.ref_genome + ".fai"):
            self.run_ok(["samtools", "faidx", self.ref_genome])

    def _count_unique_primary_qnames(self, bam_path):
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            names = {
                r.query_name
                for r in bam.fetch(until_eof=True)
                if (not r.is_secondary and not r.is_supplementary and r.reference_id >= 0)
            }
        return len(names)

    def _collect_read_lengths(self, bam_path):
        lengths = []
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for r in bam.fetch(until_eof=True):
                if r.is_secondary or r.is_supplementary or r.reference_id < 0:
                    continue
                qlen = r.query_length if r.query_length is not None else r.infer_query_length(always=False)
                if qlen is None:
                    seq = r.query_sequence
                    qlen = len(seq) if seq else None
                if qlen is not None and qlen > 0:
                    lengths.append(int(qlen))
        return lengths

    def _collect_read_coverages(self, bam_path):
        """Read-centric coverage (%) = (query_alignment_length / query_length) * 100."""
        covs = []
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for r in bam.fetch(until_eof=True):
                if r.is_secondary or r.is_supplementary or r.reference_id < 0:
                    continue
                qlen = r.query_length if r.query_length is not None else r.infer_query_length(always=False)
                if qlen is None or qlen <= 0:
                    seq = r.query_sequence
                    qlen = len(seq) if seq else 0
                qaln = r.query_alignment_length or 0  # excludes soft-clipped bases
                if qlen > 0:
                    covs.append((qaln / qlen) * 100.0)
        return covs

    def _filter_by_length_percentile(self, in_bam, out_bam, q=0.05, inclusive=True):
        """Drop the shortest fraction of primary reads according to q-quantile."""
        if q <= 0:
            print(f"[Q{int(max(q,0)*100)}] q <= 0 → no length filtering performed.")
            if in_bam != out_bam:
                self.run_ok(["samtools", "view", "-b", "-o", out_bam, in_bam])
                self.run_ok(["samtools", "index", out_bam])
            kept = self._count_unique_primary_qnames(in_bam)
            return None, kept, 0

        if q >= 1:
            q = 0.999999

        lengths = self._collect_read_lengths(in_bam)
        n_total = len(lengths)
        if n_total == 0:
            print(f"[Q{int(q*100)}] {os.path.basename(in_bam)}: no reads → nothing to filter.")
            return None, 0, None

        threshold = float(np.quantile(lengths, q))

        tmp = out_bam + ".tmp.bam"
        kept = 0
        with pysam.AlignmentFile(in_bam, "rb") as bam, pysam.AlignmentFile(tmp, "wb", header=bam.header) as out:
            for r in bam.fetch(until_eof=True):
                if r.is_secondary or r.is_supplementary or r.reference_id < 0:
                    continue
                qlen = r.query_length if r.query_length is not None else r.infer_query_length(always=False)
                if qlen is None:
                    seq = r.query_sequence
                    qlen = len(seq) if seq else 0
                keep = (qlen > threshold) if inclusive else (qlen >= threshold)
                if keep:
                    out.write(r)
                    kept += 1

        self.run_ok(["samtools", "sort", "-@", str(self.threads), "-o", out_bam, tmp])
        self.run_ok(["samtools", "index", out_bam])
        os.remove(tmp)

        removed = n_total - kept
        rule_txt = "drop ≤ threshold" if inclusive else "drop < threshold"
        print(f"[Q{int(q*100)}] threshold={threshold:.3f} bp | total={n_total} → kept={kept} removed={removed} ({rule_txt})")
        return threshold, kept, removed

    def _downsample_bam_exact_primaries(self, in_bam, out_bam, k):
        random.seed(_INTERNAL_SEED)
        with pysam.AlignmentFile(in_bam, "rb") as bam:
            prim_qnames = [
                r.query_name
                for r in bam.fetch(until_eof=True)
                if (not r.is_secondary and not r.is_supplementary and r.reference_id >= 0)
            ]
        uniq = list(set(prim_qnames))
        picked = set(random.sample(uniq, k))

        tmp = out_bam + ".tmp.bam"
        with pysam.AlignmentFile(in_bam, "rb") as bam, pysam.AlignmentFile(tmp, "wb", header=bam.header) as out:
            for r in bam.fetch(until_eof=True):
                if r.is_secondary or r.is_supplementary:
                    continue
                if r.query_name in picked:
                    out.write(r)

        self.run_ok(["samtools", "sort", "-@", str(self.threads), "-o", out_bam, tmp])
        self.run_ok(["samtools", "index", out_bam])
        os.remove(tmp)
        print(f"[DOWN] kept {len(picked)} primaries → {os.path.basename(out_bam)}")

    def _plot_and_export(self, bam_path, sample_name):
        # Lengths
        read_lengths = self._collect_read_lengths(bam_path)
        out_len_html = os.path.join(self.plots_dir, f"{sample_name}.final_read_lengths.html")
        _plot_lengths_html(read_lengths, out_len_html,
                           title=f"{sample_name} — Lengths of final primary reads", bins=500)
        pd.DataFrame({"ReadLength": read_lengths}).to_csv(
            os.path.join(self.plots_dir, f"{sample_name}.final_lengths.csv"), index=False
        )

        # Coverage
        coverages = self._collect_read_coverages(bam_path)
        out_cov_html = os.path.join(self.plots_dir, f"{sample_name}.final_coverage.html")
        _plot_coverage_html(coverages, out_cov_html,
                            title=f"{sample_name} — Coverage of final primary reads", bins=500)
        pd.DataFrame({"CoveragePct": coverages}).to_csv(
            os.path.join(self.plots_dir, f"{sample_name}.final_coverage.csv"), index=False
        )

    def map_reads(self):
        file_extension = ".fastq.gz"
        platform = "ONT"
        read_files = sorted([f for f in os.listdir(self.fastq_dir) if f.endswith(file_extension)])
        if not read_files:
            raise FileNotFoundError(f"No {file_extension} files found in {self.fastq_dir}")

        rg_counter = 0
        for read_file in read_files:
            bname = os.path.splitext(read_file)[0]
            rg_counter += 1
            rgid = f"RG{rg_counter}"
            rgsn = bname
            readgroup = f"@RG\\tID:{rgid}\\tPL:{platform}\\tPU:unit1\\tLB:lib1\\tSM:{rgsn}"

            fq_path = os.path.join(self.fastq_dir, read_file)
            print(f"[MAP] {read_file}")

            # 1) minimap2 → SAM → BAM
            output_bam = os.path.join(self.output_dir, f"{bname}.bam")
            minimap_cmd = [
                "minimap2", "-t", str(self.threads),
                "-R", readgroup,
                "-ax", "map-ont", "--secondary=no",
                self.ref_genome, fq_path
            ]
            samtools_view = ["samtools", "view", "-bS", "-"]
            p1 = subprocess.Popen(minimap_cmd, stdout=subprocess.PIPE)
            p2 = subprocess.Popen(samtools_view, stdin=p1.stdout, stdout=subprocess.PIPE)
            with open(output_bam, 'wb') as f:
                for chunk in p2.stdout:
                    f.write(chunk)
            p1.stdout.close()
            p2.communicate()

            # 2) primary only
            primary_bam = os.path.join(self.output_dir, f"{bname}.primary.bam")
            subprocess.run(["samtools", "view", "-b", "-F", "2304", "-o", primary_bam, output_bam], check=True)

            # 3) sort + index primaries
            primary_sorted = os.path.join(self.output_dir, f"{bname}.primary.sorted.bam")
            subprocess.run(["samtools", "sort", "-@", str(self.threads), "-o", primary_sorted, primary_bam], check=True)
            subprocess.run(["samtools", "index", primary_sorted], check=True)

            try:
                os.remove(output_bam)
                os.remove(primary_bam)
            except Exception:
                pass

            # 4) low-tail removal
            qfilt_sorted = os.path.join(self.output_dir, f"{bname}.primary.qfilt.sorted.bam")
            threshold, kept_after_q, removed = self._filter_by_length_percentile(
                primary_sorted, qfilt_sorted, q=self.low_tail_q, inclusive=self.inclusive
            )

            try:
                os.remove(primary_sorted)
                if os.path.exists(primary_sorted + ".bai"):
                    os.remove(primary_sorted + ".bai")
            except Exception:
                pass

            if not kept_after_q:
                print(f"[SKIP] {bname}: no reads left after q{int(max(self.low_tail_q,0)*100)} filter.")
                try:
                    if os.path.exists(qfilt_sorted):
                        os.remove(qfilt_sorted)
                    if os.path.exists(qfilt_sorted + ".bai"):
                        os.remove(qfilt_sorted + ".bai")
                except Exception:
                    pass
                continue

            # 5) exact downsampling if requested
            final_sorted = os.path.join(self.output_dir, f"{bname}.sorted.bam")
            if self.target_reads > 0:
                n_prim = self._count_unique_primary_qnames(qfilt_sorted)
                if n_prim < self.target_reads:
                    print(f"[SKIP] {bname}: {n_prim} < target {self.target_reads}. Discarded.")
                    try:
                        os.remove(qfilt_sorted)
                        if os.path.exists(qfilt_sorted + ".bai"):
                            os.remove(qfilt_sorted + ".bai")
                    except Exception:
                        pass
                    continue
                self._downsample_bam_exact_primaries(qfilt_sorted, final_sorted, self.target_reads)
                try:
                    os.remove(qfilt_sorted)
                    if os.path.exists(qfilt_sorted + ".bai"):
                        os.remove(qfilt_sorted + ".bai")
                except Exception:
                    pass
            else:
                os.replace(qfilt_sorted, final_sorted)
                os.replace(qfilt_sorted + ".bai", final_sorted + ".bai")

            # 6) HTML + CSV exports for ResourceSet (length & coverage)
            self._plot_and_export(final_sorted, bname)

        print("[DONE] All samples processed.")

    def preprocess_bams(self):
        for f in os.listdir(self.output_dir):
            if f.endswith(".bam"):
                new_name = f.replace(".fastq.sorted.bam", ".sorted.bam")
                if new_name != f:
                    os.rename(os.path.join(self.output_dir, f), os.path.join(self.output_dir, new_name))
            if f.endswith(".bai"):
                new_name = f.replace(".fastq.sorted.bam.bai", ".sorted.bam.bai")
                if new_name != f:
                    os.rename(os.path.join(self.output_dir, f), os.path.join(self.output_dir, new_name))
        print("Files have been renamed to match the desired format.")

# ---- CLI ----
if __name__ == "__main__":
    if len(sys.argv) != 8:
        print("Usage: python3 _mapping_minimap2.py <fastq_dir> <ref.fa> <threads> <out_dir> <target_reads> <plots_dir> <low_tail_q>", file=sys.stderr)
        sys.exit(1)

    fastq_dir    = sys.argv[1]
    ref_genome   = sys.argv[2]
    threads      = int(sys.argv[3])
    output_dir   = sys.argv[4]
    target_reads = int(sys.argv[5])
    plots_dir    = sys.argv[6]
    low_tail_q   = float(sys.argv[7])

    mapper = ReadsMappingMinimap(ref_genome, fastq_dir, threads, output_dir, target_reads, plots_dir,
                                 low_tail_q=low_tail_q, inclusive=True)
    mapper.create_output_directory()
    mapper.index_reference()
    mapper.map_reads()
    mapper.preprocess_bams()
