import csv, gzip, os, subprocess, sys
from typing import Dict, List, Tuple
import pandas as pd

# ordre souhaité des colonnes dans merged_summary.tsv
ORDER = ["Ref_Similar", "Type_1", "Type2_338", "Type2_369", "Type_3", "Other"]

# ───────────────────────────  OUTILS SHELL  ────────────────────────────
def run(cmd: List[str] | str) -> None:
    """Exécute *cmd* (liste ou chaîne). Arrête le script si retour ≠ 0."""
    if isinstance(cmd, str):
        print(cmd); subprocess.run(cmd, shell=True, check=True)
    else:
        cmd = list(map(str, cmd))
        print(" ".join(cmd)); subprocess.run(cmd, check=True)

def ensure_blast_db(ref_fa: str) -> str:
    """Crée la base BLAST+ si absente et renvoie le préfixe racine."""
    root, _ = os.path.splitext(ref_fa)
    if not os.path.exists(root + ".nsq"):
        run(["makeblastdb", "-in", ref_fa, "-dbtype", "nucl", "-out", root])
    else:
        print("[makeblastdb] database already present – skipped.")
    return root

# ───────────────────── FASTQ(.gz) ➜ FASTA (pour BLAST) ─────────────────────
def fastq_to_fasta(fq_in: str, fa_out: str) -> int:
    opener = gzip.open if fq_in.endswith(".gz") else open
    n = 0
    with opener(fq_in, "rt") as fin, open(fa_out, "w") as fout:
        while True:
            h = fin.readline()
            if not h:
                break
            s = fin.readline()
            fin.readline(); fin.readline()
            fout.write(f">{h[1:]}{s}")
            n += 1
    return n

# ─────────────────────────────  BLAST  ──────────────────────────────
def blast_outfmt6(fa: str, db_root: str, out_tsv: str, threads: int) -> None:
    run([
        "blastn", "-query", fa,
        "-db", db_root,
        "-evalue", "1000",          # aucune sélection à ce stade
        "-num_threads", threads,
        "-outfmt", "6 qseqid sseqid length slen pident",
        "-out", out_tsv,
    ])

# ─────────────────────  CLASSIFICATION PAR COUVERTURE  ─────────────────────
def classify_reads(tsv6: str, total_reads: int, min_cov: float
                   ) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Sélectionne, pour chaque read, le hit ayant la couverture la plus élevée
    (puis % identité en cas d’égalité).  Si coverage < *min_cov* => « Other ».
    """
    best: Dict[str, Dict[str, float | str]] = {}
    with open(tsv6) as fh:
        for q, s, length, slen, pid in csv.reader(fh, delimiter="\t"):
            cov = 100.0 * int(length) / int(slen)
            pid = float(pid)
            if (
                q not in best
                or cov > best[q]["cov"]
                or (cov == best[q]["cov"] and pid > best[q]["pid"])
            ):
                best[q] = {"ref": s, "cov": cov, "pid": pid}

    # applique le seuil de couverture
    for rec in best.values():
        if rec["cov"] < min_cov:
            rec["ref"] = "Other"

    # ---------- fichiers de sortie par read ----------
    metrics = pd.DataFrame.from_records(
        {
            "read_name": rid,
            "best_ref":  rec["ref"],
            "coverage_%": round(rec["cov"], 2),
            "identity_%": round(rec["pid"], 2),
        } for rid, rec in best.items()
    )

    # ---------- résumé par classe ----------
    counts: Dict[str, int] = {}
    for rec in best.values():
        counts[rec["ref"]] = counts.get(rec["ref"], 0) + 1
    # reads sans aucun hit BLAST
    no_hit = total_reads - len(metrics)
    if no_hit:
        counts["Other"] = counts.get("Other", 0) + no_hit

    summary = (
        pd.Series(counts, name="Count")
          .to_frame()
          .assign(Frequency=lambda d: 100 * d.Count / total_reads)
          .reset_index(names="Class")
    )
    return metrics, summary

# ───────────────────────  MERGE DES RÉSUMÉS  ────────────────────────
def merge_summaries(out_dir: str) -> None:
    files = [f for f in os.listdir(out_dir) if f.endswith("_summary.csv")]
    if not files:
        print("No *_summary.csv produced – nothing to merge."); return

    rows = []
    for fn in files:
        sample = fn.replace("_summary.csv", "")
        df = pd.read_csv(os.path.join(out_dir, fn))
        total = int(df["Count"].sum())
        row = {"Sample": sample, "Total_Reads": total}
        for _, r in df.iterrows():
            cls = r["Class"]
            row[f"{cls}_Count"] = int(r["Count"])
            row[f"{cls}_Frequency"] = float(r["Frequency"])
        rows.append(row)

    merged = pd.DataFrame(rows).fillna(0)

    cols = ["Sample", "Total_Reads"]
    for cls in ORDER:
        cols += [f"{cls}_Count", f"{cls}_Frequency"]
        if f"{cls}_Count" not in merged.columns:
            merged[f"{cls}_Count"] = 0
            merged[f"{cls}_Frequency"] = 0.0

    merged = merged[cols].sort_values("Sample")
    merged.to_csv(os.path.join(out_dir, "merged_summary.tsv"),
                  sep="\t", index=False)
    print("[INFO] merged_summary.tsv written.")

def main() -> None:
    fastq_dir   = sys.argv[1]
    reference_fa= sys.argv[2]
    out_dir     = sys.argv[3]
    threads      = int(sys.argv[4])
    min_cov      = int(sys.argv[5])

    os.makedirs(out_dir, exist_ok=True)
    db_root = ensure_blast_db(reference_fa)

    fastqs = sorted(
        f for f in os.listdir(fastq_dir)
        if f.endswith((".gz", ".fastq", ".fq"))
    )
    if not fastqs:
        sys.exit(f"No FASTQ files found in {fastq_dir}")

    for fq in fastqs:
        fq_path = os.path.join(fastq_dir, fq)
        sample  = fq.rsplit(".", 1)[0].replace(".fastq", "")
        fa_tmp     = os.path.join(out_dir, f"{sample}.fa")
        blast_tsv  = os.path.join(out_dir, f"{sample}.blast.tsv")
        metrics_tsv= os.path.join(out_dir, f"{sample}_read_metrics.tsv")
        summary_csv= os.path.join(out_dir, f"{sample}_summary.csv")

        total = fastq_to_fasta(fq_path, fa_tmp)
        blast_outfmt6(fa_tmp, db_root, blast_tsv, threads)

        metrics, summary = classify_reads(blast_tsv, total, min_cov)
        metrics.to_csv(metrics_tsv, sep="\t", index=False)
        summary.to_csv(summary_csv, index=False)

        os.remove(fa_tmp); os.remove(blast_tsv)

    merge_summaries(out_dir)

if __name__ == "__main__":
    main()
