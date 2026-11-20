#!/usr/bin/env python3
"""Genome scan pipeline: one command, five chained steps.

    python3 run_pipeline.py

Step 1  download   NCBI genome FASTA            -> data/genome.fasta
Step 2  windows    slice genome into windows    -> intermediate/windows.fasta
Step 3  gc         GC% of each window           -> intermediate/window_gc.tsv
Step 4  orfs       six-frame ORF search         -> intermediate/orfs.tsv
Step 5  summarise  report + plot                -> results/summary.md, results/gc_orf_map.png

Each step reads the file written by the step before it, and every intermediate
file is kept so any step can be inspected on its own.
"""

# =====================================================================
# PARAMETERS - change values here only; nothing below is hard-coded.
# =====================================================================
ACCESSION = "NC_001422.1"          # NCBI nucleotide record to analyse (phiX174 genome)
DOWNLOAD_URL = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
                "?db=nuccore&id={acc}&rettype=fasta")
WINDOW_SIZE = 200                  # bases per window
WINDOW_STEP = 100                  # bases between window starts
MIN_ORF_NT = 300                   # shortest ORF to report, in nucleotides (incl. stop codon)
START_CODON = "ATG"
STOP_CODONS = ("TAA", "TAG", "TGA")
TOP_N_ORFS = 10                    # ORFs listed in the summary
DATA_DIR = "data"
INTERMEDIATE_DIR = "intermediate"
RESULTS_DIR = "results"
# =====================================================================

import argparse
import hashlib
import json
import platform
import sys
import time
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


# ---------- small reusable functions (tested in tests/) ----------

def read_fasta(path):
    """Return list of (id, sequence) for every record in a FASTA file."""
    records, rec_id, chunks = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if rec_id is not None:
                    records.append((rec_id, "".join(chunks).upper()))
                rec_id, chunks = line[1:].split()[0], []
            elif rec_id is not None:
                chunks.append(line)
    if rec_id is not None:
        records.append((rec_id, "".join(chunks).upper()))
    return records


def reverse_complement(seq):
    return seq.translate(COMPLEMENT)[::-1]


def gc_percent(seq):
    """GC% over A/C/G/T bases only; None if there are none."""
    acgt = sum(seq.count(b) for b in "ACGT")
    return None if acgt == 0 else 100 * (seq.count("G") + seq.count("C")) / acgt


def find_orfs(seq, min_len=MIN_ORF_NT, start=START_CODON, stops=STOP_CODONS):
    """Six-frame ORF search. Returns dicts with 1-based forward-strand coordinates.

    An ORF runs from the first start codon after a stop (or sequence start) to
    the next in-frame stop codon, and must be at least min_len nt long including
    the stop. Circularity is ignored: ORFs spanning the end are not found.
    """
    n = len(seq)
    orfs = []
    for strand, s in (("+", seq), ("-", reverse_complement(seq))):
        for frame in range(3):
            orf_start = None
            for i in range(frame, n - 2, 3):
                codon = s[i:i + 3]
                if orf_start is None:
                    if codon == start:
                        orf_start = i
                elif codon in stops:
                    length = i + 3 - orf_start
                    if length >= min_len:
                        if strand == "+":
                            a, b = orf_start + 1, i + 3
                        else:
                            a, b = n - (i + 3) + 1, n - orf_start
                        orfs.append({"strand": strand, "frame": frame + 1, "start": a,
                                     "end": b, "length_nt": length, "length_aa": length // 3 - 1,
                                     "gc_percent": gc_percent(s[orf_start:i + 3])})
                    orf_start = None
    orfs.sort(key=lambda o: (-o["length_nt"], o["start"]))
    for k, o in enumerate(orfs, 1):
        o["orf_id"] = f"ORF{k:03d}"
    return orfs


# ---------- pipeline steps ----------

def step1_download(refresh=False):
    out = Path(DATA_DIR) / "genome.fasta"
    if out.exists() and not refresh:
        print(f"[1/5] download: using existing {out}")
        return out
    url = DOWNLOAD_URL.format(acc=ACCESSION)
    print(f"[1/5] download: {url}")
    out.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(url, timeout=60) as resp:
        data = resp.read()
    if not data.startswith(b">"):
        sys.exit("Download did not return a FASTA record; check ACCESSION / network.")
    out.write_bytes(data)
    return out


def step2_windows(genome_path):
    rec_id, seq = read_fasta(genome_path)[0]
    out = Path(INTERMEDIATE_DIR) / "windows.fasta"
    out.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with open(out, "w") as fh:
        for start in range(0, len(seq) - WINDOW_SIZE + 1, WINDOW_STEP):
            fh.write(f">{rec_id}:{start + 1}-{start + WINDOW_SIZE}\n{seq[start:start + WINDOW_SIZE]}\n")
            count += 1
    print(f"[2/5] windows: {count} windows of {WINDOW_SIZE} bp (step {WINDOW_STEP}) -> {out}")
    return out


def step3_gc(windows_path):
    out = Path(INTERMEDIATE_DIR) / "window_gc.tsv"
    with open(out, "w") as fh:
        fh.write("window_id\tstart\tend\tgc_percent\n")
        for wid, seq in read_fasta(windows_path):
            span = wid.split(":")[1]
            s, e = span.split("-")
            gc = gc_percent(seq)
            fh.write(f"{wid}\t{s}\t{e}\t{'NA' if gc is None else f'{gc:.2f}'}\n")
    print(f"[3/5] gc: window GC% -> {out}")
    return out


def step4_orfs(genome_path):
    _, seq = read_fasta(genome_path)[0]
    orfs = find_orfs(seq)
    out = Path(INTERMEDIATE_DIR) / "orfs.tsv"
    cols = ["orf_id", "strand", "frame", "start", "end", "length_nt", "length_aa", "gc_percent"]
    with open(out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for o in orfs:
            gc = "NA" if o["gc_percent"] is None else f"{o['gc_percent']:.2f}"
            fh.write("\t".join(str(o[c]) if c != "gc_percent" else gc for c in cols) + "\n")
    print(f"[4/5] orfs: {len(orfs)} ORFs >= {MIN_ORF_NT} nt -> {out}")
    return out


def read_tsv(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t"))) for line in fh if line.strip()]


def step5_summarise(genome_path, gc_path, orfs_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rec_id, seq = read_fasta(genome_path)[0]
    gc_rows = [r for r in read_tsv(gc_path) if r["gc_percent"] != "NA"]
    orfs = read_tsv(orfs_path)
    Path(RESULTS_DIR).mkdir(exist_ok=True)

    genome_gc = gc_percent(seq)
    gcs = [float(r["gc_percent"]) for r in gc_rows]
    lo = min(gc_rows, key=lambda r: float(r["gc_percent"]))
    hi = max(gc_rows, key=lambda r: float(r["gc_percent"]))
    n_plus = sum(o["strand"] == "+" for o in orfs)

    lines = [f"# Genome scan summary: {rec_id}", "",
             f"- Genome length: {len(seq)} bp; overall GC: {genome_gc:.2f}%",
             f"- Windows: {len(gc_rows)} (size {WINDOW_SIZE}, step {WINDOW_STEP})",
             f"- Window GC range: {min(gcs):.2f}% (positions {lo['start']}-{lo['end']}) "
             f"to {max(gcs):.2f}% (positions {hi['start']}-{hi['end']})",
             f"- ORFs >= {MIN_ORF_NT} nt: {len(orfs)} ({n_plus} forward, {len(orfs) - n_plus} reverse)",
             "", f"## Longest {min(TOP_N_ORFS, len(orfs))} ORFs", "",
             "| ORF | strand | start | end | nt | aa | GC% |", "|---|---|---|---|---|---|---|"]
    for o in orfs[:TOP_N_ORFS]:
        lines.append(f"| {o['orf_id']} | {o['strand']} | {o['start']} | {o['end']} | "
                     f"{o['length_nt']} | {o['length_aa']} | {o['gc_percent']} |")
    lines += ["", "ORFs are open reading frames (ATG to stop), not confirmed genes. "
              "Overlapping and nested ORFs are expected."]
    (Path(RESULTS_DIR) / "summary.md").write_text("\n".join(lines) + "\n")

    mids = [(int(r["start"]) + int(r["end"])) / 2 for r in gc_rows]
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6), sharex=True,
                                   gridspec_kw={"height_ratios": [2, 1]})
    ax1.plot(mids, gcs, lw=1)
    ax1.axhline(genome_gc, ls="--", c="grey", lw=0.8)
    ax1.set_ylabel("GC% per window")
    ax1.set_title(f"{rec_id}: GC content and ORFs >= {MIN_ORF_NT} nt")
    for o in orfs:
        y = 1 if o["strand"] == "+" else 0
        ax2.broken_barh([(int(o["start"]), int(o["end"]) - int(o["start"]))], (y + 0.1, 0.8),
                        alpha=0.4, color="tab:blue" if y else "tab:orange")
    ax2.set_yticks([0.5, 1.5]); ax2.set_yticklabels(["- strand", "+ strand"])
    ax2.set_xlabel("Position (bp)")
    fig.tight_layout(); fig.savefig(Path(RESULTS_DIR) / "gc_orf_map.png", dpi=150); plt.close(fig)
    print(f"[5/5] summarise: {RESULTS_DIR}/summary.md and gc_orf_map.png")


def write_run_log(genome_path, timings):
    params = {k: v for k, v in globals().items() if k.isupper() and not k.startswith("_")}
    params["STOP_CODONS"] = list(STOP_CODONS)
    log = {
        "run_utc": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S"),
        "python": platform.python_version(),
        "input_file": str(genome_path),
        "input_sha256": hashlib.sha256(Path(genome_path).read_bytes()).hexdigest(),
        "parameters": params,
        "step_seconds": timings,
    }
    (Path(RESULTS_DIR) / "run_log.json").write_text(json.dumps(log, indent=2) + "\n")


def main():
    p = argparse.ArgumentParser(description="Run the genome scan pipeline end to end.")
    p.add_argument("--refresh-download", action="store_true", help="download the genome again")
    args = p.parse_args()

    timings = {}
    def timed(name, fn, *a, **kw):
        t = time.time(); result = fn(*a, **kw); timings[name] = round(time.time() - t, 2); return result

    genome = timed("download", step1_download, args.refresh_download)
    windows = timed("windows", step2_windows, genome)
    gc = timed("gc", step3_gc, windows)
    orfs = timed("orfs", step4_orfs, genome)
    timed("summarise", step5_summarise, genome, gc, orfs)
    write_run_log(genome, timings)
    print("Done. Results in", RESULTS_DIR)


if __name__ == "__main__":
    main()
