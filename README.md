# genome-pipeline

A five-step genome scan that runs end to end with **one command**:

    python3 run_pipeline.py

It analyses the phiX174 genome (NCBI NC_001422.1): GC content in sliding windows
and a six-frame search for open reading frames (ORFs).

## The steps

Each step reads the file the previous step wrote, and all intermediate files are kept.

| Step | What it does | Output |
|---|---|---|
| 1 download | Fetches the genome FASTA from NCBI (skipped if already present) | `data/genome.fasta` |
| 2 windows | Slices the genome into overlapping windows | `intermediate/windows.fasta` |
| 3 gc | Computes GC% for each window | `intermediate/window_gc.tsv` |
| 4 orfs | Finds ORFs on both strands, all six frames | `intermediate/orfs.tsv` |
| 5 summarise | Combines GC and ORFs into a report and a plot | `results/summary.md`, `results/gc_orf_map.png` |

`results/run_log.json` records the parameters used, the Python version, the
SHA-256 checksum of the input genome and the time each step took.

## Parameters

All settings are in one block at the **top of `run_pipeline.py`** (accession,
window size and step, minimum ORF length, start and stop codons, folder names).
Nothing is hard-coded further down. To analyse another genome, change `ACCESSION`.

## Reproducing from a clean checkout

    git clone <this repo>
    cd genome-pipeline
    pip install -r requirements.txt
    python3 -m unittest -v tests/test_pipeline.py    # optional: unit tests
    python3 run_pipeline.py

Needs Python 3.8+, internet access for step 1, and matplotlib (for the plot only).
On GitHub, **Actions -> Run pipeline -> Run workflow** does exactly this on a
fresh machine and commits the outputs, so the files in `data/`, `intermediate/`
and `results/` were produced from a clean checkout. The pipeline is deterministic:
the same input and parameters give identical output files (except the timestamp
and timings in `run_log.json`).

## Decisions

- Plain Python with only the standard library for the analysis, so there is little to install.
- Steps pass data through files rather than memory, so any step can be inspected or rerun alone.
- Only windows that fit fully inside the genome are used; a short last window is dropped.
- An ORF runs from the first ATG after a stop codon to the next in-frame stop and must be at least `MIN_ORF_NT` long.
- The ORF finder and GC function have unit tests in `tests/`, including reverse-strand coordinates.

## Limitations

- ORFs are only candidate genes, not confirmed ones. Many overlap, especially in a compact genome like phiX174.
- The genome is circular, but ORFs that cross the end of the sequence are not found.
- Only ATG is used as a start codon.
