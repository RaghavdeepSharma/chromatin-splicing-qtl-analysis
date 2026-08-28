#!/usr/bin/env python3
import pandas as pd
from pathlib import Path
import argparse

parser = argparse.ArgumentParser(description="Clean featureCounts output and derive CPM/log2CPM matrices.")
parser.add_argument("--input", required=True, help="featureCounts output table")
parser.add_argument("--out-dir", required=True, help="output directory")
args = parser.parse_args()

infile = Path(args.input)
out_dir = Path(args.out_dir)
out_dir.mkdir(parents=True, exist_ok=True)
out_counts = out_dir / "consensus_counts.clean.tsv"
out_cpm = out_dir / "consensus_cpm.clean.tsv"
out_log2cpm = out_dir / "consensus_log2cpm.clean.tsv"

df = pd.read_csv(infile, sep="\t", comment="#")

# FeatureCounts annotation columns
anno_cols = ["Geneid", "Chr", "Start", "End", "Strand", "Length"]
count_cols = [c for c in df.columns if c not in anno_cols]

# Clean sample names from BAM paths
clean_names = [Path(c).name.replace(".merged.dedup.bam", "") for c in count_cols]

counts = df[["Geneid"] + count_cols].copy()
counts.columns = ["peak_id"] + clean_names
counts.to_csv(out_counts, sep="\t", index=False)

# CPM
mat = counts.set_index("peak_id")
lib_sizes = mat.sum(axis=0)
cpm = mat.div(lib_sizes, axis=1) * 1_000_000
cpm.to_csv(out_cpm, sep="\t")

# log2(CPM + 0.5), matching Calderon-style pseudocount idea
log2cpm = (cpm + 0.5).applymap(lambda x: __import__("math").log2(x))
log2cpm.to_csv(out_log2cpm, sep="\t")

print("Wrote:")
print(out_counts)
print(out_cpm)
print(out_log2cpm)
print("\nCounts matrix shape:", counts.shape)
print("CPM matrix shape:", cpm.shape)
print("Library size summary:")
print(lib_sizes.describe())
