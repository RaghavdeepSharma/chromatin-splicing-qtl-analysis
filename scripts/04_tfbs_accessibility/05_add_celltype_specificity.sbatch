#!/bin/bash
#SBATCH --job-name=add_spec
#SBATCH --partition=short
#SBATCH --time=08:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

set -eo pipefail

BASE="${PROJECT_ROOT:?Set PROJECT_ROOT to the project root}"
ACC_OUT=$BASE/tfbs_celltype_accessibility_20260701
BAMPE_Q05_OUT=$ACC_OUT/BAMPE_raw_processed_accessibility_q05_significant_peaks

INDIR=$BAMPE_Q05_OUT/annotated_with_accessibility
OUTDIR=$BAMPE_Q05_OUT/annotated_with_accessibility_plus_specificity

mkdir -p "$OUTDIR"

LIST=$BAMPE_Q05_OUT/tables/final_q05_files_for_specificity.list

find "$INDIR" \
  -maxdepth 1 \
  -name "*.general_annotated.celltype_accessibility.tsv.gz" \
  | grep -v '\.chr' \
  | sort > "$LIST"

echo "Files to process:"
wc -l "$LIST"

while read f; do
  base=$(basename "$f")
  out="$OUTDIR/$base"

  echo "Adding cell-type specificity metrics to: $base"

  zcat "$f" | awk '
    BEGIN {
      FS=OFS="\t"
    }

    NR==1 {
      nacc=0
      for (i=1; i<=NF; i++) {
        if ($i ~ /^accessible_/) {
          nacc++
          acc_col[nacc]=i
        }
      }

      if (nacc==0) {
        print "ERROR: no accessible_ columns found" > "/dev/stderr"
        exit 1
      }

      print $0, "n_accessible_celltypes", "accessibility_breadth_fraction"
      next
    }

    {
      n=0
      for (j=1; j<=nacc; j++) {
        v=$(acc_col[j])
        if (v==1 || v=="1" || v=="TRUE" || v=="true" || v=="True") {
          n++
        }
      }

      frac=n/nacc
      print $0, n, frac
    }
  ' | gzip -c > "$out"

  gzip -t "$out"
done < "$LIST"

echo "Done adding specificity metrics."
echo "Output directory:"
echo "$OUTDIR"
