#!/usr/bin/env bash
set -euo pipefail
CONFIG=${1:?Usage: bash 01_validate_inputs.sh config.sh}
source "$CONFIG"

mkdir -p "$BASE"/{logs,metadata,trimmed,bam/sorted,bam/filtered,bam/dedup,bam/merged,bam/tmp,peaks,qc,consensus,counts,bigwig,scripts_runtime}

# Check conda but do not fail if this script is run outside an interactive shell; SLURM steps activate it again.
if [[ "${CONDA_DEFAULT_ENV:-}" != "$CONDA_ENV" ]]; then
  echo "WARN: current conda env is '${CONDA_DEFAULT_ENV:-none}', expected '$CONDA_ENV'. Batch jobs will activate it."
fi

required=(cutadapt bowtie2 samtools picard macs2 bedtools awk sort python)
if [[ "$COUNT_METHOD" == "nucleoatac_get_count" ]]; then
  if ! conda run -n nucleoatac_py2 pyatac counts -h >/dev/null 2>&1; then
    echo "ERROR: pyatac counts not available in nucleoatac_py2 env" >&2
    exit 1
  fi
fi
for x in "${required[@]}"; do
  if ! command -v "$x" >/dev/null 2>&1; then
    if [[ "$x" == "get_count" && "$ALLOW_BEDTOOLS_COUNT_FALLBACK" == "1" ]]; then
      echo "WARN: get_count not found; fallback counting is enabled."
    else
      echo "ERROR: required command not found in PATH: $x" >&2
      exit 1
    fi
  fi
done

for f in "$MANIFEST" "$REFSEQ_ANNOTATION" "$HG19_BLACKLIST"; do
  [[ -s "$f" ]] || { echo "ERROR: missing or empty input file: $f" >&2; exit 1; }
done

# Validate Bowtie2 index basename.
if ls "${HG19_BOWTIE2_INDEX}"*.bt2 >/dev/null 2>&1; then
  echo "Found Bowtie2 small index: ${HG19_BOWTIE2_INDEX}*.bt2"
elif ls "${HG19_BOWTIE2_INDEX}"*.bt2l >/dev/null 2>&1; then
  echo "Found Bowtie2 large index: ${HG19_BOWTIE2_INDEX}*.bt2l"
else
  echo "ERROR: Bowtie2 index files not found for basename: $HG19_BOWTIE2_INDEX" >&2
  exit 1
fi

python - <<'PY' "$MANIFEST" "$BASE/metadata/libraries.tsv" "$BASE/metadata/biosamples.tsv"
import csv, os, sys
inp, libs_out, bios_out = sys.argv[1:4]
with open(inp, newline='') as fh:
    sample = fh.read(4096); fh.seek(0)
    dialect = csv.Sniffer().sniff(sample, delimiters='\t,')
    reader = csv.DictReader(fh, dialect=dialect)
    rows = list(reader)
if not rows:
    raise SystemExit('manifest has no rows')
cols = {c.lower(): c for c in rows[0].keys()}
for req in ['sample_id','fastq1','fastq2']:
    if req not in cols:
        raise SystemExit(f'manifest missing required column: {req}')
clean=[]
for r in rows:
    sid = r[cols['sample_id']].strip()
    f1 = r[cols['fastq1']].strip(); f2 = r[cols['fastq2']].strip()
    if not sid or not f1 or not f2:
        raise SystemExit(f'blank sample_id/fastq1/fastq2 in row: {r}')
    if not os.path.exists(f1):
        raise SystemExit(f'FASTQ1 missing for {sid}: {f1}')
    if not os.path.exists(f2):
        raise SystemExit(f'FASTQ2 missing for {sid}: {f2}')
    bid = r.get(cols.get('biosample_id',''), '').strip() if 'biosample_id' in cols else ''
    if not bid: bid = sid
    donor = r.get(cols.get('donor',''), '').strip() if 'donor' in cols else 'NA'
    cell = r.get(cols.get('cell_type',''), '').strip() if 'cell_type' in cols else 'NA'
    cond = r.get(cols.get('condition',''), '').strip() if 'condition' in cols else 'NA'
    trep = r.get(cols.get('technical_replicate',''), '').strip() if 'technical_replicate' in cols else 'NA'
    clean.append((sid,bid,f1,f2,donor,cell,cond,trep))
with open(libs_out,'w',newline='') as out:
    w=csv.writer(out, delimiter='\t')
    w.writerow(['sample_id','biosample_id','fastq1','fastq2','donor','cell_type','condition','technical_replicate'])
    w.writerows(clean)
seen={}
for sid,bid,f1,f2,donor,cell,cond,trep in clean:
    seen.setdefault(bid, [bid,donor,cell,cond,0])
    seen[bid][4]+=1
with open(bios_out,'w',newline='') as out:
    w=csv.writer(out, delimiter='\t')
    w.writerow(['biosample_id','donor','cell_type','condition','n_libraries'])
    for k in sorted(seen): w.writerow(seen[k])
print(f'Validated {len(clean)} libraries and {len(seen)} biological samples')
PY

cat > "$BASE/metadata/run_parameters.tsv" <<PARAMS
parameter	value
BASE	$BASE
MANIFEST	$MANIFEST
HG19_BOWTIE2_INDEX	$HG19_BOWTIE2_INDEX
REFSEQ_ANNOTATION	$REFSEQ_ANNOTATION
HG19_BLACKLIST	$HG19_BLACKLIST
CUTADAPT_MINLEN	$CUTADAPT_MINLEN
CUTADAPT_OVERLAP	$CUTADAPT_OVERLAP
BOWTIE2_MAX_INSERT	$BOWTIE2_MAX_INSERT
MAPQ_MIN	$MAPQ_MIN
SAMTOOLS_EXCLUDE_FLAGS	$SAMTOOLS_EXCLUDE_FLAGS
SAMTOOLS_REQUIRE_FLAGS	$SAMTOOLS_REQUIRE_FLAGS
MACS2_GENOME	$MACS2_GENOME
MAX_PEAK_WIDTH	$MAX_PEAK_WIDTH
CONSENSUS_MIN_SUPPORT	$CONSENSUS_MIN_SUPPORT
MIN_FILTERED_READS_PER_BIOSAMPLE	$MIN_FILTERED_READS_PER_BIOSAMPLE
MAX_MITO_FRAC	$MAX_MITO_FRAC
MAX_BLACKLIST_FRAC	$MAX_BLACKLIST_FRAC
MIN_TSS_ENRICHMENT	$MIN_TSS_ENRICHMENT
COUNT_METHOD	$COUNT_METHOD
PARAMS

echo "Validation complete. Normalized manifest: $BASE/metadata/libraries.tsv"
