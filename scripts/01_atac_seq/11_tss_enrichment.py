#!/usr/bin/env python3
"""Approximate Calderon-style TSS enrichment using RefSeq TSSs.
Counts fragment centers in +/-100 bp around RefSeq TSS and compares to centers ~2 kb away.
Requires pysam. TSS enrichment is calculated directly from fragment centers for reproducible QC.
"""
import argparse, gzip, os, re, subprocess, sys, tempfile
from collections import defaultdict

try:
    import pysam
except Exception as e:
    sys.exit("ERROR: pysam is required for TSS enrichment helper: %s" % e)


def open_maybe_gz(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def parse_refseq(path):
    tsses=[]
    with open_maybe_gz(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith('#'):
                continue
            f=line.rstrip('\n').split('\t')
            # UCSC refGene: bin name chrom strand txStart txEnd ...
            chrom=strand=txs=txe=None
            if len(f) >= 6 and f[2].startswith('chr') and f[3] in ['+','-']:
                chrom=f[2]; strand=f[3]; txs=int(f[4]); txe=int(f[5])
            # BED6-like: chrom start end name score strand
            elif len(f) >= 6 and f[0].startswith('chr') and f[5] in ['+','-']:
                chrom=f[0]; strand=f[5]; txs=int(f[1]); txe=int(f[2])
            # BED3-like fallback, assume plus strand
            elif len(f) >= 3 and f[0].startswith('chr'):
                chrom=f[0]; strand='+'; txs=int(f[1]); txe=int(f[2])
            else:
                continue
            tss = txs if strand == '+' else txe
            if chrom and chrom not in ['chrM','MT','M']:
                tsses.append((chrom, max(0, tss), strand))
    # Unique positions to reduce duplicate transcript inflation.
    return sorted(set(tsses))


def count_centers(bam_path, intervals):
    bam=pysam.AlignmentFile(bam_path, 'rb')
    total=0
    for chrom,start,end in intervals:
        try:
            for read in bam.fetch(chrom, max(0,start), max(0,end)):
                if read.is_unmapped or read.mate_is_unmapped or not read.is_proper_pair:
                    continue
                if read.is_secondary or read.is_supplementary or read.is_duplicate:
                    continue
                if read.mapping_quality < 30:
                    continue
                if not read.is_read1:
                    continue
                # template_length can be negative depending on orientation
                frag_start = min(read.reference_start, read.next_reference_start)
                frag_len = abs(read.template_length)
                if frag_len <= 0:
                    continue
                center = frag_start + frag_len // 2
                if start <= center < end:
                    total += 1
        except ValueError:
            # contig not in BAM
            continue
    bam.close()
    return total


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--refseq', required=True)
    ap.add_argument('--bam-dir', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--threads', default='1')
    args=ap.parse_args()
    tsses=parse_refseq(args.refseq)
    if not tsses:
        sys.exit('ERROR: no TSSs parsed from RefSeq annotation')
    center=[]; flank=[]
    for chrom,tss,strand in tsses:
        center.append((chrom, max(0,tss-100), tss+100))
        flank.append((chrom, max(0,tss-2100), max(0,tss-1900)))
        flank.append((chrom, tss+1900, tss+2100))
    bams=[os.path.join(args.bam_dir, x) for x in os.listdir(args.bam_dir) if x.endswith('.merged.dedup.bam')]
    with open(args.out,'w') as out:
        out.write('biosample_id\ttss_center_reads\ttss_flank_reads\ttss_enrichment\tn_tss\n')
        for bam in sorted(bams):
            bid=os.path.basename(bam).replace('.merged.dedup.bam','')
            c=count_centers(bam, center)
            f=count_centers(bam, flank)
            # flanks cover twice as many bases as center; normalize by interval width.
            enrich=(c/200.0)/(f/400.0) if f>0 else 0
            out.write(f'{bid}\t{c}\t{f}\t{enrich:.6g}\t{len(tsses)}\n')

if __name__ == '__main__':
    main()
