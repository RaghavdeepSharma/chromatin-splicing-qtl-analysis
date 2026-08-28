# Methods overview

## 1. ATAC-seq preprocessing

Paired-end ATAC-seq libraries were validated through a sample manifest, adapter-trimmed with cutadapt, aligned to hg19 with Bowtie2, filtered with samtools, deduplicated with Picard, and assessed using read-depth, mitochondrial-read, blacklist-overlap, and TSS-enrichment metrics. Technical replicates were merged at the biological-sample level before peak calling.

MACS2 was used for peak calling. A consensus accessible-region set was constructed by merging overlapping peaks and requiring support from at least two samples. Fragment counts over the consensus regions were generated for downstream analysis.

## 2. ATAC-seq validation

Consensus-region counts were summarized with featureCounts/CPM transformations. Differential accessibility between stimulated and resting conditions was modeled with edgeR/limma-voom while accounting for donor and, where available, sample-quality covariates. Reprocessed results were compared against published Calderon differential-accessibility estimates.

## 3. Motif and genomic annotation

Vierstra motif-archetype coordinates were integrated with GENCODE comprehensive annotations using GenomicRanges/GenomicFeatures. Motif sites were annotated for gene overlap, gene biotype, genic region, strand/orientation, and distance to the nearest relevant TSS. hg19 ATAC peak coordinates were lifted to hg38 before integration with hg38 motif/QTL resources.

## 4. Cell-type-specific chromatin accessibility

Significant ATAC peaks were aggregated by cell type. Accessible peak sets were intersected with annotated transcription-factor binding sites to derive motif × cell-type accessibility profiles and accessibility-breadth metrics.

## 5. sQTL processing and integration

Immune-cell LeafCutter datasets were selected from the eQTL Catalogue. SuSiE credible-set results were organized and filtered to define high-confidence (high-PIP) and comparison sQTL sets. Datasets were mapped to compatible ATAC cell types and restricted to protein-coding/lncRNA genes before overlap with accessible motif sites.

## 6. Enrichment and regression

Motif/sQTL overlap counts were summarized in contingency tables and assessed with Fisher's exact tests. Variant-gene-level logistic regression modeled high-PIP status as a function of motif overlap while adjusting for splice-junction distance, gene-TSS distance, and local GC content. Variance inflation factors and convergence/separation diagnostics were used to assess model stability.

Firth penalized logistic regression was used for sparse or separated motif × cell-type combinations. Multiple-testing correction was performed using Benjamini-Hochberg FDR.

## 7. Dataset-pooling audit

A later audit identified cases where biologically distinct source datasets had been pooled within a shared cell-type label. Source-preserving inventories and corrected pools were constructed, followed by re-estimation of affected Firth models and comparison with the original pooled analysis.

## HPC execution

The original analysis ran on the Harvard Medical School O2 cluster using Bash, SLURM job arrays, R, and Python. This repository contains the computational workflow while excluding account-specific O2 paths and large generated outputs.
