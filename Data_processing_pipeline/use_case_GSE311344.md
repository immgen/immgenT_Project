# Use case: GSE311344 (ImmgenT IGT35 - B16 melanoma, ACT + checkpoint blockade)

This walks through the  pipeline end to end on a real public dataset:
[GSE311344](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE311344), *"T
cells in B16 tumor model with ACT and checkpoint blockade - scRNA-seq,
CITE-seq and TCRab"* (part of the ImmgenT Open-Source Project, IGT35).

This dataset is a good worked example precisely because it uses all three
optional data types this pipeline supports: **GEX** (RNA), **FBC**
(CITE-seq/TotalSeq-C **ADT** + **HTO** hashtag multiplexing, combined in one
Feature Barcode library), and **VDJ/TCR** (paired TCR alpha/beta). So this
run uses `ADT=YES`, `HTO=YES`, `TCR=YES`.

The 5 biological samples (hashtag-multiplexed pools of sorted CD3+ T cells)
are:

| GSM sample prefix | Description |
|---|---|
| I35H1_spleen | Spleen, untreated |
| I35H2_B16_ACT | B16 tumor, adoptive cell transfer (ACT) only |
| I35H3_LNinguinal_B16_ACT | Draining inguinal LN, ACT only |
| I35H4_B16_ACTaPD1a41BB | B16 tumor, ACT + anti-PD-1/anti-4-1BB |
| I35H5_LNinguinal_B16_ACTaPD1a41BB | Draining inguinal LN, ACT + anti-PD-1/anti-4-1BB |

Each of the 5 samples has 3 GEO records (one per 10x library type): `_GEX`,
`_FBC`, `_TCR` - 15 GSMs total.

## 1. Download the FASTQs

GEO itself only hosts the processed matrices for this series (cell metadata,
gene/protein lists, the aggregated RNA/protein `.mtx` files). The **raw
FASTQs live in SRA**, under BioProject
[PRJNA1368984](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1368984). Use the
[SRA Run Selector](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA1368984)
to see the exact SRR run accession(s) behind each GSM (a given GSM/library
can be split across more than one SRR run if it was sequenced on multiple
flow cells/lanes).

With [sra-tools](https://github.com/ncbi/sra-tools) installed:

```bash
# repeat per SRR accession you find in the Run Selector for each GSM
prefetch SRRxxxxxxx
fasterq-dump --split-files -O fastqs_raw/ SRRxxxxxxx
```

`fasterq-dump` will not produce 10x-style filenames
(`<sample>_S1_L00X_R1_001.fastq.gz`) on its own - rename/gzip the output so
`cellranger mkfastq`/`cellranger count` can find them, e.g.:

```bash
mv SRRxxxxxxx_1.fastq I35H2_B16_ACT_GEX_S1_L001_R1_001.fastq
mv SRRxxxxxxx_2.fastq I35H2_B16_ACT_GEX_S1_L001_R2_001.fastq
gzip I35H2_B16_ACT_GEX_S1_L001_R*_001.fastq
```

Do this for all three library types for all 5 samples (GEX, FBC, TCR).

## 2. Concatenate FASTQs per library type

To create GSMs for separate samples, we split the original FASTQs into their
corresponding sample based on the downstream demultiplexing using HTODemux. 
To process these FASTQs in the same was as the original, concatenate the 
corresponding library files per read before running Cell Ranger - keep GEX, 
FBC, and TCR fastqs in three separate folders since they go through different 
Cell Ranger steps:

```bash
mkdir -p fastqs/GEX fastqs/FBC fastqs/TCR

cat I35*_GEX_S1_L00*_R1_001.fastq.gz > fastqs/GEX/I35_GEX_S1_L001_R1_001.fastq.gz
cat I35*_GEX_S1_L00*_R2_001.fastq.gz > fastqs/GEX/I35_GEX_S1_L001_R2_001.fastq.gz
# ...repeat for FBC and TCR
```

## 3. Run `CellRanger_Run.sh` for GEX + FBC + VDJ together

You need:

- **`libraries.csv`** pointing at the GEX and FBC fastq folders:
  ```
  fastqs,sample,library_type
  fastqs/GEX,I35_GEX,Gene Expression
  fastqs/FBC,I35_FBC,Antibody Capture
  ```
- **`feature_refs_totalseq_2023.csv`** - GEO ships the exact panel used for
  this study as `GSE311344_feature_reference.csv.gz` (download it and use it
  directly, or cross-check it against `feature_refs_totalseq_2023.csv` in
  this pipeline folder - it should contain both the HTO hashtag rows and the
  TotalSeq-C ADT antibody rows used for this run, all tagged
  `Antibody Capture`).
- A mouse transcriptome reference for `--transcriptome=` (e.g. 10x's
  prebuilt `refdata-gex-mm10-2020-A`, or `refdata-cellranger-arc-GRCm39` if
  running fresher references) - swap this in for the paththat's hardcoded 
  in the script.
- A VDJ reference for `--reference=` in the `cellranger vdj` step (10x
  publishes `refdata-cellranger-vdj-GRCm38-alts-ensembl` for mouse; the
  original script uses an IMGT-derived custom reference - either works, but
  it must match the species/build of your transcriptome reference).

Then:

```bash
sbatch CellRanger_Run.sh /path/to/run_dir/I35
```

This produces `combined_sample/outs/{raw,filtered}_feature_bc_matrix` (GEX +
FBC together, `combined_sample` being the `--id=` Cell Ranger was run with)
and `TCR/outs/filtered_contig_annotations.csv` +
`TCR/outs/all_contig.fasta` (and `all_contig_annotations.csv`) - exactly what
`filter_empty_droplets.R`, `build_seurat_object.R`, and
`tcr_annotate_all_chains.R` expect (give the run its own `EXP.txt`, e.g.
`IGT35`).

## 4. Run `Ranalysis_wrapper.sh` for QC and the final Seurat object

Start here if you already have the cellranger count matrices and only wish to use 
the R/ Seurat processing pipeline.

The `Ranalysis_wrapper.sh` pipeline starts from the cellranger output, which has
two count matrices containing the combined GEX and FBC (if present) libraries, 
with the following structure:

raw_feature_matrix (unfiltered count matrix):
'barcodes.tsv(.gz)' - unique identifiers for each cell:
```
ATGCCGTCCCTGGATT-1
CCGCCGTATCTGGATC-1
GGGCCGTACCGTTATG-1
...
GAACCGTACCGTTTGA-1
```
'features.tsv(.gz)' - Feature names:
If FBC present:
```
ENSMUSG00000051951	Xkr4	Gene Expression
ENSMUSG00000089699	Gm1992	Gene Expression
ENSMUSG00000102343	Gm37381	Gene Expression
ENSMUSG00000025900	Rp1	Gene Expression
ENSMUSG00000025902	Sox17	Gene Expression
ENSMUSG00000104328	Gm37323	Gene Expression
...
HT1	HT1	Antibody Capture
HT2	HT2	Antibody Capture
HT3	HT3	Antibody Capture
...
TCRVG1.1	TCRVG1.1	Antibody Capture
CD28	CD28	Antibody Capture
CD38	CD38	Antibody Capture
CD16_CD32	CD16_CD32	Antibody Capture
```

If FBC not present:
```
ENSMUSG00000051951	Xkr4	Gene Expression
ENSMUSG00000089699	Gm1992	Gene Expression
ENSMUSG00000102343	Gm37381	Gene Expression
ENSMUSG00000025900	Rp1	Gene Expression
ENSMUSG00000025902	Sox17	Gene Expression
ENSMUSG00000104328	Gm37323	Gene Expression
...
ENSMUSG00000094855	AC133095.1	Gene Expression
ENSMUSG00000095019	AC234645.1	Gene Expression
ENSMUSG00000095041	AC149090.1	Gene Expression
```


'matrix.mtx(.gz)' - count matrix in condensed mtx format
number of rows in 'matrix.mtx(.gz)' = number of features in 'features.tsv(.gz)'
number of columns in 'matrix.mtx(.gz)' = number of barcodes in 'barcodes.tsv(.gz)'

In the run directory (which now has `GEX/`, `FBC/`, `TCR/`,
`EXP.txt`, `samples.csv` mapping each of that pool's hashtags to a mouse ID,
`T_cutoffs.csv`, `feature_refs`/`adt_hash_seq_file` matching your panel, and
`LineageSpecGns072018_top27.csv`):

```bash
# at the top of Ranalysis_wrapper.sh:
#   ADT=YES
#   HTO=YES
#   TCR=YES

sbatch Ranalysis_wrapper.sh /path/to/scripts /path/to/run_dir/I35
```

This runs empty-droplet removal -> Seurat object + HTO demultiplexing -> TCR
annotation -> RNA QC -> ADT QC -> T cell lineage
scoring and non-T cell removal -> slim/final Seurat object
(`dataset.Rds`/`dataset_clean.Rds`) -> QC summary tables. See
`README.md` for exactly what each of those steps needs and produces.

## 5. TCR clonotype annotation via IMGT HighV-QUEST

`tcr_annotate_all_chains.R` (run automatically by the wrapper, step 2c above)
writes, per sample:

- `<EXP>_allCells_productiveANDnonproductiveTCR_pairedANDunpaired_trv_cdr3alpha.fasta`
- `<EXP>_allCells_productiveANDnonproductiveTCR_pairedANDunpaired_trv_cdr3beta.fasta`

These are per-cell alpha and beta chain nucleotide FASTAs pulled straight
from `cellranger vdj`'s `all_contig.fasta` (rows with no sequence for that
chain are dropped), with each sequence's header carrying `<cell ID> V-Gene:
<v_gene> CDR3:<cdr3>`.

To get IMGT's own V/D/J/CDR3 annotation and germline gene calls, then fold
that back onto your Seurat object:

1. Go to [IMGT/HighV-QUEST](https://www.imgt.org/HighV-QUEST/).
2. Create/log into an IMGT account - this is required before you can submit
   anything for batch annotation.
3. Submit the alpha FASTA and beta FASTA as separate jobs, selecting
   *Mus musculus (Strain C57BL/6J)* and the TRA/TRB locus for each respectively.
4. Wait for IMGT to finish and email you - this is typically anywhere from
   ~6 hours to a full day depending on queue and file size (split very large
   FASTAs into smaller batches if you hit IMGT's per-submission sequence
   limit).
5. Download each job's result archive and unzip it with `tar xvf`, e.g.:
   ```bash
   mkdir -p TCR_annotations/TRA TCR_annotations/TRB
   tar xvf <alpha_job_results>.txz -C TCR_annotations/TRA
   tar xvf <beta_job_results>.txz -C TCR_annotations/TRB
   ```
   Each folder now contains IMGT's numbered output files (`1_Summary.txt`,
   `6_Junction.txt`, etc.) - `6_Junction.txt` is the one the next step reads.
6. Run `Create_IMGT_annotation_summary_table.R` to join the IMGT annotation
   back onto your Seurat object's metadata by cell ID (the same
   `<EXP>.<barcode>` identifiers used throughout this pipeline). See
   `README.md` for its exact arguments and outputs.

## Notes specific to this dataset

- This series pools multiple mice per 10x lane via hashtag, so `HTO=YES` is
  required to demultiplex `I35H1`...`I35H5` back into per-mouse samples -
  check `GSE311344_cell_metadata.tsv.gz` on GEO for the ground-truth
  hashtag-to-mouse mapping if you want to sanity-check your own
  `build_seurat_object.R` demultiplexing output against it.
