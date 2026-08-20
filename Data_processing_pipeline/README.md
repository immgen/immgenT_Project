# scRNA-seq QC pipeline (RNA / ADT / HTO / TCR)

This pipeline is intended for the processing and QC of Single-cell 
RNA/ADT/HTO/TCR sequencing data. The pipeline assumes that you
are starting from raw BCL files, but you may start from FASTQs,
or even a count matrix if that is what you prefer.
To start from FASTQs, edit `CellRanger_Run.sh` before running.
To start from a count matrix, skip `CellRanger_Run.sh` altogether and
start directly from`'Ranalysis_wrapper.sh`.

The `CellRanger_Run.sh` and  `Ranalysis_wrapper` scripts require several 
supporting files to run,  which are heavily dependent on  the context of 
your data. Please check the description of these files below and make sure 
they are  present in your working directory whenrunning the pipeline. 
We have provided templates of each file along with the scripts.

Every script in `Ranalysis_wrapper` takes `ADT=YES/NO` and `HTO=YES/NO` 
(RNA is always assumed to be yes) and behaves in the following way for 
all four combinations:

| ADT | HTO | What you get |
|-----|-----|---------------|
| YES | YES | RNA + ADT + HTO demultiplexing (the "full" pipeline) |
| YES | NO  | RNA + ADT, no hashtag demultiplexing |
| NO  | YES | RNA + HTO demultiplexing, no antibody assay |
| NO  | NO  | RNA only |

ADT and HTO are fully independent of each other - HTO does not require ADT
to be present, and vice versa. Both are delivered through the same 10x
"Feature Barcode / Antibody Capture" library, so `filter_empty_droplets.R`
treats them as one combined pool (called "FBC" internally); the split into a
dedicated ADT assay and/or HTO assay happens afterwards, in
`build_seurat_object.R`.

**All "spleen" control-hashtag handling has been removed.** 
If you have a control sample (e.g. a spleen reference) that you want to exclude 
from downstream T cell analysis, do it yourself as a `subset()` call on 
`HTO_classification.simplified` after `Tcell_filter_part2.R` runs.

## Pipeline order

1. `CellRanger_Run.sh` - mkfastq, count, and (optionally) vdj
2. `Ranalysis_wrapper.sh` - orchestrates every R analysis step below (2a-2j).
    These are sub-steps of step 2, not separate top-level steps - they only
    ever run as part of `Ranalysis_wrapper.sh`:
    - 2a. `filter_empty_droplets.R` - empty-droplet removal
    - 2b. `build_seurat_object.R` - build the Seurat object, demultiplex HTO if present
    - 2c. `tcr_annotate_all_chains.R` - TCR, only if you ran cellranger vdj
      Resulting FASTAs are then submitted to IMGT High-V Quest for TCR annotation
    - 2d. `rna_qc.R` - RNA QC/filtering
    - 2e. `adt_qc.R` - ADT QC/filtering, only run if ADT=YES
    - 2f. `Tcell_filter_part1.R` - RNA lineage module scoring
    - 2g. `Tcell_filter_part2.R` - drop non-T cells, RNA clustering, optional ADT clustering
    - 2h. `slim_seurat_object.R` - slim the object down for distribution
    - 2i. `add_sample_names.R` - add human-readable sample names to seurat object
    - 2j. `Generating_QC_table.R` - final QC tables

`Ranalysis_wrapper.sh` runs steps 2a-2j end to end. Steps 2c and 2e are
skipped automatically (with a plain file-copy "pass-through" so every later
step still finds the file it expects) when TCR/ADT are set to NO.

3. IMGT High-V Quest annotation (manual, external - only if TCR=YES). Submit
   the FASTAs written by step 2c to IMGT, wait for it to finish, then run
   `Create_IMGT_annotation_summary_table.R` on the downloaded results. This
   step happens outside `Ranalysis_wrapper.sh` - see "IMGT High-V Quest
   annotation" below for the full walkthrough.

See `pipeline_schematic.pdf` for a visual, step-by-step diagram of the whole
pipeline (CellRanger_Run.sh through Generating_QC_table.R), including each
step's inputs and outputs.

The stage tags (`alldata`, `singlet`, `preRNAfiltering`, `postADTfiltering`,
`postTfiltering`...), they describe what's actually been done to the object 
regardless of which flags were used.

## Description of each pipeline-generated

Every script below reads/writes plain files in the run directory (`$2` in
the wrapper) - nothing is passed between steps except through these files.
Files marked "(ADT only)" or "(HTO only)" are only produced when that flag is
YES; everything else is produced regardless of flags.

**`CellRanger_Run.sh`**
- `fastqs/outs/fastq_path/` - demultiplexed FASTQs from `cellranger mkfastq`
    (if you are starting from BCL files).
- `combined_sample/outs/raw_feature_bc_matrix`,
  `combined_sample/outs/filtered_feature_bc_matrix` - `cellranger count`'s raw
  and filtered UMI count matrices (GEX rows always, "Antibody Capture" rows
  too if ADT and/or HTO=YES).
  These matrix folders have the following contents:
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
  - `TCR/outs/...` - `cellranger vdj` outputs (only if you ran that step).

**`filter_empty_droplets.R`** (empty-droplet removal)
- `GEX/filtered_feature_bc_matrix` - RNA UMI matrix for the cells kept after
  the cutoff.
- `FBC/filtered_feature_bc_matrix(2)` (ADT and/or HTO only) - the matching
  ADT/HTO "Feature Barcode" UMI matrix for the same kept cells.
- `rna_elbow_plot.png` / `fbc_elbow_plot.png` - knee/inflection diagnostic
  plots (`Method=Automatic` only).
- `FBC_RNA_cutoff_plot.png` / `RNA_cutoff_plot.png` - RNA vs. FBC UMI-count
  scatter (or RNA rank plot if no FBC data) showing where the cutoff falls.

**`build_seurat_object.R`** (build the Seurat object)
- `seuratobject_alldata.Rds` - every barcode that survived empty-droplet
  removal, with the ADT assay and/or HTO demultiplexing results attached,
  before dropping doublets/negatives.
- `seuratobject_alldata_QCstats.txt` - per-hashtag cell counts/means for the
  object above.
- `seuratobject_singlet.Rds` - `alldata` with only Singlet/Rescued cells kept
  (if HTO=NO, this is just every cell, since there's nothing to demultiplex).
- `seuratobject_singlet_QCstats.txt` - same summary as above, for the singlet
  object.
- `Demultiplexing_QC_stats.txt` - counts of Singlet/Doublet/Negative/Rescued
  cells.
- `make_seurat-1_HTO_ridgeplot.png`, `make_seurat-2_HTO_pairs.png`,
  `make_seurat-3_HTO_tsne_heatmap.png`, `hto_counts_table.csv`,
  `cells_with_no_hashtag_counts.txt` (HTO only) - hashtag demultiplexing
  diagnostics.

**`tcr_annotate_all_chains.R`** (only if TCR=YES; supplementary - doesn't
feed back into the steps below, which keep reading `seuratobject_singlet.Rds`
directly)
- `seuratobject_singlet_fullTCRinfo_productiveANDnonproductive_pairedANDunpaired.Rds`
  - singlet object plus per-cell clonotype info for every chain (productive
    AND nonproductive, paired AND unpaired).
- `<EXP>_allCells_productiveANDnonproductiveTCR_pairedANDunpaired_trv_cdr3alpha/beta.fasta`
  (+ matching `.csv`) - per-cell alpha/beta contig sequences, rows with no
  sequence for that chain dropped. Submit these to IMGT High-V Quest for TCR
  annotation (see "IMGT High-V Quest annotation" below).

**`Create_IMGT_annotation_summary_table.R`** (manual, external - only if
TCR=YES; run yourself after downloading IMGT's results, not called by
`Ranalysis_wrapper.sh` - see "IMGT High-V Quest annotation" below)
- `<EXP>_TRA.txt` / `<EXP>_TRB.txt` - IMGT's per-chain junction calls,
  reduced to the columns this pipeline uses.
- `<EXP>_summary_table.tsv` / `_xl.csv` - one row per cell: hashtag,
  sample_name, alpha/beta V/J/CDR3 calls (both contigs where present), and
  inferred `clonotype_alpha_beta*` IDs.
- `<EXP>_Enhanced_summary_table_xl.csv` - the summary table above plus each
  cell's RNA cluster and UMAP coordinates.

**`rna_qc.R`** (RNA QC/filtering)
- `seuratobject_singlet_preRNAfiltering.Rds`, `..._QC_stats.txt`,
  `..._preRNAfiltering-1.png` through `-4.png` - object and diagnostics
  before RNA filtering.
- `seuratobject_singlet_postRNAfiltering.Rds`, `..._QC_stats.txt`, `.pdf` -
  object and diagnostics after dropping low-quality/high-mitochondrial cells.

**`adt_qc.R`** (ADT QC/filtering; only runs if ADT=YES - otherwise the
wrapper just copies the RNA-QC'd file forward under the `postADTfiltering`
name so later steps see a consistent filename either way)
- `seuratobject_singlet_postRNAfiltering_preADTfiltering.Rds` / `.txt` -
  object/stats before ADT-specific filtering (low ADT counts,
  isotype/autofluorescence flags).
- `seuratobject_singlet_postRNAfiltering_postADTfiltering.Rds` / `.txt` -
  object/stats after ADT filtering.

**`Tcell_filter_part1.R`** (RNA lineage module scoring)
- `Tcell_cleanup_Step1_RNA_CellScoring.pdf` - T/MNP/B/ILC module-score
  scatter plots.
- `seuratobject_singlet_postRNAfiltering_postADTfiltering.Rds` - the same
  object, now with module scores added (overwrites the file from the
  previous step).

**`Tcell_filter_part2.R`** (drop non-T cells)
- `T_labeled_RNA_CellScoring-1.png`, `Tcell_cleanup_final_numbers.pdf` -
  cutoff diagnostics and a before/after cell-count table.
- `T_labeled_RNA_CellScoring-2.png`, `Tcell_cleanup_ADT_post_filtering.png`,
  `adt_qc_postTfilt_Dotplot_ADT_ADTclusters.png` (ADT only) - ADT-count
  cutoff and RNA/ADT backgating diagnostics.
- `seuratobject_singlet_postRNAfiltering_postADTfiltering_postTfiltering.csv`
  - per-hashtag QC summary after T-cell filtering.
- `seuratobject_EXP_singlet_postRNAfiltering_postADTfiltering_postTfiltering.csv`
  - the same summary collapsed to one row for this EXP/run.
- `seuratobject_singlet_postRNAfiltering_postADTfiltering_postTfiltering.Rds`
  - the final T-cells-only object, RNA-(and ADT-)clustered.
- `seuratobject_singlet_postRNAfiltering_postADTfiltering.Rds` - the full
  object again (all cell types, not just T cells), now labeled with
  `cell_type`/`Tcell` metadata (overwrites the file from the previous step).

**`slim_seurat_object.R`** (slim down for distribution)
- `dataset_clean.Rds` - slimmed version of the T-cells-only object.
- `dataset.Rds` - slimmed version of the full, all-cell-type object.

**`add_sample_names.R`** (add human-readable sample names)
- `dataset.Rds` / `dataset_clean.Rds` - the same two files, now with a
  `sample_name` column added (overwrites).

**`Generating_QC_table.R`** (final tables)
- `seuratobject_singlet_postRNAfiltering_postADTfiltering_postTfiltering_FinalTable.csv`
  - one row per hashtag/sample, with QC stats pulled from every stage above.
- `seuratobject_EXP_singlet_postRNAfiltering_postADTfiltering_postTfiltering_FinalTable.csv`
  - one row summarizing the whole EXP/run.

## IMGT High-V Quest annotation

If `TCR=YES`, step 2c (`tcr_annotate_all_chains.R`) writes alpha and beta
FASTAs but does not annotate them - that part happens on IMGT's website, then
gets folded back in by `Create_IMGT_annotation_summary_table.R`. This whole
step is manual and happens outside `Ranalysis_wrapper.sh`:

1. Go to [IMGT/HighV-QUEST](https://www.imgt.org/HighV-QUEST/) and
   create/log into an IMGT account - you need one before you can submit
   anything for batch annotation.
2. Submit the alpha FASTA and beta FASTA as two separate jobs, selecting the
   correct species/strain and the TRA/TRB locus for each respectively.
3. Wait for IMGT to email you that each job is done - this typically takes
   anywhere from ~6 hours to a full day depending on queue and file size
   (split very large FASTAs into smaller batches if you hit IMGT's
   per-submission sequence limit).
4. Download each job's result archive and unzip it with `tar xvf`, e.g.:
   ```bash
   mkdir -p TCR_annotations/TRA TCR_annotations/TRB
   tar xvf <alpha_job_results>.txz -C TCR_annotations/TRA
   tar xvf <beta_job_results>.txz -C TCR_annotations/TRB
   ```
   Each folder now has IMGT's numbered output files -
   `Create_IMGT_annotation_summary_table.R` only reads `6_Junction.txt` from
   each.
5. Run the summary script:
   ```bash
   Rscript Create_IMGT_annotation_summary_table.R \
     TCR_annotations/TRA TCR_annotations/TRB \
     seuratobject_singlet_fullTCRinfo_productiveANDnonproductive_pairedANDunpaired.Rds \
     seuratobject_singlet_postRNAfiltering_postADTfiltering_postTfiltering.Rds \
     samples.csv EXP.txt
   ```
   This produces the `<EXP>_TRA.txt`/`<EXP>_TRB.txt`, `<EXP>_summary_table*`,
   and `<EXP>_Enhanced_summary_table_xl.csv` files described above.

## Supporting Files needed for each run - context dependent

Templates for all of these now live in this folder - copy and edit them per
run/experiment rather than starting from scratch. They describe *your*
experiment, so even with a template you'll need to fill in your own paths,
hashtags, and sample names etc.

**`layoutfile.csv`** (template provided) - the standard Illumina sample
sheet used by `cellranger mkfastq --csv=`. Columns: `Lane,Sample,Index`. The
template has one row per lane per library (GEX/FBC/VDJ share the same lanes
but different 10x dual-index codes) - add/remove lanes to match your run.

**`libraries.csv`** (template provided) - tells `cellranger count` which
fastqs are GEX vs. FBC. Columns: `fastqs,sample,library_type`. The template
points at a GEX row and an Antibody Capture (FBC) row from the same fastq
folder; update the `fastqs` path and `sample` names for your run, and omit
the `Antibody Capture` row entirely if ADT=NO and HTO=NO.

**`feature_refs_totalseq_2023.csv`** (template provided) - the
`--feature-ref` CSV for `cellranger count`. Columns:
`id,name,read,pattern,sequence,feature_type`. Only include the rows for the
ADT antibodies and/or HTO hashtags you actually used in this run (delete the
rest, or keep the whole panel - extra unused rows are harmless, they just
won't have any reads). Only needed if ADT=YES or HTO=YES.

**`samples.csv`** (template provided) - maps each hashtag to a
human-readable sample name. Used by `build_seurat_object.R` (to decide which
hashtag rows are real, if HTO=YES) and `add_sample_names.R` (to label the
final object). Two columns, header row required (column names themselves
don't matter - the scripts read by position); the template looks like:
```
hashtag,sample_name
HT1,Control.Spleen
HT2,B16.ACT.Tumor
HT3,B16.ACT.dLN
HT4,B16.ACT/CBI.Tumor
HT5,B16.ACT/CBI.dLN
```
If HTO=NO, this pipeline still labels every cell's
`HTO_classification.simplified` as the placeholder `no_hashing`, so give it
one row so `add_sample_names.R` can still assign a sample name:
```
hashtag,sample_name
no_hashing,MySampleName
```

**`adt_hash_seq_file_2022_panel.csv`** (template provided) - antibody
metadata (clone, isotype, barcode sequence, etc.) that `adt_qc.R` joins onto
the ADT assay. Only needed if ADT=YES; edit it down to the antibodies in
your panel (column `Protein_Symbol` must match your ADT assay's feature
names).

**`EXP.txt`** (template provided) - one line, the dataset/run ID used to
prefix cell barcodes (e.g. `IGT35`). One per sequencing run.

**`T_cutoffs.csv`** (template provided) - module-score cutoffs used to call
a cell "T" vs. B/MNP/ILC. Four rows, no header:
```
cutoff_T,0.005
cutoff_B,0.05
cutoff_ILC,0.12
cutoff_MNP,0.05
```
There is no separate ADT-count cutoff here - low-`nCount_ADT` cells are
already dropped upstream in `adt_qc.R` (`th_nCount_ADT_lo`), so
`Tcell_filter_part2.R` doesn't filter on ADT counts a second time.

## Static reference files (reused across runs, already provided)

- **`LineageSpecGns072018_top27.csv`** - the lineage marker gene list used by
  `Tcell_filter_part1.R`/`Tcell_filter_part2.R` to score T/MNP/B/ILC
  identity. Reuse as-is for mouse datasets; point `genelist=` at wherever you
  keep this file.
- **CellRanger transcriptome reference** (`--transcriptome=`) and **VDJ
  reference** (`--reference=` for `cellranger vdj`) - these are external
  10x-format references, not part of this repo. See the use-case walkthrough
  for where to get public equivalents (e.g. 10x's mm10/GRCm39 references and
  IMGT-based VDJ references) if you do not have the reference paths that are
  hardcoded in `CellRanger_Run.sh`.

## Portability notes (cluster-specific paths)

This pipeline was originally written for the HMS O2 cluster, and a couple of
hardcoded paths remain in the **unchanged** scripts:

- `rna_qc.R` calls
  `source("/path/to/file/seurat_rna_umap_clustering.R")`.
  Off that cluster, this will fail. `seurat_rna_umap_clustering.R` in this
  folder is a drop-in copy of that function - either point the `source()`
  line at it, or place a copy at that exact path.
- `CellRanger_Run.sh` hardcodes a `--transcriptome=` and `--reference=`
  (VDJ) path from the original cluster - swap these for wherever you keep
  your own Cell Ranger references (see `use_case_GSE311344.md` for public
  alternatives).
- The wrapper's default `genelist=/path/to/file/LineageSpecGns072018_top27.csv`
  path is also cluster-specific - point it at wherever you keep
  `LineageSpecGns072018_top27.csv` (a copy is in this folder).

## Flags (set at the top of `Ranalysis_wrapper.sh`)

- `ADT=YES/NO`
- `HTO=YES/NO`
- `TCR=YES/NO` - whether you ran `cellranger vdj` for this sample
- `Method=Automatic/Manual` - empty-droplet cutoff method (`filter_empty_droplets.R`)
