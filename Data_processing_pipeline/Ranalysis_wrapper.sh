#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-04:30:00
#SBATCH --mem 10G
#SBATCH -c 4
#SBATCH -o wrapper.log
#SBATCH -e wrapper.err
#SBATCH --mail-type=END,FAIL
#
# Single, combined R analysis wrapper (RNA+ADT+HTO)
#
# Run as: sbatch Ranalysis_wrapper.sh <path_to_scripts> <run_dir>
# (run from the directory that has GEX/FBC/TCR cellranger outputs, or that
#  will be cd'ed into - see $2 below)

### Load Relevant Modules  ------------------------------------------------
module load gcc/14.2.0 R/4.4.2 hdf5/1.14.5 boost/1.87.0 openmpi/4.1.8 fftw/3.3.10 java/jdk-23.0.1 conda/miniforge3/24.11.3-0
conda activate /n/groups/cbdm_lab/odc180/Python/conda/R_4.4.2_clean_env

### Change to target directory  -------------------------------------------
cd $2

### ============================================================
### EDIT THESE FOR YOUR RUN
### ============================================================

# Does this dataset have an antibody-derived tag (ADT) library? (YES/NO)
ADT=YES
# Does this dataset have a hashtag oligo (HTO) library? (YES/NO)
# ADT and HTO are independent - HTO=YES/ADT=NO and ADT=YES/HTO=NO are both
# valid, as is ADT=NO/HTO=NO (plain RNA-only).
HTO=YES
# Did you run cellranger vdj for this dataset? (YES/NO)
TCR=YES

# Empty-droplet removal method: "Automatic" (DropletUtils knee/inflection) or "Manual"
Method=Automatic
# Only used if Method=Manual
RNA_Cutoff=300
FBC_Cutoff=500

# Information for antibody-derived tags - only used by adt_qc.R, so only
# needed if ADT=YES
adt_hash_seqs=adt_hash_seq_file_2022_panel.csv
# sample names per hashtag (or one row mapping "no_hashing" -> your sample
# name if HTO=NO) - see README for the expected format
samples=samples.csv
# Lineage marker gene list used to identify/keep T cells and drop non-T cells
genelist=/n/groups/cbdm_lab/immgen_t/LineageSpecGns072018_top27.csv
# One-line file with this run's dataset/EXP ID
EXP=EXP.txt
# Cutoffs used for T cell filtering (module-score thresholds + ADT count cutoff)
T_cutoffs=T_cutoffs.csv

# Thresholds for RNA QC (rna_qc.R)
th_nFeature_RNA_lo=300
th_nFeature_RNA_hi=25000
th_percent_mito=5

# Thresholds for ADT QC (adt_qc.R) - only used if ADT=YES
th_nCount_ADT_lo=500
n_isotype_ctrl_signal_to_flag=2

# Where the R scripts live (usually same as $1)
SCRIPTS=$1

### ============================================================
### PIPELINE - shouldn't need to edit below this line
### ============================================================

echo "starting seurat analysis (ADT=$ADT, HTO=$HTO, TCR=$TCR)"

mkdir -p GEX
mkdir -p FBC

## 2a. Empty droplet removal (writes GEX/filtered_feature_bc_matrix, and
##     FBC/filtered_feature_bc_matrix if ADT and/or HTO are present)
Rscript $SCRIPTS/filter_empty_droplets.R combined_sample/outs/raw_feature_bc_matrix combined_sample/outs/filtered_feature_bc_matrix $2 $EXP $Method $RNA_Cutoff $FBC_Cutoff $ADT $HTO > filter.log 2> filter.err

## 2b. Build the Seurat object (RNA + optional ADT assay + optional HTO
##     demultiplexing)
Rscript $SCRIPTS/build_seurat_object.R GEX/filtered_feature_bc_matrix FBC/filtered_feature_bc_matrix $2 $ADT $HTO > make_seurat.log 2> make_seurat.err

## 2c. TCR analysis (only if this dataset has a VDJ/TCR library). We only
##     use tcr_annotate_all_chains.R - it keeps every chain (productive AND
##     nonproductive, paired AND unpaired), which is a superset of what
##     tcr_annotate_productive_paired.R produces, so there's no need to run
##     both.
if [ "$TCR" == "YES" ]; then
  Rscript $SCRIPTS/tcr_annotate_all_chains.R $EXP seuratobject_singlet.Rds > tcr_info.log 2> tcr_info.err
fi

## 2d. RNA QC and cleaning (always runs)
Rscript $SCRIPTS/rna_qc.R seuratobject_singlet.Rds $th_nFeature_RNA_lo $th_nFeature_RNA_hi $th_percent_mito > rna_qc.log 2> rna_qc.err

## 2e. ADT QC and cleaning (only if ADT=YES). If skipped, carry the RNA-QC'd
##     file forward under the name the rest of the pipeline expects.
if [ "$ADT" == "YES" ]; then
  Rscript $SCRIPTS/adt_qc.R seuratobject_singlet_postRNAfiltering.Rds $adt_hash_seqs $th_nCount_ADT_lo $n_isotype_ctrl_signal_to_flag > adt_qc.log 2> adt_qc.err
else
  cp seuratobject_singlet_postRNAfiltering.Rds seuratobject_singlet_postRNAfiltering_postADTfiltering.Rds
fi

## 2f. RNA-based lineage module scoring (T/MNP/B/ILC); adds one extra
##     diagnostic plot if ADT=YES
Rscript $SCRIPTS/Tcell_filter_part1.R seuratobject_singlet_postRNAfiltering_postADTfiltering.Rds $genelist $ADT > T1.log 2> T1.err

## 2g. Drop non-T cells; RNA and ADT clustering of the surviving T cells if
##     ADT=YES. No spleen handling - if you need to exclude a control
##     hashtag, subset() on HTO_classification.simplified yourself after
##     this step.
Rscript $SCRIPTS/Tcell_filter_part2.R seuratobject_singlet_postRNAfiltering_postADTfiltering.Rds $T_cutoffs $EXP $ADT > T2.log 2> T2.err

## 2h. Slim the Seurat object down for distribution (auto-detects which
##     assays/reductions actually exist)
Rscript $SCRIPTS/slim_seurat_object.R seuratobject_singlet_postRNAfiltering_postADTfiltering_postTfiltering.Rds seuratobject_singlet_postRNAfiltering_postADTfiltering.Rds $2 $EXP

## 2i. Add human-readable sample names
Rscript $SCRIPTS/add_sample_names.R dataset.Rds dataset_clean.Rds $2 $EXP

## 2j. Final QC tables
preADTfiltering_file=seuratobject_singlet_postRNAfiltering_preADTfiltering.txt
if [ "$ADT" != "YES" ]; then
  preADTfiltering_file=NA
fi
Rscript $SCRIPTS/Generating_QC_table.R seuratobject_singlet_preRNAfiltering_QC_stats.txt $preADTfiltering_file seuratobject_singlet_postRNAfiltering_postADTfiltering_postTfiltering.csv seuratobject_EXP_singlet_postRNAfiltering_postADTfiltering_postTfiltering.csv $ADT > Tables.log 2> Tables.err

echo "Done."
