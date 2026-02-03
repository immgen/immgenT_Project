#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-04:00:00
#SBATCH --mem 10G
#SBATCH -c 4
#SBATCH -o wrapper.log
#SBATCH -e wrapper.err
#SBATCH --mail-type=END,FAIL
# Run as sbatch Ranalysis_wrapper_RNA.R /n/groups/cbdm_lab/immgen_t .
# from Broad Run directory

### Load Relevant Modules
#module load bcl2fastq/2.19.1.403
#module load fastqc/0.11.5
#module load cellranger
module load gcc/9.2.0 R/4.1.2 hdf5/1.10.1 boost/1.62.0 openmpi/3.1.0 fftw/3.3.7 java/jdk-1.8u112 geos/3.10.2
export R_LIBS="/n/groups/cbdm-db/ly82/scRNAseq_DB/R-4.1.2/library"
export R_LIBS_USER="/n/groups/cbdm-db/ly82/scRNAseq_DB/R-4.1.2/library"


### Change to target directory
cd $2

### Make sure to move and modify these files to the target directory
# Change Layout File
#layoutfile=layoutfile.csv
# FASTQ information
#libraries=libraries.csv
# Hashtag file
#feature_refs=feature_refs_totalseq_2023.csv
# adds experimental conditions and HT descriptions; Taken from Ian's Single-Cel$
#project_description=project_description.csv
# Hashes that are not used
#hashes_to_remove=hash_adt_rows_to_remove_from_analysis.txt
# Information for various antibody-derived tags, only used in adt_qc.R !!!
adt_hash_seqs=adt_hash_seq_file_2022_panel.csv
#sample names per Hashtag
samples=samples.csv
#Location of BCL run  folder
genelist=/n/groups/cbdm_lab/immgen_t/LineageSpecGns072018_top27.csv
# IGT Number
IGT=IGT.txt
T_cutoffs=T_cutoffs.csv
spleen=spleen.txt
RNA_Cutoff=700
ADT_Cutoff=0



echo "starting seurat analysis"
mkdir -p GEX
mkdir -p FBC
Rscript $1/filter2.R combined_sample/outs/raw_feature_bc_matrix combined_sample/outs/filtered_feature_bc_matrix $2 $IGT Manual $RNA_Cutoff $ADT_Cutoff > filter.log 2> filter.err

Rscript $1/make_seurat_rna_adt_hto_noProjwSpleen_RNA_only.R GEX/filtered_feature_bc_matrix FBC/filtered_feature_bc_matrix $2 $spleen  > make_seurat_rna_adt_hto.log 2> make_seurat_rna_adt_hto.err 

# TCR Analysis
Rscript $1/tcr_info_all.R $IGT seuratobject_withHTOADT_singlet.Rds > tcr_info.log 2> tcr_info.err
Rscript $1/tcr_info_nonprod_nonpaired.R $IGT seuratobject_withHTOADT_singlet.Rds > tcr_info2.log 2> tcr_info2.err

th_nFeature_RNA_lo=500
th_nFeature_RNA_hi=10000
th_percent_mito=5


Rscript $1/rna_qc_DM.R seuratobject_withHTOADT_singlet.Rds $th_nFeature_RNA_lo $th_nFeature_RNA_hi $th_percent_mito > rna_qc.log 2> rna_qc.err

#Not used
#Rscript $1/singler.R seuratobject_withHTOADT_singlet_postRNAfiltering.Rds > singler.log 2> singler.err


Rscript $1/Tcell_filt_part1_RNA.R seuratobject_withHTOADT_singlet_postRNAfiltering.Rds $genelist > T1.log 2> T1.err

Rscript $1/naming_samples_RNA.R seuratobject_withHTOADT_singlet_postRNAfiltering_withSingleR.Rds $2


## TCR Analysis
#Rscript $1/tcr_info.R $IGT seuratobject_withHTOADT_singlet_postRNAfiltering_withSingleR.Rds > tcr_info.log 2> tcr_info.err


Rscript $1/Tcell_filt_part2wspleen_RNA_OC.R seuratobject_withHTOADT_singlet_postRNAfiltering_withSingleR.Rds $T_cutoffs $IGT $spleen  > T2.log 2> T2.err

Rscript $1/dietSeurat_RNA.r seuratobject_withHTOADT_singlet_postRNAfiltering_withSingleR_postADTfiltering_postTfiltering.Rds seuratobject_withHTOADT_singlet_postRNAfiltering_withSingleR_postADTfiltering.Rds $2 $IGT

Rscript $1/naming_samples2.R dataset.Rds dataset_clean.Rds $2 $IGT
Rscript $1/Generating_QC_table_RNA.R seuratobject_withHTOADT_singlet_preRNAfiltering_QC_stats.txt seuratobject_withHTOADT_singlet_postRNAfiltering_postADTfiltering_postTfiltering.csv seuratobject_IGT_singlet_postRNAfiltering_postADTfiltering_postTfiltering.csv Demultiplexing_QC_stats.txt
