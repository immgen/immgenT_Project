#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-04:30:00
#SBATCH --mem 6G
#SBATCH -c 1
#SBATCH -o wrapper.log
#SBATCH -e wrapper.err
#SBATCH --mail-type=END,FAIL
# Run as sbatch Ranalysis_wrapper.sh /n/groups/cbdm_lab/immgen_t .
# from Broad Run directory

### Load Relevant Modules
module load bcl2fastq/2.20.0.422
module load fastqc/0.11.5
module load cellranger/8.0.1
module load gcc/9.2.0 R/4.1.2 hdf5/1.10.1 boost/1.62.0 openmpi/3.1.0 fftw/3.3.7 java/jdk-1.8u112 geos/3.10.2
export R_LIBS="/n/groups/cbdm-db/ly82/scRNAseq_DB/R-4.1.2/library"
export R_LIBS_USER="/n/groups/cbdm-db/ly82/scRNAseq_DB/R-4.1.2/library"


### Change to target directory
cd $2

### Make sure to move and modify these files to the target directory
# Information for various antibody-derived tags, only used in adt_qc.R !!!
adt_hash_seqs=adt_hash_seq_file_2022_panel.csv
#sample names per Hashtag
samples=samples.csv
#Location of gene markers used for cell annotation (only for identifying T cells and cleaning out non-Ts)
genelist=/n/groups/cbdm_lab/immgen_t/LineageSpecGns072018_top27.csv
# IGT Number
IGT=IGT.txt
#File containing cutoffs used for T cell filtering (usually standard, may change for specific datasets)
T_cutoffs=T_cutoffs.csv
# File containing hashtag associated with Control Spleen
spleen=spleen.txt

echo "starting seurat analysis"
#mkdir -p GEX
#mkdir -p FBC

#R Script for empty droplet removal. The fourth argument in this command should be either 'Automatic' or 'Manual'
    # 'Automatic' = Cell ranger for empty droplet filtering
    # 'Manual' = manual thresholds for empty droplet filtering
    
    #Cutoffs for droplet removal if done manually
#RNA_Cutoff=200
#ADT_Cutoff=725

#Rscript $1/filter2.R combined_sample/outs/raw_feature_bc_matrix combined_sample/outs/filtered_feature_bc_matrix $2 $IGT Automatic $RNA_Cutoff $ADT_Cutoff > filter.log 2> filter.err
 
 #R Script for creating a Seurat Object and Demultiplexing.
    #This most current R script does not filter out spleen from the demultiplexing and cleaning pipeline
    # Spleen will be filtered before the removal of non-T cells and the new umap generation
Rscript $1/make_seurat_rna_adt_hto_noProjwSpleen.R GEX/filtered_feature_bc_matrix FBC/filtered_feature_bc_matrix $2 $spleen  > make_seurat_rna_adt_hto.log 2> make_seurat_rna_adt_hto.err


## TCR Analysis
    # using all singlets (pre-filtering) to generate 4 fasta files:
        # cells with paired and productive data only: 1 alpha fasta and 1 beta fasta and seurat object
        # all cells including those with nonpaired and nonproductive data: 1 alpha fasta and 1 beta fasta
Rscript $1/tcr_info_all.R $IGT seuratobject_withHTOADT_singlet.Rds > tcr_info.log 2> tcr_info.err
Rscript $1/tcr_info_nonprod_nonpaired.R $IGT seuratobject_withHTOADT_singlet.Rds > tcr_info2.log 2> tcr_info2.err

# thresholds for RNA cleaning
th_nFeature_RNA_lo=300
th_nFeature_RNA_hi=25000
th_percent_mito=5

# R script for RNA QC and cleaning
Rscript $1/rna_qc_DM.R seuratobject_withHTOADT_singlet.Rds $th_nFeature_RNA_lo $th_nFeature_RNA_hi $th_percent_mito > rna_qc.log 2> rna_qc.err

# R script for cell annotation using singler (currently, output is ignored)
Rscript $1/singler.R seuratobject_withHTOADT_singlet_postRNAfiltering.Rds > singler.log 2> singler.err

# R script for ADT QC and cleaning
Rscript $1/adt_qc.R seuratobject_withHTOADT_singlet_postRNAfiltering_withSingleR.Rds $adt_hash_seqs 500 2 > adt_qc.log 2> adt_qc.err

# R script for getting module scores for different cell types. Does no filtering in this script, but adds modulescore to metadata
Rscript $1/Tcell_filt_part1.R seuratobject_withHTOADT_singlet_postRNAfiltering_withSingleR_postADTfiltering.Rds $genelist


# R script for cleaning nonT cells based on module scores and thresholds
Rscript $1/Tcell_filt_part2wSpleen.R seuratobject_withHTOADT_singlet_postRNAfiltering_withSingleR_postADTfiltering.Rds $T_cutoffs $IGT $spleen  > T2.log 2> T2.err

# R script for making seurat smaller (removing unnecessary metadata, and HTO assay)
Rscript $1/dietSeurat.R seuratobject_withHTOADT_singlet_postRNAfiltering_withSingleR_postADTfiltering_postTfiltering.Rds seuratobject_withHTOADT_singlet_postRNAfiltering_withSingleR_postADTfiltering.Rds $2 $IGT

#naming samples 
Rscript $1/naming_samples2.R dataset.Rds dataset_clean.Rds $2 $IGT

# R script for generating final QC tables
Rscript $1/Generating_QC_table.R seuratobject_withHTOADT_singlet_preRNAfiltering_QC_stats.txt seuratobject_withHTOADT_singlet_postRNAfiltering_withSingleR_preADTfiltering.txt seuratobject_withHTOADT_singlet_postRNAfiltering_postADTfiltering_postTfiltering.csv seuratobject_IGT_singlet_postRNAfiltering_postADTfiltering_postTfiltering.csv Demultiplexing_QC_stats.txt  > Tables.log 2> Tables.err
