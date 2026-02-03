#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-11:59:59
#SBATCH --mem 32G
#SBATCH -c 8
#SBATCH -o wrapper_totalvi.log
#SBATCH -e wrapper_totalvi.err
#SBATCH --mail-type=END,FAIL
# Run as sbatch run_totalvi_wrapper_OC_20250213.sh .
# from working directory

module load gcc/14.2.0 python/3.13.1 hdf5/1.14.5 boost/1.87.0 openmpi/4.1.8 fftw/3.3.10 java/jdk-23.0.1 conda/miniforge3/24.11.3-0
#module load python
#source $(conda info --base)/etc/profile.d/conda.sh

export NUMBA_CACHE_DIR="/tmp"
#source activate /project/zemmour/david/envs/scvi120_20241008
#source ~/scvi-tools_20241105/bin/activate
#source /n/groups/cbdm_lab/odc180/Python/envs/scvi_tools_py3.9_250318/bin/activate
source /n/groups/cbdm_lab/odc180/Python/envs/RHEL_env/scvi_tools_pymde_annoy_py3.9_20250522/bin/activate


SCRIPT_DIR=/n/groups/cbdm_lab/odc180/ImmgenT_workshop/Integration_Webpage/Scripts_TRBI/
## Changes based on the dataset/ run
working_dir=/n/groups/cbdm_lab/odc180/ImmgenT_workshop/Integration_Webpage/Trial_251211/
path_to_matrix=/n/groups/cbdm_lab/odc180/ImmgenT_workshop/Integration_Webpage/Trial_251211/matrix/
path_to_samples=/n/groups/cbdm_lab/odc180/ImmgenT_workshop/Integration_Webpage/Trial_251211/cell_level_samples.csv
path_to_ImmgenT=/n/groups/cbdm_lab/odc180/ImmgenT_workshop/Query_Integration/SCVI_Integration/Webpage_Trial/ImmgenT_downsampled_20250505_adata.h5ad
sample_sheet=/n/groups/cbdm_lab/odc180/ImmgenT_workshop/Integration_Webpage/Trial_251211/cell_level_samples.csv
## Stays constant
path_to_spikein=/n/groups/cbdm_lab/odc180/ImmgenT_workshop/Query_Integration/SCVI_Integration/Webpage_Trial/RNA_Integration_15k_Downsampled.h5ad
prefix=Not_classified_Full_Trial
#prefix_SCANVI=Trial_SCANVI_Combined_Fix
batchkey=IGT
categorical_covariate_keys=IGTHT
corrected_counts=False
denoised_data=False
cd $working_dir

mkdir $prefix

#python3 $SCRIPT_DIR/run_SCVI_SCANVI.py --working_dir=$working_dir --path_to_matrix=$path_to_matrix --path_to_ImmgenT=$path_to_ImmgenT --path_to_spikein=$path_to_spikein --query_sample_path=$sample_sheet --prefix=$prefix --batchkey=$batchkey --categorical_covariate_keys=$categorical_covariate_keys --corrected_counts=$corrected_counts --denoised_data=$denoised_data >totalvi.log 2>totalvi.err

path_to_anndata_not_classified=$prefix/adata_RNA.h5ad
predictions_output_file=$prefix/predictions_output_file.csv
#python3 $SCRIPT_DIR/run_SCVI_SCANVI_not_classified.py --working_dir=$working_dir --path_to_anndata=$path_to_anndata_not_classified --metadata_query=$predictions_output_file --prefix=$prefix --batchkey=$batchkey --categorical_covariate_keys=$categorical_covariate_keys --corrected_counts=$corrected_counts --denoised_data=$denoised_data >totalvi_not_classified.log 2>totalvi_not_classified.err

predictions_output_file_not_classified=$prefix/predictions_output_file_not_classified.csv
latent_level1_abT=$prefix/latent_level1_abT.csv
latent_df_scvi_abT=$prefix/latent_level2_abT.csv
latent_df_scvi_gdT=$prefix/latent_level2_gdT.csv
mde_ref_file=/n/groups/cbdm_lab/odc180/ImmgenT_workshop/ImmgenT_freeze_20250109/igt1_104_withtotalvi20250505_downsampled_level1_and_level2_mde.csv
#python3 $SCRIPT_DIR/run_mde_incremental_OC.py $working_dir $prefix $mde_ref_file $latent_df_scvi_abT $latent_df_scvi_gdT $latent_level1_abT $predictions_output_file_not_classified >mde.log 2>mde.err

conda activate /n/groups/cbdm_lab/odc180/Python/conda/R_4.4.2_clean_env

query_IGTHT=$prefix/query_all_metadata.csv
annotation_column=IGTHT
output_file=$prefix/output_file.csv
mde_plot=$SCRIPT_DIR/mde_plot.csv
Rscript $SCRIPT_DIR/Create_final_seurat_object_multipage.R $path_to_matrix $query_IGTHT $output_file $mde_plot $annotation_column $prefix >make_final_seurat_objectFinal_separate_SCVIs.log 2>make_final_seurat_objectFinal_separate_SCVIs.err
