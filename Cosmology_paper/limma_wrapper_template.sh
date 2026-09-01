#!/bin/bash
#SBATCH -A pi-zemmour ##SBATCH -q jfkfloor2 --exclusive 
#SBATCH --partition=caslake 
#SBATCH --nodes=1 
#SBATCH --mem=96GB
#SBATCH -J limma_%j            
#SBATCH -o limma_%j.log
#SBATCH -t 24:00:00              ##SBATCH --mem=8G
#SBATCH --mail-type=ALL         # Type of email notification- BEGIN,END,FAIL,ALL ; 
#SBATCH --mail-user=zemmour@rcc.uchicago.edu   # Email to which notifications will be sent

#run as: sbatch $SCRIPT_DIR/limma_wrapper.sh

module load python/anaconda-2022.05 
source $(conda info --base)/etc/profile.d/conda.sh
module load hdf5/1.12.0 #for BPCells but also others
module load openblas/0.3.13
conda activate /project/zemmour/david/envs/scvi_20240315
module load openblas/0.3.13

SCRIPT_DIR=/project/jfkfloor2/zemmourlab/david/immgent/immgent_integration_git
path_to_wd=/project/zemmour/david/ImmgenT/analysis/data_integration/IGT1_96/CD4

cd $path_to_wd

path_to_seurat_object=/project/zemmour/david/ImmgenT/analysis/data_integration/IGT1_96/CD4/igt1_96_CD4_20241113.Rds
output_dir=DGE_limma/20241215
so_file_name=igt1_96_CD4_20241113_sampled_so.Rds

# Rscript $SCRIPT_DIR/limma_sample_seuratobject.R $path_to_seurat_object $output_dir $so_file_name

path_to_seurat_object=$output_dir/$so_file_name
tmm_file_name=igt1_96_CD4_20241113_sampled_tmm.Rds
# Rscript $SCRIPT_DIR/limma_make_tmm.R $path_to_seurat_object $output_dir $tmm_file_name

path_to_tmm_object=$output_dir/$tmm_file_name
fit_file_name=igt1_96_CD4_20241113_sampled_fit.Rds
# Rscript $SCRIPT_DIR/limma_fit_level2.IGTHT.R $path_to_seurat_object $path_to_tmm_object $output_dir $fit_file_name

path_to_fit_object=$output_dir/$fit_file_name
prefix_file_name=in_Resting
# Rscript $SCRIPT_DIR/limma_contrasts_inResting.R $path_to_seurat_object $path_to_tmm_object $path_to_fit_object $output_dir $prefix_file_name

path_to_fit_object=$output_dir/$fit_file_name
prefix_file_name=in_Activated
Rscript $SCRIPT_DIR/limma_contrasts_inActivated.R $path_to_seurat_object $path_to_tmm_object $path_to_fit_object $output_dir $prefix_file_name



