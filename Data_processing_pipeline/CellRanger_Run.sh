#!/bin/bash
#SBATCH -p medium
#SBATCH -t 23:59:59
#SBATCH --mem 64G
#SBATCH -c 4
#SBATCH -o wrapper.log
#SBATCH -e wrapper.err
#SBATCH --mail-type=END,FAIL


# Run as sbatch /n/groups/cbdm_lab/immgen_t/CellRanger_Run.sh .


### Load Relevant Modules
module load bcl2fastq/2.20.0.422
module load fastqc/0.11.5
module load cellranger/8.0.1
module load gcc/9.2.0 R/4.1.2 hdf5/1.10.1 boost/1.62.0 openmpi/3.1.0 fftw/3.3.7 java/jdk-1.8u112 geos/3.10.2



### Change to target directory
cd $1

### Make sure to move and modify these files to the target directory
# Change Layout File
layoutfile=layoutfile.csv
# FASTQ information
libraries=libraries.csv
# Hashtag file
feature_refs=feature_refs_totalseq_2023.csv


# CellRanger Pipeline

    #____1. Generate fastq files . #for --run, put location of BroadRun folder with Data
echo "Running cellranger mkfastq".
cellranger mkfastq --id=fastqs --run=/n/groups/cbdm_lab/scRNAseq/BroadRun240919/1726591964 --csv=$layoutfile


    #____2. Generate counts data from fastq files
echo "Running cellranger count"
cellranger count  --id=combined_sample \
                  --libraries=$libraries \
                  --feature-ref=$feature_refs \
                  --transcriptome=/n/groups/cbdm-db/bv43/M25/M25 \
		  --chemistry=SC5P-R2 \
                  --create-bam=true \
                  --expect-cells=100000


    #____3. Generate processed data for TCR (using IMGT B6 reference)
echo "Running cellranger vdj"
cellranger vdj --id=TCR \
                 --reference=/n/groups/cbdm_lab/scRNAseq/files/vdj_IMGT_B6_ref \
                 --fastqs=fastqs/outs/fastq_path/ \
                 --sample=VDJ

