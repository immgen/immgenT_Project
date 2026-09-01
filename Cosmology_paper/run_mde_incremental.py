#!/project/zemmour/david/envs/scvi120_20241008/bin/python
"""This script updates a reference MDE with a new integration with a query dataset"""
#author: David Zemmour
#date: 10/08/2024
# run_mde_incremental.py [cwd] [prefix] [mde_ref_file] [totalvi_integrated_file]

import warnings; warnings.simplefilter('ignore')
import argparse
parser = argparse.ArgumentParser(description="run TOTALVI from a MuData object specifying batch_key and categorical_covariate_keys if needed, can also return corrected counts and denoised data")
parser.add_argument('--working_dir', help='Working directory')
# parser.add_argument('--path_to_mudata', help='Path to MuData object')
parser.add_argument('--prefix', default='myprefix', 
                    help='Prefix for the output files (default: myprefix)')
parser.add_argument('--mde_ref_file', default=None, help='csv file with reference MDE')
parser.add_argument('--totalvi_integrated_file', default=None, help='csv file with reference and query integrated latent space')

print("Arguments")
args = parser.parse_args()
working_dir = args.working_dir
# path_to_mudata = args.path_to_mudata
prefix = args.prefix
mde_ref_file = args.mde_ref_file
totalvi_integrated_file = args.totalvi_integrated_file

# working_dir='/project/zemmour/david/ImmgenT/analysis/data_integration/IGT1_96/CD4'
# #path_to_mudata='/project/zemmour/david/ImmgenT/analysis/data_integration/IGT1_96/export_data/igt1_96_20241006.h5mu'
# path_to_mudata='/project/zemmour/david/ImmgenT/analysis/data_integration/IGT1_96/export_data/igt1_96_20241113_CD4.h5mu'
# prefix='totalvi_20241113_CD4_rmIGTsample'
# mde_ref_file='totalvi_20241014_CD4_rmIGTsample/mde2.csv'
# totalvi_integrated_file='totalvi_20241113_CD4_rmIGTsample/latent.csv'

print(f"Working Directory: {working_dir}")
# print(f"Path to AnnData: {path_to_mudata}")
print(f"Prefix: {prefix}")
print(f"MDE reference file: {mde_ref_file}")
print(f"totalvi_integrated_file: {totalvi_integrated_file}")


print("Importing libraries")
import warnings; warnings.simplefilter('ignore')
import scvi
import os
import sys
import scanpy as sc
import muon as mu
import mudata as mu

import numpy as np
import mplscience
import matplotlib.pyplot as plt
import pickle
import pandas as pd
import seaborn as sns
import torch
import pymde #to run MDE
pymde.seed(0)

print("Global configurations")
pd.set_option('display.max_rows', 1000)
pd.set_option('display.max_columns', 1000)
# sc.set_figure_params(figsize=(6, 6), frameon=False)
scvi.settings.seed = 0  # optional: ensures reproducibility
pymde.seed(0)
#sns.set_theme()
if torch.cuda.is_available():
    print("CUDA is available")
    print("Using device:", torch.cuda.get_device_name())
    torch.set_float32_matmul_precision("high")

print("Reading integrated ref+query latent space")  
totalvi_integrated = pd.read_csv(totalvi_integrated_file, index_col = 0, header = 0)

print("anchor MDE")
##test dataset to troubleshoot
## mde_ref_embedding_test = mde_ref_embedding[:5000]
## rownames_tensor_test = rownames_tensor[:5000]
## common_indices = totalvi_integrated.index.intersection(mde_ref_embedding_test.index)
## other_indices = totalvi_integrated.index.difference(mde_ref_embedding_test.index)
## random_sample_indices = np.random.choice(other_indices, size=1000, replace=False)
## combined_indices = common_indices.union(random_sample_indices)
## totalvi_integrated_test = totalvi_integrated.loc[combined_indices]

# match mde_ref_embedding with index in totalvi_integrated
mde_ref_embedding = pd.read_csv(mde_ref_file, index_col = 0, header = 0)
mde_ref_embedding = mde_ref_embedding.loc[mde_ref_embedding.index.intersection(totalvi_integrated.index)] #subset cells present in the integrated space other wise mismatch between anchors and values in pymde.Anchored
rownames = mde_ref_embedding.index.tolist()  # Convert index to list
index_positions = [totalvi_integrated.index.get_loc(item) for item in rownames if item in totalvi_integrated.index]
len(index_positions)
len(mde_ref_embedding.index)
rownames_tensor = torch.tensor(index_positions, device='cuda:0')

print("pymde.Anchored")
anchor_constraint = pymde.Anchored(
    anchors=rownames_tensor,
    values=torch.tensor(mde_ref_embedding.values, dtype=torch.float32,device='cuda:0'), #device='cpu' 'cuda:0'
)

print("pymde.preserve_neighbors")
incremental_mde = pymde.preserve_neighbors(
    torch.tensor(totalvi_integrated.values, dtype=torch.float32,device='cuda:0'), #device='cpu'
    embedding_dim=2,
    constraint=anchor_constraint,
    repulsive_fraction=0.7,
    n_neighbors = 15,
    verbose=True,
    device='cuda:0'
)

#if needed to run without anchoring
## incremental_mde = pymde.preserve_neighbors(
##     torch.tensor(totalvi_integrated.values, dtype=torch.float32), #device='cpu'
##     constraint=anchor_constraint,
##     init='random',
##     verbose=True #,device = 'cpu'
##     )

print("incremental_mde.embed")
incremental_mde.embed(eps=1e-6, verbose=True)

pymde.plot(mde_ref_embedding.values)
fig = plt.gcf()
fig.savefig(prefix+"/mde_incremental_original.png")
pymde.plot(incremental_mde.X)
fig = plt.gcf()
fig.savefig(prefix+"/mde_incremental.png")

model_pkl_file = prefix+'/mde_incremental_model.pkl'
with open(model_pkl_file, 'wb') as file:
    pickle.dump(incremental_mde, file)
 
# model_pkl_file = prefix+'/mde_incremental_model.pkl'
# with open(model_pkl_file, 'rb') as file:
#     incremental_mde = pickle.load(file)

incremental_mde_df = pd.DataFrame(incremental_mde.X.cpu().numpy(), index = totalvi_integrated.index)
incremental_mde_df.to_csv(prefix+'/mde_incremental.csv' , index=True)
