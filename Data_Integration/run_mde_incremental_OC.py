#!/project/zemmour/david/envs/scvi120_20241008/bin/python
"""This script updates a reference MDE with a new integration with a query dataset"""
#author: David Zemmour
#date: 10/08/2024
# run_mde_incremental.py [cwd] [prefix] [mde_ref_file] [totalvi_integrated_file]

import warnings; warnings.simplefilter('ignore')
#import argparse
import sys
#parser = argparse.ArgumentParser(description="run TOTALVI from a MuData object specifying batch_key and categorical_covariate_keys if needed, can also return corrected counts and denoised data")
#parser.add_argument('--working_dir', help='Working directory')
## parser.add_argument('--path_to_mudata', help='Path to MuData object')
#parser.add_argument('--prefix', default='myprefix', 
#                    help='Prefix for the output files (default: myprefix)')
#parser.add_argument('--mde_ref_file', default=None, help='csv file with reference MDE')
#parser.add_argument('--totalvi_integrated_file', default=None, help='csv file with reference and query integrated latent space')

print("Arguments")
#args = parser.parse_args()
working_dir = sys.argv[1] #  args.working_dir
# path_to_mudata = args.path_to_mudata
prefix = sys.argv[2] # args.prefix
mde_ref_file = sys.argv[3] # args.mde_ref_file
latent_df_scvi_abT = sys.argv[4] # args.totalvi_integrated_file
latent_df_scvi_gdT = sys.argv[5] # args.totalvi_integrated_file
latent_level1_abT = sys.argv[6]
output_file = sys.argv[7]

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
#print(f"totalvi_integrated_file: {totalvi_integrated_file}")


print("Importing libraries")
import warnings; warnings.simplefilter('ignore')
#import scvi
import os
#import sys
#import scanpy as sc
#import muon as mu
#import mudata as mu

import numpy as np
#import mplscience
import matplotlib.pyplot as plt
import pickle
import pandas as pd
#import seaborn as sns
import torch
import pymde #to run MDE
import re
pymde.seed(0)

print("Global configurations")
pd.set_option('display.max_rows', 1000)
pd.set_option('display.max_columns', 1000)
# sc.set_figure_params(figsize=(6, 6), frameon=False)
#scvi.settings.seed = 0  # optional: ensures reproducibility
pymde.seed(0)
#sns.set_theme()
if torch.cuda.is_available():
    print("CUDA is available")
    print("Using device:", torch.cuda.get_device_name())
    torch.set_float32_matmul_precision("high")

print("run pyMDE - first on all data, second on subgroups")
print("pyMDE on all data")
##latent_df = mdata.obsm[LEVEL1_SCANVI_KEY]
#latent_df = latent_df_scvi
## Load the full MDE reference embedding
latent_level1_abT = pd.read_csv(latent_level1_abT, index_col=0)
latent_df_scvi_abT = pd.read_csv(latent_df_scvi_abT, index_col=0)
latent_df_scvi_gdT = pd.read_csv(latent_df_scvi_gdT, index_col=0)

## Remove All_IGT duplicates
#latent_df_scvi_abT = latent_df_scvi_abT[~latent_df_scvi_abT.index.duplicated(keep='first')]
#latent_df_scvi_gdT = latent_df_scvi_gdT[~latent_df_scvi_gdT.index.duplicated(keep='first')]
dup_rows_abT_level1 = latent_level1_abT.duplicated()
latent_level1_abT = latent_level1_abT[~dup_rows_abT_level1]

dup_rows_abT = latent_df_scvi_abT.duplicated()
latent_df_scvi_abT =latent_df_scvi_abT[~dup_rows_abT]
dup_rows_gdT = latent_df_scvi_gdT.duplicated()
latent_df_scvi_gdT = latent_df_scvi_gdT[~dup_rows_gdT]

cells = latent_df_scvi_gdT.index.union(latent_df_scvi_abT.index)


mde_ref_embedding = pd.read_csv(mde_ref_file, index_col=0)
mde_ref_embedding = mde_ref_embedding.loc[mde_ref_embedding.index.intersection(cells)]

## Make level1 and level2 mde files
## level1
mde_ref_embedding_level1 = mde_ref_embedding[["level1_MDE1","level1_MDE2"]]
mde_ref_embedding_level1.index = mde_ref_embedding.index.copy()

## level2
mde_ref_embedding_level2 = mde_ref_embedding[["level2_MDE1","level2_MDE2"]]
mde_ref_embedding_level2.index = mde_ref_embedding.index.copy()

#latent_df = latent_df_scvi

## abT
## level1 mde
mde_ref_embedding_level1_abT = mde_ref_embedding_level1.loc[mde_ref_embedding_level1.index.intersection(latent_level1_abT.index)]
rownames = mde_ref_embedding_level1_abT.index.tolist()  # Convert index to list
index_positions = [latent_level1_abT.index.get_loc(item) for item in rownames if item in latent_level1_abT.index]
print(index_positions)
print(len(index_positions))
print(len(mde_ref_embedding_level1_abT.index))
rownames_tensor = torch.tensor(index_positions, device='cpu') # , device='cuda:0')

print("pymde.Anchored")
anchor_constraint = pymde.Anchored(
    anchors=rownames_tensor,
    values=torch.tensor(mde_ref_embedding_level1_abT.values, dtype=torch.float32, device='cpu'), # device='cuda:0),' #device='cpu' 'cuda:0'
)

print("pymde.preserve_neighbors")
incremental_mde_abT = pymde.preserve_neighbors(
    torch.tensor(latent_level1_abT.values, device='cpu'), # dtype=torch.float32, device='cpu'), #device='cuda:0'), #device='cpu'
    embedding_dim=2,
    constraint=anchor_constraint,
    repulsive_fraction=0.7,
    n_neighbors = 15,
    verbose=True,
    device='cpu'  # 'cuda:0'
)

print("incremental_mde.embed")
incremental_mde_abT.embed(eps=1e-6, verbose=True)
level1_abT_mde_df = pd.DataFrame(incremental_mde_abT.X.cpu().numpy(), index = latent_level1_abT.index)

output_file = pd.read_csv(output_file, index_col=0)
common_ids = output_file.index.intersection(level1_abT_mde_df.index)
level1_abT_mde_df = level1_abT_mde_df.loc[common_ids, :]

print(f"Found {len(common_ids)} common IDs")

output_file['allT_MDE1'] =  -3.0
#if len(common_ids) > 0:
#    output_file.loc[common_ids, 'allT_MDE1'] = level1_abT_mde_df.loc[common_ids, 0].values
#else:
#    print("⚠️ No matching cell IDs between output_file and MDE data.")
output_file.loc[common_ids, 'allT_MDE1'] =  level1_abT_mde_df.loc[common_ids, 0].values

#output_file.loc[latent_df_scvi_abT.index, 'level2_MDE1'] =  level1_abT_mde_df[0]
output_file['level2_MDE1'] =  -3.0
output_file['allT_MDE2'] =  -3.0
output_file.loc[common_ids, 'allT_MDE2'] =  level1_abT_mde_df.loc[common_ids, 1].values
#output_file.loc[latent_df_scvi_abT.index, 'level2_MDE2'] =  level1_abT_mde_df[1]
output_file['level2_MDE2'] =  -3.0

## gdT
## level1 mde
mde_ref_embedding_level1_gdT = mde_ref_embedding_level1.loc[mde_ref_embedding_level1.index.intersection(latent_df_scvi_gdT.index)]
rownames = mde_ref_embedding_level1_gdT.index.tolist()  # Convert index to list
index_positions = [latent_df_scvi_gdT.index.get_loc(item) for item in rownames if item in latent_df_scvi_gdT.index]
print(index_positions)
print(len(index_positions))
print(len(mde_ref_embedding_level1_gdT.index))
rownames_tensor = torch.tensor(index_positions, device='cpu') # , device='cuda:0')

print("pymde.Anchored")
anchor_constraint = pymde.Anchored(
    anchors=rownames_tensor,
    values=torch.tensor(mde_ref_embedding_level1_gdT.values, dtype=torch.float32, device='cpu'), # device='cuda:0),' #device='cpu' 'cuda:0'
)

print("pymde.preserve_neighbors")
incremental_mde_gdT = pymde.preserve_neighbors(
    torch.tensor(latent_df_scvi_gdT.values, device='cpu'), # dtype=torch.float32, device='cpu'), #device='cuda:0'), #device='cpu'
    embedding_dim=2,
    constraint=anchor_constraint,
    repulsive_fraction=0.7,
    n_neighbors = 15,
    verbose=True,
    device='cpu'  # 'cuda:0'
)

print("incremental_mde.embed")
incremental_mde_gdT.embed(eps=1e-6, verbose=True)
level1_gdT_mde_df = pd.DataFrame(incremental_mde_gdT.X.cpu().numpy(), index = latent_df_scvi_gdT.index)
common_ids = output_file.index.intersection(level1_gdT_mde_df.index)
level1_gdT_mde_df = level1_gdT_mde_df.loc[common_ids, :]

#output_file.loc['allT_MDE1'] =  -3
output_file.loc[common_ids, 'allT_MDE1'] =  level1_gdT_mde_df.loc[common_ids, 0].values
#output_file.loc[latent_df_scvi_abT.index, 'level2_MDE1'] =  level1_abT_mde_df[0]
#output_file.loc['level2_MDE1'] =  -3
#output_file.loc['allT_MDE2'] =  -3
output_file.loc[common_ids, 'allT_MDE2'] =  level1_gdT_mde_df.loc[common_ids, 1].values
#output_file.loc[latent_df_scvi_abT.index, 'level2_MDE2'] =  level1_abT_mde_df[1]
#output_file.loc['level2_MDE2'] =  -3

output_file.to_csv(prefix+"/predictions_level1_output_file.csv", index = True)

print("pyMDE on subgroups")
## Subset to level1 annotations
##level1_annotations = mdata.obs["level1_transfer_labels"].unique()
##level1_annotations = ["CD4", "CD8", "Treg", "gdT", "CD8aa", "nonconv", "DN", "DP"]
threshold = 20
LEVEL1_SCANVI_PREDICTIONS_KEY="level1_final"
counts = output_file[LEVEL1_SCANVI_PREDICTIONS_KEY].value_counts()
level1_annotations = counts[counts > threshold].index.tolist()
level1_annotations = [x for x in level1_annotations if x not in ["unclear", "nonT", "remove", "not classified", "thymocyte"]]
#latent_df = level2_latent_df

## Directory to save outputs
##os.makedirs(f"{prefix}/mde_by_level1", exist_ok=True)

for annotation in level1_annotations:
    print(f"\n--- Processing group: {annotation} ---")
    
    ## Cells in the current group
    ##group_cells = mdata.obs_names[mdata.obs[LEVEL1_SCANVI_PREDICTIONS_KEY] == annotation]
    group_cells = output_file.index[output_file[LEVEL1_SCANVI_PREDICTIONS_KEY] == annotation]
    
    latent_df = latent_df_scvi_abT
    if annotation == "gdT":
        latent_df = latent_df_scvi_gdT

    ## Subset the totalVI embedding to group
    
    #group_integrated = latent_df.loc[latent_df.index.intersection(group_cells)]
    
    ## Also subset reference anchors to those in this group
    #group_ref_embedding = mde_ref_embedding_level2.loc[mde_ref_embedding_level2.index.intersection(group_integrated.index)]
    
    ## Get all IGT-prefixed cells (IGT1..IGT96 but not All_IGT)
    IGT_pattern = re.compile(r"^(?!All_IGT)(IGT[1-9][0-9]?)")  # matches IGT1–IGT99 but not All_IGT
    IGT_cells = [cell for cell in latent_df.index if IGT_pattern.match(cell)]

    ## Union of both
    union_cells = set(group_cells).union(IGT_cells)

    ## Subset embeddings
    group_integrated = latent_df.loc[latent_df.index.intersection(union_cells)]
    group_ref_embedding = mde_ref_embedding_level2.loc[mde_ref_embedding_level2.index.intersection(group_integrated.index)]

    ## Check if we have enough anchors to continue
    if group_ref_embedding.shape[0] < 5:
        print(f"Skipping {annotation} (too few anchors: {group_ref_embedding.shape[0]})")
        continue
    
    if group_integrated.shape[0] < 1:
        print(f"Skipping {annotation} (too few {annotation} cells: {group_integrated.shape[0]})")
        continue


    ## Get positions of anchor rows in SCVI
    index_positions = [
        group_integrated.index.get_loc(cell)
        for cell in group_ref_embedding.index if cell in group_integrated.index
    ]
    
    if len(index_positions) != group_ref_embedding.shape[0]:
        print("Mismatch between index positions and anchors — skipping.")
        continue

    ## Prepare PyTorch tensors
    rownames_tensor = torch.tensor(index_positions, device='cpu')
    anchor_values = torch.tensor(group_ref_embedding.values, dtype=torch.float32, device='cpu')
    full_tensor = torch.tensor(group_integrated.values, dtype=torch.float32, device='cpu')

    ## Anchored MDE constraint
    print("  > Building anchor constraint")
    anchor_constraint = pymde.Anchored(anchors=rownames_tensor, values=anchor_values)

    ## Run MDE
    if annotation in ["DN", "DP", "nonconv"]:
        init_method = "random"
    else:
        init_method = "quadratic"

    print("  > Running pymde.preserve_neighbors")
    mde_model = pymde.preserve_neighbors(
        full_tensor,
        embedding_dim=2,
        constraint=anchor_constraint,
        repulsive_fraction=0.7,
        n_neighbors=15,
        verbose=True,
        device='cpu',
        init=init_method
    )

    print("  > Embedding")
    mde_model.embed(eps=1e-6, verbose=True)

    ## Add result to output file
    level2_subgroup_mde_df = pd.DataFrame(mde_model.X.cpu().numpy(), index=group_integrated.index)
    common_ids = output_file.index.intersection(group_integrated.index)
    output_file.loc[common_ids, 'level2_MDE1'] = level2_subgroup_mde_df.loc[common_ids, 0].values
    output_file.loc[common_ids, 'level2_MDE2'] = level2_subgroup_mde_df.loc[common_ids, 1].values
    ##out_path = f"{prefix}/mde_by_level1/mde_{annotation.replace(' ', '_')}.csv"
    ##output_df.to_csv(out_path)
    ##print(f"  > Saved to {out_path}")


print("Prepare output files - plotting and annotations for user")
##user_output_file = output_file.loc[query_mask.obs.index, :] 
output_file.to_csv(prefix+"/output_file.csv", index=True)
##user_output_file.index = query_mask.obs.index
##user_output_file.to_csv(prefix+"/user_output_file.csv", index=True)
