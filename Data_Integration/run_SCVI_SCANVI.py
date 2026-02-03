#!/usr/bin/env python3
"""This script run TOTALVI from a Mudata object specifying batch_key and categorical_covariate_keys if needed, """
#author: David Zemmour
#date: 10/08/2024
#run_totalvi_v2.py [cwd] [path to mudata .h5mu] [prefix] [batchkey] [categorical_covariate_keys] [corrected_counts] [denoised_data]

import warnings; warnings.simplefilter('ignore')
import argparse
parser = argparse.ArgumentParser(description="run TOTALVI from a MuData object specifying batch_key and categorical_covariate_keys if needed, can also return corrected counts and denoised data")
parser.add_argument('--working_dir', help='Working directory')
parser.add_argument('--path_to_ImmgenT', help='Path to ImmgenT MuData object')
parser.add_argument('--path_to_matrix', help='Path to query matrix')
parser.add_argument('--path_to_spikein', help='Path to combined AnnData object')
parser.add_argument('--query_sample_path', help='Path to query sample table')
parser.add_argument('--prefix', default='myprefix',
                    help='Prefix for the output files (default: myprefix)')
parser.add_argument('--batchkey', default=None, help='Batch key for analysis')
parser.add_argument('--categorical_covariate_keys', default=None, help='categorical_covariate_keys variables (default: None)')
parser.add_argument('--corrected_counts', default=False, help='Returns corrected counts, aka posterior_predictive_sample() (default: False)')
parser.add_argument('--denoised_data', default=False, help='Returns denoised data, aka get_normalized_expression()  (default: False)')
parser.add_argument('--mde_ref_file', help='mde_ref_file  (default: None)')
#parser.add_argument('--latent_key', help='Key for latent space')

print("Arguments")
args = parser.parse_args()
working_dir = args.working_dir
path_to_ImmgenT = args.path_to_ImmgenT
path_to_matrix = args.path_to_matrix
path_to_spikein = args.path_to_spikein
query_sample_path = args.query_sample_path
#path_to_query = args.path_to_query
prefix = args.prefix
batchkey = args.batchkey
#confoundings = args.confoundings
if args.categorical_covariate_keys:
        # Split the string into a list by commas
        categorical_covariate_keys = args.categorical_covariate_keys.split(',')
        #print(categorical_covariate_keys)
else:
    categorical_covariate_keys = None
corrected_counts = args.corrected_counts
denoised_data = args.denoised_data
#mde_ref_file = args.mde_ref_file
#latent_key = args.latent_key
# working_dir = '/project/jfkfloor2/zemmourlab/david/immgent/analysis/integration/IGT1_56/'
# path_to_mudata = '/project/jfkfloor2/zemmourlab/david/immgent/analysis/integration/IGT1_56/export_data/totalvi_igt1_56_20231030_allgenes_mdata.h5mu'
# prefix = 'totalvi_igt1_56_allgenes_20240526_igtsampleregressedout'
# batchkey = 'IGT'
# categorical_covariate_keys = ['sample_id']

# working_dir='/project/zemmour/david/ImmgenT/analysis/data_integration/IGT1_96'
# path_to_mudata='/project/zemmour/david/ImmgenT/analysis/data_integration/IGT1_96/export_data/igt1_96_20241006.h5mu'
# prefix='totalvi_20241006' #prefix = "totalvi_20241008_rmIGTsample"
# batchkey='IGT'
# categorical_covariate_keys='IGTHT'
# corrected_counts=False
# denoised_data=False


print(f"Working Directory: {working_dir}")
print(f"Path to ImmgenT AnnData: {path_to_ImmgenT}")
print(f"Path to matrix: {path_to_matrix}")
print(f"Path to spike-in AnnData: {path_to_spikein}")
#print(f"Path to query AnnData: {path_to_query}")
#print(f"Path to combined AnnData: {path_to_anndata}")
print(f"Prefix: {prefix}")
print(f"Batch Key: {batchkey}")
print(f"categorical_covariate_keys: {categorical_covariate_keys}")
print(f"corrected_counts: {corrected_counts}")
print(f"denoised_data: {denoised_data}")
#print(f"mde_ref_file: {mde_ref_file}")
#print(f"Latent Key: {latent_key}")

print("Importing libraries")
import warnings; warnings.simplefilter('ignore')
import scvi
import os
import sys
import scanpy as sc
#import muon as mu
import mudata as mu
import anndata as AnnData
import numpy as np
#import mplscience
import scipy
import scipy.io
from scipy.sparse import csc_matrix
import matplotlib.pyplot as plt
import pickle
import pandas as pd
import seaborn as sns
import torch
import fast_matrix_market as fmm
import pymde #to run MDE
from annoy import AnnoyIndex
from collections import Counter
import re
from pathlib import Path


print("Global configurations")
pd.set_option('display.max_rows', 1000)
pd.set_option('display.max_columns', 1000)
# sc.set_figure_params(figsize=(6, 6), frameon=False)
scvi.settings.seed = 0  # optional: ensures reproducibility
#pymde.seed(0)
#sns.set_theme()
if torch.cuda.is_available():
    print("CUDA is available")
    print("Using device:", torch.cuda.get_device_name())
    torch.set_float32_matmul_precision("high")

print("Reading mudata")
os.chdir(working_dir)
##mdata = mu.read(path_to_mudata) #totalvi_igt1_56_allgenes_Treg_20240306/adata.h5mu") #mu.read("totalvi_igt1_56_allgenes_Treg_20240306/adata.h5mu
##print(mdata)
#mdata = AnnData.read(path_to_anndata) #load this one (results of this block)
#mdata.raw = None
## ImmgenT
ImmgenT_mdata = AnnData.read(path_to_ImmgenT)
## Keep only IGT1_96
# Define the values to remove
exclude_IGT_values = [f'IGT{i}' for i in range(97, 105)]  # 105 is exclusive
All_IGT_values = [f'IGT{i}' for i in range(1, 97)]
# Filter out those cells
ImmgenT_mdata = ImmgenT_mdata[~ImmgenT_mdata.obs['IGT'].isin(exclude_IGT_values)].copy()

ImmgenT_mdata.X = ImmgenT_mdata.X.copy()
ImmgenT_mdata.layers["counts"] = ImmgenT_mdata.X.copy()
print(ImmgenT_mdata)
## Read in query Anndata object
#query_mdata = AnnData.read(path_to_query)

print("Create anndata object from matrix")
##matrix = scipy.io.mmread(path_to_matrix+"/matrix.mtx.gz").tocsr()
##matrix = matrix.transpose().tocsr()
##genes = pd.read_csv(path_to_matrix+"/features.tsv.gz", sep = "\t", header=None, compression='infer')
##cells = pd.read_csv(path_to_matrix+"/barcodes.tsv.gz", sep = "\t", header=None, compression='infer')

p = Path(path_to_matrix)

def find_file(stem, exts=("", ".gz")):
    for ext in exts:
        f = p / f"{stem}{ext}"
        if f.exists():
            return f
    raise FileNotFoundError(f"Could not find {stem}[.gz] in {p}")

# locate files
matrix_file   = find_file("matrix.mtx")
#barcodes_file = find_file("barcodes.tsv")
barcodes_file = (
    find_file("barcodes.tsv") if (p / "barcodes.tsv").exists() or (p / "barcodes.tsv.gz").exists()
    else find_file("cells.tsv")
)

features_file = (
    find_file("features.tsv") if (p / "features.tsv").exists() or (p / "features.tsv.gz").exists()
    else find_file("genes.tsv")
)

# read data
matrix = scipy.io.mmread(matrix_file).tocsr().T
cells  = pd.read_csv(barcodes_file, sep="\t", header=None)
genes  = pd.read_csv(features_file, sep="\t", header=None)

genes.columns = ['gene_id', 'gene_name'] #, 'feature_type']
# Make genes column unique
def make_unique(names):
    counts = {}
    result = []
    for name in names:
        if name not in counts:
            counts[name] = 1
            result.append(name)
        else:
            new_name = f"{name}_{counts[name]}"
            while new_name in counts:
                counts[name] += 1
                new_name = f"{name}_{counts[name]}"
            result.append(new_name)
            counts[new_name] = 1
    return result

# Apply to gene_name column
genes['gene_name'] = make_unique(genes['gene_name'])

cells.columns = ['barcode']

query_mdata = AnnData.AnnData(X=matrix, 
        obs=pd.DataFrame(index=cells["barcode"].values), 
        var=pd.DataFrame(index=genes["gene_name"].values))
#del mdata
query_mdata.X = query_mdata.X.copy()
query_mdata.layers["counts"] = query_mdata.X.copy()
print("Pre gene filtering: ")
print(query_mdata)
cells_nGenes = mdata.obs.index[mdata.obs['n_genes_by_counts' < 300]
sc.pp.filter_cells(query_mdata, min_genes=300)
print("Post gene filtering: ")
print(query_mdata)
## Add sample and hashtag information to query mdata before spike-in merge
query_sample_table = pd.read_csv(query_sample_path, index_col=0)
#query_sample_table.index = query_sample_table["Cell_ID"]
joint_bcs = query_sample_table.index.intersection(query_mdata.obs_names)
query_sample_table = query_sample_table.loc[joint_bcs, :].copy()
query_mdata = query_mdata[joint_bcs, :].copy()
query_mdata.obs["IGT"] = query_sample_table["dataset"]
query_mdata.obs["HT"] = query_sample_table["sample"]
query_mdata.obs["IGTHT"] = query_mdata.obs["IGT"].astype(str) + "_" + query_mdata.obs["HT"].astype(str)
query_mdata.obs.to_csv(prefix+"/query_all_metadata.csv", index = True)

spikein_mdata = AnnData.read(path_to_spikein)
query_mdata = AnnData.concat(
    [query_mdata, spikein_mdata],
    axis=0,
    join="inner"
)

print("query with spike-in: ")
print(query_mdata)

##full.mod['RNA'].obs_names = full_RNA.obs_names.astype(str)
##full.mod['RNA'].var_names = genes['x'].astype(str)
##full.mod['ADT'].obs_names = full_ADT.obs_names.astype(str)
##full.mod['ADT'].var_names = proteins.astype(str)

##ImmgenT_IGTHT = pd.read_csv(working_dir+"ImmgenT_IGTHT.csv", delimiter = ",")
##ImmgenT_IGTHT.index = ImmgenT_IGTHT['Unnamed: 0'].values
##ImmgenT_IGTHT = ImmgenT_IGTHT.drop(columns = "Unnamed: 0")
##ImmgenT_mdata.obs['IGTHT'] = ImmgenT_IGTHT["IGTHT"].astype(str)
##ImmgenT_mdata.obs['IGT'] = ImmgenT_IGTHT["IGT"].astype(str)


##query_IGTHT = pd.read_csv(working_dir+"query_IGTHT.csv", delimiter = ",")
##query_IGTHT.index = query_IGTHT['Unnamed: 0'].values
##query_IGTHT = query_IGTHT.drop(columns = "Unnamed: 0")
##query_IGTHT['Unknown'] = "Unknown"

print("Concat anndata - keep only common genes")
mdata = AnnData.concat(
    [ImmgenT_mdata, query_mdata],
    axis=0,
    join="inner",
    label="origin",
    keys=["ImmgenT", "query"]
)

#print("Cells pre gene filtering:")
#print(mdata)

## Light QC filtering for nFeatures
##sc.pp.calculate_qc_metrics(mdata, inplace=True)
##mdata = mdata[mdata.obs['n_genes_by_counts'] > 300].copy()
#sc.pp.filter_cells(mdata, min_genes=300)

#print("Cells post gene filtering:")
#print(mdata)


## Add T Cell score and filter based on this
## Read in gene signature list - select T cell genes
signature_list = pd.read_csv("/n/groups/cbdm_lab/odc180/ImmgenT_workshop/Integration_Webpage/Scripts_TRBI/LineageSpecGns072018_top27.csv", header = 0)
Tcell_genes = signature_list.loc[signature_list["Type"] == "T", "Marker"].astype(str).tolist()
#Tcell_genes_clean = [g for g in Tcell_genes if g in mdata.var_names]
#print(Tcell_genes_clean)

## Add T cell score based on expression of genes
## Extract expression for those genes only
#Tcell_subset = mdata[:, Tcell_genes_clean].X.copy()

## If sparse, convert to dense
#if not isinstance(Tcell_subset, np.ndarray):
#    Tcell_subset = Tcell_subset.toarray()

## Total score per cell (sum over genes)
#mdata.obs['Tcell_gene_score'] = Tcell_subset.sum(axis=1)

# Subset based on this score
query_cells = query_mdata.obs_names
ImmgenT_cells = ImmgenT_mdata.obs_names
#Tcell_mask = (mdata.obs_names.isin(query_cells)) & (mdata.obs['Tcell_gene_score'] > 10)

## Normalize first to use scoring 
for b in mdata.obs['IGT'].unique():
    sc.pp.normalize_total(mdata[mdata.obs['IGT'] == b], target_sum=1e4)
    sc.pp.log1p(mdata[mdata.obs['IGT'] == b])

# 3. Compute gene scores (after normalization but before batch correction)
sc.tl.score_genes(mdata, gene_list=Tcell_genes, score_name='Tcell_gene_module_score')

print("Pre Tcell filter: ")
print(mdata)
Tcell_mask = (
    (mdata.obs_names.isin(query_cells) & (mdata.obs['Tcell_gene_module_score'] > 0))
    | mdata.obs_names.isin(ImmgenT_cells)
)

cells_Tcell_score = mdata.obs.index[mdata.obs['Tcell_gene_module_score'] <= 0]
mdata = mdata[Tcell_mask, :].copy()

print("save QC filtered cells (if any)")
dfs = []

if cells_Tcell_score is not None:
    dfs.append(pd.DataFrame({"cell_id": cells_Tcell_score, "QC_filter": "non T cell"}))

if cells_nGenes is not None:
    dfs.append(pd.DataFrame({"cell_id": cells_nGenes, "QC_filter": "low gene count"}))

if dfs:
    filtered_cells = pd.concat(dfs, axis=0, ignore_index=True)
    filtered_cells.to_csv("QC_filtered_cells.csv", index=False)

print("Remaining Tcells: ")
print(mdata)

## Remove TCR genes from anndata object
genes = mdata.var_names

# Identify TCR-related genes using regex
##Trav = [g for g in genes if re.match(r'^Trav', g)]
##Traj = [g for g in genes if re.match(r'^Traj', g)]
##Trac = [g for g in genes if re.match(r'^Trac', g)]
##Trbv = [g for g in genes if re.match(r'^Trbv', g)]
##Trbd = [g for g in genes if re.match(r'^Trbd', g)]
##Trbj = [g for g in genes if re.match(r'^Trbj', g)]
Trd  = [g for g in genes if re.match(r'^Trd',  g)]
Trg  = [g for g in genes if re.match(r'^Trg',  g)]

# Combine all TCR-related genes
##TCR_genes = Trav + Traj + Trac + Trbv + Trbd + Trbj + Trd + Trg
gdT_genes = Trd + Trg + ["Sox13"]

# List of sex-specific genes
##sex_specific_genes = ['Ddx3y', 'Uty', 'Xist', 'Eif2s3y', 'Kdm5d', 'Tsix']

# Full list of genes to remove
##remove_var_genes = set(TCR_genes + sex_specific_genes)

# Ensure remove_var_genes is a flat list of strings
##remove_var_genes = list(map(str, remove_var_genes))  # guarantees string type
##remove_var_genes = [g for g in remove_var_genes if g in mdata.var_names]

# Subset columns (genes/features) — note the [:, ...] for genes
##mdata = mdata[:, ~mdata.var_names.isin(remove_var_genes)].copy()
print(mdata)

## Add gdT score - two methods - 1. add counts, 2. score based on signature list
## 1. Add counts
## Extract expression for those genes only
#X_subset = mdata[:, gdT_genes].X.copy()
present_genes = [g for g in gdT_genes if g in mdata.var_names]
X_subset = mdata[:, present_genes].X.copy()

# If sparse, convert to dense
if not isinstance(X_subset, np.ndarray):
    X_subset = X_subset.toarray()

# Total score per cell (sum over genes)
mdata.obs['gdT_gene_score'] = X_subset.sum(axis=1)

## 2. Score based on signature list
sc.tl.score_genes(mdata, gene_list=gdT_genes, score_name='gdT_gene_module_score')


print("Creating RNA layer counts")
mdata.layers["counts"] = mdata.X.copy()

##mdata.obs['IGTHT'] = pd.concat([ImmgenT_IGTHT["IGTHT"],query_IGTHT["IGTHT"]]).astype(str)
##mdata.obs['IGT'] = pd.concat([ImmgenT_IGTHT["IGT"],query_IGTHT["IGT"]]).astype(str)
mdata.layers["counts"] = mdata.X.copy()
mdata.obs_names = mdata.obs_names.astype(str)
mdata.var_names = mdata.var_names.astype(str)
mdata.X = mdata.X.copy()
##mdata.var = mdata.var.reset_index(drop=True)
##mdata.mod['RNA'].var = mdata.mod['RNA'].var.reset_index(drop=True)
##mdata.mod['protein'].var = mdata.mod['protein'].var.reset_index(drop=True)
print(mdata)
##mdata.update()

print("Batch key list")
print(mdata.obs[batchkey].unique().tolist())
print("Batch key list: any NA?")
print((mdata.obs[batchkey].isna()).any())


print("categorical_covariate_keys")
if categorical_covariate_keys is not None:
    for c in categorical_covariate_keys:
        print(c)
        print(mdata.obs[c].unique().tolist())
        print("any NA?")
        print((mdata.obs[c].isna()).any())
        ##print(mdata.mod["RNA"].obs[c].head(10))
        ##mdata.mod["RNA"].obs[c] = mdata.mod["RNA"].obs[c].str.replace('.', '', regex=False)
        ##print(mdata.mod["RNA"].obs[c].head(10))

print("Train new SCVI model")
## Change path when updated annotation file is ready
annotation_level2 = pd.read_csv("/n/groups/cbdm_lab/odc180/ImmgenT_workshop/ImmgenT_freeze_20250109/igt1_104_withtotalvi20250505_annotation_table.csv", delimiter = ",", index_col=0)
annotation_level2 = annotation_level2.loc[annotation_level2.index.intersection(mdata.obs.index)]
mdata.obs['level1'] = "Unknown"
mdata.obs['level2'] = "Unknown"
##mdata.obs['level2.group'] = "Unknown"
mdata.obs.loc[annotation_level2.index, 'level1'] = annotation_level2['level1'].values
mdata.obs.loc[annotation_level2.index, 'level2'] = annotation_level2['level2'].values
##mdata.obs.loc[annotation_level2.index, 'level2.group'] = annotation_level2['level2.group'].values
query_cells = query_mdata.obs_names
gdT_mask = (mdata.obs_names.isin(query_cells)) & (mdata.obs['gdT_gene_module_score'] > 0)
mdata.obs.loc[gdT_mask, 'level1'] = "gdT"
print("query gdT labelled cells")
print(mdata.obs.loc[(mdata.obs_names.isin(query_cells)), 'level1'].value_counts())
print(mdata.obs['IGT'].value_counts())

## Separate gdT and abT mdata
abT_mdata = mdata[mdata.obs['level1'].astype(str) != "gdT" , :].copy()
gdT_mdata = mdata[mdata.obs['level1'] == "gdT" ,:].copy()

## Train separate SCVI models for gdT and abT
## abT
# Comment if already trained
scvi.model.SCVI.setup_anndata(abT_mdata, layer = "counts", batch_key = batchkey, categorical_covariate_keys = categorical_covariate_keys)
scvi_model_abT = scvi.model.SCVI(abT_mdata, n_latent=30, n_layers=2)
scvi_model_abT.train()
scvi_model_abT.save(prefix+"/scvi_model_abT/") #, save_anndata=True)
#scvi_model_abT = scvi.model.SCVI.load(prefix+"/scvi_model_abT/", adata=abT_mdata)

SCANVI_LATENT_KEY = "X_scANVI"
print("Save latent_representation.csv")
latent_representation = scvi_model_abT.get_latent_representation()
SCVI_LATENT_KEY = "X_scVI"
##mdata.mod['RNA'].obsm[SCVI_LATENT_KEY] = latent_representation
abT_mdata.obsm[SCVI_LATENT_KEY] = latent_representation
##latent_df = pd.DataFrame(latent_representation, index = mdata.mod['RNA'].obs.index)
latent_df_scvi_abT = pd.DataFrame(latent_representation, index = abT_mdata.obs.index)

latent_df_scvi_abT.to_csv(prefix+"/latent_df_scvi_abT.csv", index=True)

#print("Save umap.csv")
#sc.pp.neighbors(abT_mdata, use_rep=SCVI_LATENT_KEY)
#sc.tl.umap(abT_mdata, min_dist=0.4)
##umap_df = pd.DataFrame(mdata.obsm['X_umap'], index = mdata.mod['RNA'].obs.index)
#umap_df = pd.DataFrame(abT_mdata.obsm['X_umap'], index = abT_mdata.obs.index)
#umap_df.to_csv(prefix+"/umap_python.csv", index=True)


## gdT
# Comment if already trained
scvi.model.SCVI.setup_anndata(gdT_mdata, layer = "counts", batch_key = batchkey, categorical_covariate_keys = categorical_covariate_keys)
scvi_model_gdT = scvi.model.SCVI(gdT_mdata, n_latent=30, n_layers=2)
scvi_model_gdT.train()
scvi_model_gdT.save(prefix+"/scvi_model_gdT/") #, save_anndata=True)
#scvi_model_gdT = scvi.model.SCVI.load(prefix+"/scvi_model_gdT/", adata=gdT_mdata)

SCANVI_LATENT_KEY = "X_scANVI"
print("Save latent_representation.csv")
latent_representation = scvi_model_gdT.get_latent_representation()
SCVI_LATENT_KEY = "X_scVI"
##mdata.mod['RNA'].obsm[SCVI_LATENT_KEY] = latent_representation
gdT_mdata.obsm[SCVI_LATENT_KEY] = latent_representation
##latent_df = pd.DataFrame(latent_representation, index = mdata.mod['RNA'].obs.index)
latent_df_scvi_gdT = pd.DataFrame(latent_representation, index = gdT_mdata.obs.index)

latent_df_scvi_gdT.to_csv(prefix+"/latent_df_scvi_gdT.csv", index=True)

print("Save umap.csv")
#sc.pp.neighbors(gdT_mdata, use_rep=SCVI_LATENT_KEY)
#sc.tl.umap(gdT_mdata, min_dist=0.4)
##umap_df = pd.DataFrame(mdata.obsm['X_umap'], index = mdata.mod['RNA'].obs.index)
#umap_df = pd.DataFrame(gdT_mdata.obsm['X_umap'], index = gdT_mdata.obs.index)
##umap_df.to_csv(prefix+"/umap_python.csv", index=True)

print("Save Anndata")
mdata.write_h5ad(prefix+"/adata_RNA.h5ad")

##mdata = AnnData.read(prefix+"_orig/adata_RNA.h5ad")
##latent_df = pd.read_csv(prefix+"_orig/latent.csv", index_col=0)

print("Predict level1 and 2 annotations with SCANVI")
## Define ref and query masks (future use)
ref_mask = mdata.obs["origin"] == "ImmgenT"
query_mask = mdata.obs["origin"] == "query"

print("Pre level1 training: ")
print(mdata)

## abT only
##abT_mask = (mdata.obs_names.isin(query_cells)) & (mdata.obs['gdT_gene_module_score'] < 0.1) | (mdata.obs_names.isin(ImmgenT_cells)) & (mdata.obs['level1'] != "gdT")
##abT_mdata = mdata[mdata.obs['level1'].astype(str) != "gdT" , :].copy()
##abT_mdata = mdata[abT_mask, :].copy()
print("abT mdata:")
print(abT_mdata)

## level1
## Training scanvi model on scvi model
abT_labels_to_keep = abT_mdata.obs["level1"].unique().astype(str)
abT_mdata.obs["level1_abT"] = abT_mdata.obs["level1"].where(
    abT_mdata.obs["level1"].isin(abT_labels_to_keep),
    other="Unknown"
)

## Fix batch norm error - cant have a predicting label with < cells - 5 cells is the more stable count
label_counts = abT_mdata.obs["level1_abT"].value_counts()
drop = label_counts[label_counts < 5].index.tolist()     # labels with < 5 cells are considered drop
abT_mdata.obs.loc[abT_mdata.obs["level1_abT"].isin(drop), "level1_abT"] = "Unknown"

abT_mdata.obs["level1_abT"] = abT_mdata.obs["level1_abT"].astype("category")
# Setup anndata fresh so the label encoder is rebuilt
scvi.model.SCVI.setup_anndata(
    abT_mdata,
    layer="counts",
    batch_key=batchkey,
    categorical_covariate_keys=categorical_covariate_keys,
    labels_key="level1_abT"
)
level1_model_abT = scvi.model.SCANVI.from_scvi_model(scvi_model_abT, "Unknown", labels_key="level1_abT")
level1_model_abT.train(25)
level1_model_abT.save(prefix+"/abT_scanvi_level1_model/") #, save_anndata=True)
#level1_model_abT = scvi.model.SCANVI.load(prefix+"/abT_scanvi_level1_model/", adata=abT_mdata)


## level2
## Training scanvi model on scvi model
abT_labels_to_keep = abT_mdata.obs["level2"].unique().astype(str)
abT_mdata.obs["level2_abT"] = abT_mdata.obs["level2"].where(
    abT_mdata.obs["level2"].isin(abT_labels_to_keep),
    other="Unknown"
)

## Fix batch norm error - cant have a predicting label with < cells - 5 cells is the more stable count
label_counts = abT_mdata.obs["level2_abT"].value_counts()
drop = label_counts[label_counts < 5].index.tolist()     # labels with < 5 cells are considered drop
abT_mdata.obs.loc[abT_mdata.obs["level2_abT"].isin(drop), "level2_abT"] = "Unknown"

abT_mdata.obs["level2_abT"] = abT_mdata.obs["level2_abT"].astype("category")
# Setup anndata fresh so the label encoder is rebuilt
scvi.model.SCVI.setup_anndata(
    abT_mdata,
    layer="counts",
    batch_key=batchkey,
    categorical_covariate_keys=categorical_covariate_keys,
    labels_key="level2_abT"
)
level2_model_abT = scvi.model.SCANVI.from_scvi_model(scvi_model_abT, "Unknown", labels_key="level2_abT")
level2_model_abT.train(25)
level2_model_abT.save(prefix+"/abT_scanvi_level2_model/") #, save_anndata=True)
#level2_model_abT = scvi.model.SCANVI.load(prefix+"/abT_scanvi_level2_model/", adata=abT_mdata)

## gdT only
#gdT_mdata = mdata[mdata.obs['level1'] == "gdT" ,:].copy()
print("gdT mdata:")
print(gdT_mdata)
## level1
## Training scanvi model on scvi model
gdT_labels_to_keep = gdT_mdata.obs["level2"].unique().astype(str)
gdT_mdata.obs["level2_gdT"] = gdT_mdata.obs["level2"].where(
    gdT_mdata.obs["level2"].isin(gdT_labels_to_keep),
    other="Unknown"
)

## Fix batch norm error - cant have a predicting label with < cells - 5 cells is the more stable count
label_counts = gdT_mdata.obs["level2_gdT"].value_counts()
print(label_counts)
print(drop)
drop = label_counts[label_counts < 5].index.tolist()     # labels with < 5 cells are considered drop
gdT_mdata.obs.loc[gdT_mdata.obs["level2_gdT"].isin(drop), "level2_gdT"] = "Unknown"

gdT_mdata.obs["level2_gdT"] = gdT_mdata.obs["level2_gdT"].astype("category")
# Setup anndata fresh so the label encoder is rebuilt
scvi.model.SCVI.setup_anndata(
    gdT_mdata,
    layer="counts",
    batch_key=batchkey,
    categorical_covariate_keys=categorical_covariate_keys,
    labels_key="level2_gdT"
)
level2_model_gdT = scvi.model.SCANVI.from_scvi_model(scvi_model_gdT, "Unknown", labels_key="level2_gdT")
level2_model_gdT.train(max_epochs=25, train_size=1.0, validation_size=0) # Add last parameters to stop batchNorm error - 1 value per channel when training
level2_model_gdT.save(prefix+"/gdT_scanvi_level2_model/") #, save_anndata=True)
#level2_model_gdT = scvi.model.SCANVI.load(prefix+"/gdT_scanvi_level2_model/", adata=gdT_mdata)


## Predictions and scores - create output file
## level1
##SCVI_LATENT_KEY = "X_SCVI"
LEVEL1_SCANVI_LATENT_KEY = "level1_X_scANVI"
LEVEL1_SCANVI_PREDICTIONS_KEY = "level1_C_scANVI"

##abT_mdata.obsm[SCVI_LATENT_KEY] = scvi_model.get_latent_representation(abT_mdata)
abT_mdata.obsm[LEVEL1_SCANVI_LATENT_KEY] = level1_model_abT.get_latent_representation(abT_mdata)
abT_mdata.obs[LEVEL1_SCANVI_PREDICTIONS_KEY]= level1_model_abT.predict(abT_mdata)
abT_output_file = pd.DataFrame(abT_mdata.obs[LEVEL1_SCANVI_PREDICTIONS_KEY], index = abT_mdata.obs.index)
level1_abT_latent_df = pd.DataFrame(abT_mdata.obsm[LEVEL1_SCANVI_LATENT_KEY], index = abT_mdata.obs.index)
level1_abT_latent_df.to_csv(prefix+"/latent_level1_abT.csv", index=True)

## Get posterior probabilities for all labels
level1_probs = level1_model_abT.predict(abT_mdata, soft=True)
## Get max probability per cell (i.e., model confidence)
level1_confidence = level1_probs.max(axis=1)
## Add to AnnData
abT_output_file["level1_scanvi_confidence"] = level1_confidence
##abT_output_file.to_csv(prefix_SCANVI+"/predicted_celltypes.csv")

## Add final annotation  with unclear below a threshold
confidence_threshold = 0.85
abT_output_file["level1_final"] = abT_output_file[LEVEL1_SCANVI_PREDICTIONS_KEY]
abT_output_file.loc[abT_output_file["level1_scanvi_confidence"] < confidence_threshold, "level1_final"] = "not classified"

## level2
##SCVI_LATENT_KEY = "X_SCVI"
LEVEL2_SCANVI_LATENT_KEY = "level2_X_scANVI"
LEVEL2_SCANVI_PREDICTIONS_KEY = "level2_C_scANVI"

##abT_mdata.obsm[SCVI_LATENT_KEY] = scvi_model.get_latent_representation(abT_mdata)
abT_mdata.obsm[LEVEL2_SCANVI_LATENT_KEY] = level2_model_abT.get_latent_representation(abT_mdata)
abT_mdata.obs[LEVEL2_SCANVI_PREDICTIONS_KEY]= level2_model_abT.predict(abT_mdata)
abT_output_file[LEVEL2_SCANVI_PREDICTIONS_KEY] = abT_mdata.obs[LEVEL2_SCANVI_PREDICTIONS_KEY]
level2_abT_latent_df = pd.DataFrame(abT_mdata.obsm[LEVEL2_SCANVI_LATENT_KEY], index = abT_mdata.obs.index)
level2_abT_latent_df.to_csv(prefix+"/latent_level2_abT.csv", index=True)


## Get posterior probabilities for all labels
level2_probs = level2_model_abT.predict(abT_mdata, soft=True)
## Get max probability per cell (i.e., model confidence)
level2_confidence = level2_probs.max(axis=1)
## Add to AnnData
abT_output_file["level2_scanvi_confidence"] = level2_confidence
##abT_output_file.to_csv(prefix_SCANVI+"/predicted_celltypes.csv")

## Add final annotation  with unclear below a threshold
abT_output_file["level2_final"] = abT_output_file[LEVEL2_SCANVI_PREDICTIONS_KEY]
abT_output_file.loc[abT_output_file["level2_scanvi_confidence"] < confidence_threshold, "level2_final"] = "not classified"

## Add gdT level1 and level2 to output file
## Get posterior probabilities for all labels
#gdT_probs = gdT_model.predict(gdT_mdata, soft=True)
## Get max probability per cell (i.e., model confidence)
#gdT_confidence = gdT_probs.max(axis=1)
## Add to AnnData
#output_file.loc[output_file.index.isin(gdT_mdata.obs_names), "level2_scanvi_confidence"] = gdT_confidence
#output_file.loc[output_file.index.isin(gdT_mdata.obs_names), "level2_final"] = gdT_model.predict(gdT_mdata)
#output_file.loc[gdT_mdata.obs_names, "level2_scanvi_confidence"] = gdT_confidence
#output_file.loc[gdT_mdata.obs_names, LEVEL2_SCANVI_PREDICTIONS_KEY] = gdT_model.predict(gdT_mdata)
#output_file.loc[gdT_mdata.obs_names, "level2_final"] = output_file.loc[gdT_mdata.obs_names, LEVEL2_SCANVI_PREDICTIONS_KEY]


## gdT output file
## level1
##SCVI_LATENT_KEY = "X_SCVI"
LEVEL1_SCANVI_LATENT_KEY = "level1_X_scANVI"
LEVEL1_SCANVI_PREDICTIONS_KEY = "level1_C_scANVI"

##gdT_mdata.obsm[SCVI_LATENT_KEY] = scvi_model.get_latent_representation(gdT_mdata)
#gdT_mdata.obsm[LEVEL1_SCANVI_LATENT_KEY] = gdT_level1_model.get_latent_representation(gdT_mdata)
gdT_mdata.obs[LEVEL1_SCANVI_PREDICTIONS_KEY]= "gdT"
gdT_output_file = pd.DataFrame(gdT_mdata.obs[LEVEL1_SCANVI_PREDICTIONS_KEY], index = gdT_mdata.obs.index)
#gdT_level1_latent_df = pd.DataFrame(gdT_mdata.obsm[LEVEL1_SCANVI_LATENT_KEY], index = gdT_mdata.obs.index)
#gdT_level1_latent_df.to_csv(prefix+"/latent_gdT_level1.csv", index=True)

## Get posterior probabilities for all labels
#gdT_level1_probs = gdT_level1_model.predict(gdT_mdata, soft=True)
## Get max probability per cell (i.e., model confidence)
#gdT_level1_confidence = gdT_level1_probs.max(axis=1)
## Add to AnnData
gdT_output_file["level1_scanvi_confidence"] = 1
##gdT_output_file.to_csv(prefix_SCANVI+"/predicted_celltypes.csv")

## Add final annotation  with unclear below a threshold
confidence_threshold = 0.85
gdT_output_file["level1_final"] = gdT_output_file[LEVEL1_SCANVI_PREDICTIONS_KEY]
gdT_output_file.loc[gdT_output_file["level1_scanvi_confidence"] < confidence_threshold, "level1_final"] = "not classified"

## level2
##SCVI_LATENT_KEY = "X_SCVI"
LEVEL2_SCANVI_LATENT_KEY = "level2_X_scANVI"
LEVEL2_SCANVI_PREDICTIONS_KEY = "level2_C_scANVI"

##gdT_mdata.obsm[SCVI_LATENT_KEY] = scvi_model.get_latent_representation(gdT_mdata)
gdT_mdata.obsm[LEVEL2_SCANVI_LATENT_KEY] = level2_model_gdT.get_latent_representation(gdT_mdata)
gdT_mdata.obs[LEVEL2_SCANVI_PREDICTIONS_KEY]= level2_model_gdT.predict(gdT_mdata)
gdT_output_file[LEVEL2_SCANVI_PREDICTIONS_KEY] = gdT_mdata.obs[LEVEL2_SCANVI_PREDICTIONS_KEY]
level2_gdT_latent_df = pd.DataFrame(gdT_mdata.obsm[LEVEL2_SCANVI_LATENT_KEY], index = gdT_mdata.obs.index)
level2_gdT_latent_df.to_csv(prefix+"/latent_level2_gdT.csv", index=True)


## Get posterior probabilities for all labels
gdT_level2_probs = level2_model_gdT.predict(gdT_mdata, soft=True)
## Get max probability per cell (i.e., model confidence)
gdT_level2_confidence = gdT_level2_probs.max(axis=1)
## Add to AnnData
gdT_output_file["level2_scanvi_confidence"] = gdT_level2_confidence
##gdT_output_file.to_csv(prefix_SCANVI+"/predicted_celltypes.csv")

## Add final annotation  with unclear below a threshold
gdT_output_file["level2_final"] = gdT_output_file[LEVEL2_SCANVI_PREDICTIONS_KEY]
gdT_output_file.loc[gdT_output_file["level2_scanvi_confidence"] < confidence_threshold, "level2_final"] = "not classified"

## combine two output files
output_file = pd.concat([abT_output_file, gdT_output_file], axis=0)

groups = sorted(output_file["level2_C_scANVI"].unique())
groups = [g for g in groups if g not in ["unclear", "nonT", "remove"]]
  
for annot in groups:
  mask1 = (output_file["level1_final"] == "not classified") & (output_file["level2_final"] == annot)
  output_file.loc[mask1, "level1_final"] = annot.split("_")[0]

## Save user annotations to return to the user
##output_file = mdata.obs[['level1_transfer_labels','level2_transfer_labels']] #,'level2.group_transfer_labels']]
##output_file.index = mdata.obs.index.copy()
##output_file.to_csv(prefix+"/output_annotations.csv", index=True)
output_file.to_csv(prefix+"/predictions_output_file.csv", index=True)
user_output_file = output_file.loc[query_mask, :]
user_output_file.index = mdata.obs.loc[query_mask, :].index
##user_output_file = user_output_file[~user_output_file.index.str.contains("IGT")].copy()
user_output_file.to_csv(prefix+"/user_predictions_output_file.csv", index=True)

##output_file = pd.read_csv(prefix+"/predictions_output_file.csv", index_col=0)
##latent_df = pd.read_csv(prefix+"/latent_level1.csv", index_col=0)
##user_output_file = pd.read_csv(prefix+"/user_output_file.csv", index_col=0)
