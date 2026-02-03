#!/usr/bin/env python3
"""This script run TOTALVI from a Mudata object specifying batch_key and categorical_covariate_keys if needed, """
#author: David Zemmour
#date: 10/08/2024
#run_totalvi_v2.py [cwd] [path to mudata .h5mu] [prefix] [batchkey] [categorical_covariate_keys] [corrected_counts] [denoised_data]

import warnings; warnings.simplefilter('ignore')
import argparse
parser = argparse.ArgumentParser(description="run TOTALVI from a MuData object specifying batch_key and categorical_covariate_keys if needed, can also return corrected counts and denoised data")
parser.add_argument('--working_dir', help='Working directory')
parser.add_argument('--path_to_anndata', help='Path to combined AnnData object')
parser.add_argument('--metadata_query', help='Path to metadata query')
parser.add_argument('--prefix', default='myprefix', 
                    help='Prefix for the output files (default: myprefix)')
parser.add_argument('--batchkey', default=None, help='Batch key for analysis')
parser.add_argument('--categorical_covariate_keys', default=None, help='categorical_covariate_keys variables (default: None)')
parser.add_argument('--corrected_counts', default=False, help='Returns corrected counts, aka posterior_predictive_sample() (default: False)')
parser.add_argument('--denoised_data', default=False, help='Returns denoised data, aka get_normalized_expression()  (default: False)')
#parser.add_argument('--latent_key', help='Key for latent space')

print("Arguments")
args = parser.parse_args()
working_dir = args.working_dir
path_to_anndata = args.path_to_anndata
metadata_query_path = args.metadata_query
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
#print(f"Path to ImmgenT AnnData: {path_to_ImmgenT}")
#print(f"Path to query AnnData: {path_to_query}")
print(f"Path to combined AnnData: {path_to_anndata}")
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
from scipy.sparse import csc_matrix
import matplotlib.pyplot as plt
import pickle
import pandas as pd
import seaborn as sns
import torch
import fast_matrix_market as fmm
#import pymde #to run MDE
#from annoy import AnnoyIndex
#from collections import Counter

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
#mdata = mu.read(path_to_mudata) #totalvi_igt1_56_allgenes_Treg_20240306/adata.h5mu") #mu.read("totalvi_igt1_56_allgenes_Treg_20240306/adata.h5mu
#print(mdata)
mdata = AnnData.read(path_to_anndata) #load this one (results of this block)
#mdata = mdata[~mdata.obs_names.str.contains("IGT")].copy()
All_IGT_values = [f'IGT{i}' for i in range(1, 105)]
mdata = mdata[~mdata.obs['IGT'].isin(All_IGT_values)].copy()
mdata.raw = None
print(mdata)

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
#annotation_level2 = pd.read_csv("/n/groups/cbdm_lab/odc180/ImmgenT_workshop/ImmgenT_freeze_20250109/igt1_104_withtotalvi20250505_annotation_table.csv", delimiter = ",", index_col=0)
#annotation_level2 = annotation_level2.loc[annotation_level2.index.intersection(mdata.obs.index)]
#mdata.obs['level1'] = "Unknown"
#mdata.obs['level2'] = "Unknown"
##mdata.obs['level2.group'] = "Unknown"
#mdata.obs.loc[annotation_level2.index, 'level1'] = annotation_level2['level1'].values
#mdata.obs.loc[annotation_level2.index, 'level2'] = annotation_level2['level2'].values
##mdata.obs.loc[annotation_level2.index, 'level2.group'] = annotation_level2['level2.group'].values

## add final classifications from medata query
metadata_query = pd.read_csv(metadata_query_path, delimiter = ",", index_col=0)
print("metadata_query shape pre filter: ",  metadata_query.shape)
print("mdata pre filter: ")
print(mdata)
common_idx = mdata.obs_names.intersection(metadata_query.index)
mdata = mdata[common_idx, :].copy()
metadata_query = metadata_query.loc[common_idx]
print("metadata_query shape post filter: ",  metadata_query.shape)
print("mdata post filter: ")
print(mdata)
mdata.obs.loc[metadata_query.index, 'level1_C_scANVI'] = metadata_query.loc[metadata_query.index, 'level1_C_scANVI'].values
mdata.obs.loc[metadata_query.index, 'level1_scanvi_confidence'] = metadata_query.loc[metadata_query.index, 'level1_scanvi_confidence'].values
mdata.obs.loc[metadata_query.index, 'level1_final'] = metadata_query.loc[metadata_query.index, 'level1_final'].values
mdata.obs.loc[metadata_query.index, 'level2_C_scANVI'] = metadata_query.loc[metadata_query.index, 'level2_C_scANVI'].values
mdata.obs.loc[metadata_query.index, 'level2_scanvi_confidence'] = metadata_query.loc[metadata_query.index, 'level2_scanvi_confidence'].values
mdata.obs.loc[metadata_query.index, 'level2_final'] = metadata_query.loc[metadata_query.index, 'level2_final'].values

## use original output file

## add not classified percentage

## start for loop

## not just level1 and level2


## check percentages
#...
## < 1% not classified stop
## if > 1% and if % not classified is not prev % /2 (not 50% improvement), then stop
# if < 1% or no 50% improvement, stop
# if > 1% and 50% improvement, continue
initial_percentage = (mdata.obs['level2_final'].eq('not classified').sum() / len(mdata.obs)) * 100
original_percentage = 100
initial_difference = (original_percentage/initial_percentage)

print("original percentage and difference: ")
print(initial_percentage)
print(initial_difference)

if (initial_percentage < 1 or initial_difference < 1.5):
    output_file = metadata_query

percentage = initial_percentage
difference = 100

## Separate gdT and abT mdata
#abT_mdata = mdata[mdata.obs['level1_final'].astype(str) != "gdT" , :].copy()
#gdT_mdata = mdata[mdata.obs['level1_final'] == "gdT" ,:].copy()

i = 1
while ((percentage > 1) and (difference > 1.5)):
  
  ## Separate gdT and abT mdata
  abT_mdata = mdata[mdata.obs['level1_final'].astype(str) != "gdT" , :].copy()
  gdT_mdata = mdata[mdata.obs['level1_final'].astype(str) == "gdT" ,:].copy()

  ## Train separate SCVI models for gdT and abT
  ## abT
  print("loop start: ", i)
  print(percentage)
  print(difference)
  
  #try:
  #    scvi_model_abT
  #except NameError:
  # model does not exist — run your code
  # Comment if already trained
  scvi.model.SCVI.setup_anndata(abT_mdata, layer = "counts", batch_key = batchkey, categorical_covariate_keys = categorical_covariate_keys)
  scvi_model_abT = scvi.model.SCVI(abT_mdata, n_latent=30, n_layers=2)
  scvi_model_abT.train(100)
  #scvi_model_abT.save(prefix+"/scvi_model_abT_not_classified/", save_anndata=True)
  #scvi_model_abT = scvi.model.SCVI.load(prefix+"/scvi_model_abT_not_classified/", adata=abT_mdata)

  # Comment if already trained
  #scvi.model.SCVI.setup_anndata(abT_mdata, layer = "counts", batch_key = batchkey, categorical_covariate_keys = categorical_covariate_keys)
  #scvi_model_abT = scvi.model.SCVI(abT_mdata, n_latent=30, n_layers=2)
  #scvi_model_abT.train()
  #scvi_model_abT.save(prefix+"/scvi_model_abT/", save_anndata=True)
  #scvi_model_abT = scvi.model.SCVI.load(prefix+"/scvi_model_abT/", adata=abT_mdata)
  
  SCANVI_LATENT_KEY = "X_scANVI"
  print("Save latent_representation.csv")
  #latent_representation = scvi_model_abT.get_latent_representation()
  #SCVI_LATENT_KEY = "X_scVI"
  ##mdata.mod['RNA'].obsm[SCVI_LATENT_KEY] = latent_representation
  #abT_mdata.obsm[SCVI_LATENT_KEY] = latent_representation
  ##latent_df = pd.DataFrame(latent_representation, index = mdata.mod['RNA'].obs.index)
  #latent_df_scvi_abT = pd.DataFrame(latent_representation, index = abT_mdata.obs.index)
  
  #latent_df_scvi.to_csv(prefix+"/latent_scvi.csv", index=True)
  
  print("Save umap.csv")
  #sc.pp.neighbors(abT_mdata, use_rep=SCVI_LATENT_KEY)
  #sc.tl.umap(abT_mdata, min_dist=0.4)
  ##umap_df = pd.DataFrame(mdata.obsm['X_umap'], index = mdata.mod['RNA'].obs.index)
  #umap_df = pd.DataFrame(abT_mdata.obsm['X_umap'], index = abT_mdata.obs.index)
  #umap_df.to_csv(prefix+"/umap_python.csv", index=True)
  
  #try:
  #    scvi_model_gdT
  #except NameError:
  # model does not exist — run your code
  # Comment if already trained
  if gdT_mdata.obs.shape[0] > 1 :
      scvi.model.SCVI.setup_anndata(gdT_mdata, layer = "counts", batch_key = batchkey, categorical_covariate_keys = categorical_covariate_keys)
      scvi_model_gdT = scvi.model.SCVI(gdT_mdata, n_latent=30, n_layers=2)
      scvi_model_gdT.train(100)
      #scvi_model_gdT.save(prefix+"/scvi_model_gdT_not_classified/", save_anndata=True)
      #scvi_model_gdT = scvi.model.SCVI.load(prefix+"/scvi_model_gdT_not_classified/", adata=gdT_mdata)


  ## gdT
  # Comment if already trained
  #scvi.model.SCVI.setup_anndata(gdT_mdata, layer = "counts", batch_key = batchkey, categorical_covariate_keys = categorical_covariate_keys)
  #scvi_model_gdT = scvi.model.SCVI(gdT_mdata, n_latent=30, n_layers=2)
  #scvi_model_gdT.train()
  #scvi_model_gdT.save(prefix+"/scvi_model_gdT/", save_anndata=True)
  #scvi_model_gdT = scvi.model.SCVI.load(prefix+"/scvi_model_gdT/", adata=gdT_mdata)
  
  SCANVI_LATENT_KEY = "X_scANVI"
  print("Save latent_representation.csv")
  #latent_representation = scvi_model_gdT.get_latent_representation()
  #SCVI_LATENT_KEY = "X_scVI"
  ##mdata.mod['RNA'].obsm[SCVI_LATENT_KEY] = latent_representation
  #gdT_mdata.obsm[SCVI_LATENT_KEY] = latent_representation
  ##latent_df = pd.DataFrame(latent_representation, index = mdata.mod['RNA'].obs.index)
  #latent_df_scvi_gdT = pd.DataFrame(latent_representation, index = gdT_mdata.obs.index)
  
  #latent_df_scvi.to_csv(prefix+"/latent_scvi.csv", index=True)
  
  print("Save umap.csv")
  #sc.pp.neighbors(gdT_mdata, use_rep=SCVI_LATENT_KEY)
  #sc.tl.umap(gdT_mdata, min_dist=0.4)
  ##umap_df = pd.DataFrame(mdata.obsm['X_umap'], index = mdata.mod['RNA'].obs.index)
  #umap_df = pd.DataFrame(gdT_mdata.obsm['X_umap'], index = gdT_mdata.obs.index)
  #umap_df.to_csv(prefix+"/umap_python.csv", index=True)
  
  print("Save Anndata")
  #mdata.write_h5ad(prefix+"/adata_RNA.h5ad")
  
  ##mdata = AnnData.read(prefix+"_orig/adata_RNA.h5ad")
  ##latent_df = pd.read_csv(prefix+"_orig/latent.csv", index_col=0)
  
  print("Predict level1 and 2 annotations with SCANVI")
  ## Define ref and query masks (future use)
  #ref_mask = mdata.obs["origin"] == "ImmgenT"
  #query_mask = mdata.obs["origin"] == "query"
  
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
  abT_mdata.obs["level1_final"] = abT_mdata.obs["level1_final"].astype(str)
  abT_labels_to_keep = abT_mdata.obs["level1_final"].unique().astype(str)
  abT_labels_to_keep = [x for x in abT_labels_to_keep if x != "not classified"]
  abT_mdata.obs["level1_abT"] = abT_mdata.obs["level1_final"].where(
      abT_mdata.obs["level1_final"].isin(abT_labels_to_keep),
      other="not classified"
  )
  
  ## Fix batch norm error - cant have a predicting label with < cells - 5 cells is the more stable count
  label_counts = abT_mdata.obs["level1_abT"].value_counts()
  drop = label_counts[label_counts < 5].index.tolist()     # labels with < 5 cells are considered drop
  abT_mdata.obs.loc[abT_mdata.obs["level1_abT"].isin(drop), "level1_abT"] = "not classified"
  
  abT_mdata.obs["level1_abT"] = abT_mdata.obs["level1_abT"].astype("category")
  # Setup anndata fresh so the label encoder is rebuilt
  scvi.model.SCVI.setup_anndata(
      abT_mdata,
      layer="counts",
      batch_key=batchkey,
      categorical_covariate_keys=categorical_covariate_keys,
      labels_key="level1_abT"
  )

  n_labeled = (abT_mdata.obs["level1_abT"] != "not classified").sum()
  print("Labeled abT cells:", n_labeled)
  if n_labeled < 2:
      print("Skipping abT SCANVI: <2 labeled cells after filtering")
  else:
      level1_model_abT = scvi.model.SCANVI.from_scvi_model(scvi_model_abT, "not classified", labels_key="level1_abT")
      level1_model_abT.train(30, train_size=0.9,batch_size=128, validation_size=0.1)
      #level1_model_abT.save(prefix+"/abT_scanvi_level1_model/", save_anndata=True)
  
  #level1_model_abT = scvi.model.SCANVI.from_scvi_model(scvi_model_abT, "not classified", labels_key="level1_abT")
  #level1_model_abT.train(30, train_size=1.0, validation_size=0)
  #level1_model_abT.save(prefix+"/abT_scanvi_level1_model/", save_anndata=True)
  #level1_model_abT = scvi.model.SCANVI.load(prefix+"/abT_scanvi_level1_model/", adata=abT_mdata)
  
  
  ## level2
  ## Training scanvi model on scvi model
  abT_mdata.obs["level2_final"] = abT_mdata.obs["level2_final"].astype(str)
  abT_labels_to_keep = abT_mdata.obs["level2_final"].unique().astype(str)
  abT_labels_to_keep = [x for x in abT_labels_to_keep if x != "not classified"]
  abT_mdata.obs["level2_abT"] = abT_mdata.obs["level2_final"].where(
      abT_mdata.obs["level2_final"].isin(abT_labels_to_keep),
      other="not classified"
  )
  
  label_counts = abT_mdata.obs["level2_abT"].value_counts()
  drop = label_counts[label_counts < 5].index.tolist()     # labels with < 5 cells are considered drop
  abT_mdata.obs.loc[abT_mdata.obs["level2_abT"].isin(drop), "level2_abT"] = "not classified"

  abT_mdata.obs["level2_abT"] = abT_mdata.obs["level2_abT"].astype("category")
  # Setup anndata fresh so the label encoder is rebuilt
  scvi.model.SCVI.setup_anndata(
      abT_mdata,
      layer="counts",
      batch_key=batchkey,
      categorical_covariate_keys=categorical_covariate_keys,
      labels_key="level2_abT"
  )
  level2_model_abT = scvi.model.SCANVI.from_scvi_model(scvi_model_abT, "not classified", labels_key="level2_abT")
  level2_model_abT.train(30, train_size=1.0, validation_size=0)
  #level2_model_abT.save(prefix+"/abT_scanvi_level2_model/", save_anndata=True)
  #level2_model_abT = scvi.model.SCANVI.load(prefix+"/abT_scanvi_level2_model/", adata=abT_mdata)
  

  if gdT_mdata.obs.shape[0] > 1 :
      ## gdT only
      #gdT_mdata = mdata[mdata.obs['level1'] == "gdT" ,:].copy()
      print("gdT mdata:")
      print(gdT_mdata)
      ## level1
      ## Training scanvi model on scvi model
      gdT_mdata.obs["level2_final"] = gdT_mdata.obs["level2_final"].astype(str)
      gdT_labels_to_keep = gdT_mdata.obs["level2_final"].unique().astype(str)
      gdT_labels_to_keep = [x for x in gdT_labels_to_keep if x != "not classified"]
      gdT_mdata.obs["level2_gdT"] = gdT_mdata.obs["level2_final"].where(
              gdT_mdata.obs["level2_final"].isin(gdT_labels_to_keep),
              other="not classified"
              )

      label_counts = gdT_mdata.obs["level2_gdT"].value_counts()
      drop = label_counts[label_counts < 5].index.tolist()     # labels with < 5 cells are considered drop
      gdT_mdata.obs.loc[gdT_mdata.obs["level2_gdT"].isin(drop), "level2_gdT"] = "not classified"

      gdT_mdata.obs["level2_gdT"] = gdT_mdata.obs["level2_gdT"].astype("category")
      # Setup anndata fresh so the label encoder is rebuilt
      scvi.model.SCVI.setup_anndata(
              gdT_mdata,
              layer="counts",
              batch_key=batchkey,
              categorical_covariate_keys=categorical_covariate_keys,
              labels_key="level2_gdT"
              )

      n_labeled = (gdT_mdata.obs["level2_gdT"] != "not classified").sum()
      if n_labeled < 2:
          print("Skipping gdT SCANVI: <2 labeled cells")
      else:
          level2_model_gdT = scvi.model.SCANVI.from_scvi_model(scvi_model_gdT, "not classified",labels_key="level2_gdT")
          level2_model_gdT.train(30, train_size=1.0, validation_size=0)
      #level2_model_gdT = scvi.model.SCANVI.from_scvi_model(scvi_model_gdT, "not classified", labels_key="level2_gdT")
      #level2_model_gdT.train(30, train_size=1.0, validation_size=0)
      #level2_model_gdT.save(prefix+"/gdT_scanvi_level2_model/", save_anndata=True)
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
  #level1_abT_latent_df = pd.DataFrame(abT_mdata.obsm[LEVEL1_SCANVI_LATENT_KEY], index = abT_mdata.obs.index)
  #level1_abT_latent_df.to_csv(prefix+"/latent_level1_abT.csv", index=True)
  
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
  #level2_abT_latent_df = pd.DataFrame(abT_mdata.obsm[LEVEL2_SCANVI_LATENT_KEY], index = abT_mdata.obs.index)
  #level2_abT_latent_df.to_csv(prefix+"/latent_level2_abT.csv", index=True)
  
  
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
  
  
  if gdT_mdata.obs.shape[0] > 1 :
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

      #try level2_model_gdT:
      if level2_model_gdT is not None:
          ##gdT_mdata.obsm[SCVI_LATENT_KEY] = scvi_model.get_latent_representation(gdT_mdata)
          gdT_mdata.obsm[LEVEL2_SCANVI_LATENT_KEY] = level2_model_gdT.get_latent_representation(gdT_mdata)
          gdT_mdata.obs[LEVEL2_SCANVI_PREDICTIONS_KEY]= level2_model_gdT.predict(gdT_mdata)
          gdT_output_file[LEVEL2_SCANVI_PREDICTIONS_KEY] = gdT_mdata.obs[LEVEL2_SCANVI_PREDICTIONS_KEY]
          #level2_gdT_latent_df = pd.DataFrame(gdT_mdata.obsm[LEVEL2_SCANVI_LATENT_KEY], index = gdT_mdata.obs.index)
          #level2_gdT_latent_df.to_csv(prefix+"/latent_level2_gdT.csv", index=True)

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

      #except NameError:
      else:
          gdT_output_file[LEVEL2_SCANVI_PREDICTIONS_KEY] = "not classified"
          gdT_output_file["level2_scanvi_confidence"] = 0
          gdT_output_file["level2_final"] = gdT_output_file[LEVEL2_SCANVI_PREDICTIONS_KEY]
          #gdT_output_file.loc[gdT_output_file["level2_scanvi_confidence"] < confidence_threshold, "level2_final"] = "not classified"

      ## combine two output files
      output_file = pd.concat([abT_output_file, gdT_output_file], axis=0)

  else:
      output_file = abT_output_file

  ## Set not classified based on level1 and level2 classifications
  ## Get unique level2_C_scANVI groups, excluding some
  groups2 = sorted(output_file["level2_C_scANVI"].unique())
  groups2 = [g for g in groups2 if g not in ["unclear", "nonT", "remove"]]
  #groups2 = [g for g in groups2 if g not in ["unclear"]]

  
  #for annot2 in groups2:
    #mask2 = (output_file["level1_final"] == "not classified") & (output_file["level2_final"] == annot)
    #output_file.loc[mask1, "level1_final"] = annot.split("_")[0]

   # mask2 = (output_file["level2_final"] == annot2)
   # output_file.loc[mask2, "level1_final"] = annot2.split("_")[0]

  #mdata.obs.loc[output_file.index, ["level1_final", "level2_final"]] = output_file.loc[output_file.index, ["level1_final", "level2_final"]]
  
  percentage_new = (output_file["level2_final"].eq("not classified").sum() / len(output_file.index)) * 100
  difference = (percentage/percentage_new)
  percentage = percentage_new


  ## Check if model has updated previously known labels, or correctly updated not classified labels
  ## Only update those who have been updated properly
  # Only apply predictions that reduce NC or are high-confidence non-NC
  overlap = output_file.index.intersection(mdata.obs.index)

  ## candidates that convert NC to class
  ## level1
  convert_level1 = overlap[
      (mdata.obs.loc[overlap,"level1_final"]=="not classified") &
      (output_file.loc[overlap,"level1_final"]!="not classified") &
      (output_file.loc[overlap, "level1_scanvi_confidence"] > 0.85)
  ]

  ## Only overwrite accepted updates
  mdata.obs.loc[convert_level1, ["level1_C_scANVI", "level1_scanvi_confidence", "level1_final"]] = output_file.loc[convert_level1, ["level1_C_scANVI", "level1_scanvi_confidence", "level1_final"]].values


  ##level2
  convert_level2 = overlap[
      (mdata.obs.loc[overlap,"level2_final"]=="not classified") &
      (output_file.loc[overlap,"level2_final"]!="not classified") &
      (output_file.loc[overlap, "level2_scanvi_confidence"] > 0.85)
  ]

  ## Only overwrite accepted updates
  mdata.obs.loc[convert_level2, ["level2_C_scANVI", "level2_scanvi_confidence", "level2_final"]] = output_file.loc[convert_level2, ["level2_C_scANVI", "level2_scanvi_confidence", "level2_final"]].values

  output_file = mdata.obs[["level1_C_scANVI", "level1_scanvi_confidence", "level1_final", "level2_C_scANVI", "level2_scanvi_confidence", "level2_final"]]

  #mdata.obs[["level1_final", "level2_final"]] = output_file[["level1_final", "level2_final"]]
  #mdata.obs.loc[output_file.index, ["level1_final", "level1_scanvi_confidence", "level2_final", "level2_scanvi_confidence"]] = output_file.loc[output_file.index, ["level1_final", "level1_scanvi_confidence", "level2_final", "level2_scanvi_confidence"]]

  groups2 = sorted(output_file["level2_C_scANVI"].unique())
  groups2 = [g for g in groups2 if g not in ["unclear", "nonT", "remove"]]
  
  for annot2 in groups2:
    #mask2 = (output_file["level1_final"] == "not classified") & (output_file["level2_final"] == annot)
    #output_file.loc[mask1, "level1_final"] = annot.split("_")[0]

    mask2 = (output_file["level2_final"] == annot2)
    output_file.loc[mask2, "level1_final"] = annot2.split("_")[0]


  print("loop end: ", i)
  print(percentage)
  print(difference)

  i = i + 1
  
  del gdT_mdata
  del abT_mdata


if level1_model_abT is not None:
    level1_model_abT.save(prefix+"/abT_scanvi_level1_model_not_classified/", save_anndata=False)

if level2_model_abT is not None:
    level2_model_abT.save(prefix+"/abT_scanvi_level2_model_not_classified/", save_anndata=False)

gdT_mdata = mdata[mdata.obs['level1_final'] == "gdT" ,:].copy()
if level2_model_gdT is not None:
    level2_model_gdT.save(prefix+"/gdT_scanvi_level2_model_not_classified/", save_anndata=False)

groups1 = sorted(output_file["level1_C_scANVI"].dropna().unique())
#groups1 = sorted(output_file["level1_C_scANVI"].unique())
groups1 = [g for g in groups1 if g not in ["unclear", "nonT", "remove"]]

for annot1 in groups1:
    mask1 = (output_file["level1_final"] == annot1) & (output_file["level2_final"] == "not classified")
    output_file.loc[mask1, "level2_final"] = annot1

## group remove and unclear into nonT and treat nonT as a subgroup - same not classified settings
output_file["level1_final"] = output_file["level1_final"].str.replace(r"(remove|unclear).*", "nonT", case=False, regex=True)
output_file["level2_final"] = output_file["level2_final"].str.replace(r"(remove|unclear).*", "nonT", case=False, regex=True)

mask_nonT1 = (output_file["level1_final"] == "nonT") & (output_file["level2_final"] == "not classified")
output_file.loc[mask_nonT1, "level2_final"] = "nonT"

mask_nonT2 = (output_file["level2_final"] == "nonT")
output_file.loc[mask_nonT2, "level1_final"] = "nonT"


## Save user annotations to return to the user
##output_file = mdata.obs[['level1_transfer_labels','level2_transfer_labels']] #,'level2.group_transfer_labels']]
##output_file.index = mdata.obs.index.copy()
##output_file.to_csv(prefix+"/output_annotations.csv", index=True)
output_file.to_csv(prefix+"/predictions_output_file_not_classified.csv", index=True)
user_output_file = output_file.copy()
##user_output_file.index = query_mdata.obs.index
user_output_file.to_csv(prefix+"/user_predictions_output_file_not_classified.csv", index=True)


