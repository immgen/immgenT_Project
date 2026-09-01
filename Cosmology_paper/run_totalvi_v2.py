#!/project/zemmour/david/envs/scvi120_20241008/bin/python
"""This script run TOTALVI from a Mudata object specifying batch_key and categorical_covariate_keys if needed, """
#author: David Zemmour
#date: 10/08/2024
#run_totalvi_v2.py [cwd] [path to mudata .h5mu] [prefix] [batchkey] [categorical_covariate_keys] [corrected_counts] [denoised_data]

import warnings; warnings.simplefilter('ignore')
import argparse
parser = argparse.ArgumentParser(description="run TOTALVI from a MuData object specifying batch_key and categorical_covariate_keys if needed, can also return corrected counts and denoised data")
parser.add_argument('--working_dir', help='Working directory')
parser.add_argument('--path_to_mudata', help='Path to MuData object')
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
path_to_mudata = args.path_to_mudata
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
print(f"Path to AnnData: {path_to_mudata}")
print(f"Prefix: {prefix}")
print(f"Batch Key: {batchkey}")
print(f"categorical_covariate_keys: {categorical_covariate_keys}")
print(f"corrected_counts: {corrected_counts}")
print(f"denoised_data: {denoised_data}")
#print(f"Latent Key: {latent_key}")

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

print("Reading mudata")
os.chdir(working_dir)
mdata = mu.read(path_to_mudata) #totalvi_igt1_56_allgenes_Treg_20240306/adata.h5mu") #mu.read("totalvi_igt1_56_allgenes_Treg_20240306/adata.h5mu
print(mdata)

# print("Filter out genes with no data")
# total_counts_per_genes = np.array(mdata.mod['RNA'].X.sum(axis=0)).flatten()
# low_count_genes = total_counts_per_genes < 1
# print(f"Number of genes with no counts: {low_count_genes.sum()}")
# mdata.mod["RNA"] = mdata.mod["RNA"][:, ~low_count_genes]

print("Creating RNA layer counts")
mdata.mod["RNA"].layers["counts"] = mdata.mod["RNA"].X.copy()

#if sparse matrix in protein data: print("Converting  mdata.mod['protein'].X to array")
#mdata.mod['protein'].X = mdata.mod['protein'].X.toarray()

mdata.update()

print("Batch key list")
print(mdata.mod["RNA"].obs[batchkey].unique().tolist())
print("Batch key list: any NA?")
print((mdata.mod["RNA"].obs[batchkey].isna()).any())


print("categorical_covariate_keys")
if categorical_covariate_keys is not None:
    for c in categorical_covariate_keys:
        print(c)
        print(mdata.mod["RNA"].obs[c].unique().tolist())
        print("any NA?")
        print((mdata.mod["RNA"].obs[c].isna()).any())
        #print(mdata.mod["RNA"].obs[c].head(10))
        #mdata.mod["RNA"].obs[c] = mdata.mod["RNA"].obs[c].str.replace('.', '', regex=False)
        #print(mdata.mod["RNA"].obs[c].head(10))

print("Setup mu_data")
scvi.model.TOTALVI.setup_mudata(
    mdata,
    rna_layer="counts",
    categorical_covariate_keys = categorical_covariate_keys, #make sure this is a list!
    protein_layer=None,
    batch_key=batchkey,
    modalities={
        "rna_layer": "RNA",
        "protein_layer": "protein",
        "batch_key": "RNA",
        "categorical_covariate_keys":"RNA"
    },
)

model = scvi.model.TOTALVI(mdata, n_latent = 30, gene_likelihood = "nb") #extra_decoder_kwargs={'activation_function_bg': 'exp'}

print("Train model")
model.train()

print("Save model")
model.save(prefix, save_anndata=True)
# model = scvi.model.TOTALVI.load(prefix)
# mdata = model.adata

fig, ax = plt.subplots(1, 1)
model.history["elbo_train"].plot(ax=ax, label="train")
model.history["elbo_validation"].plot(ax=ax, label="validation")
ax.set(title="Negative ELBO over training epochs", ylim=(1200, 1400))
ax.legend()
fig = plt.gcf() 
fig.savefig(prefix+"/training_elbo_plot.pdf")

print("Save latent_representation.csv")
latent_representation = model.get_latent_representation()
TOTALVI_LATENT_KEY = "X_totalVI"
#mdata.mod['RNA'].obsm[TOTALVI_LATENT_KEY] = latent_representation
mdata.obsm[TOTALVI_LATENT_KEY] = latent_representation
latent_df = pd.DataFrame(latent_representation, index = mdata.mod['RNA'].obs.index)

latent_df.to_csv(prefix+"/latent.csv", index=True)
# latent_df = pd.read_csv(prefix+"/latent.csv", index_col = 0)
# mdata.obsm[TOTALVI_LATENT_KEY] = latent_df
# mdata.obsm[TOTALVI_LATENT_KEY] = mdata.obsm[TOTALVI_LATENT_KEY].values

print("Save umap.csv")
sc.pp.neighbors(mdata, use_rep=TOTALVI_LATENT_KEY)
sc.tl.umap(mdata, min_dist=0.4)
umap_df = pd.DataFrame(mdata.obsm['X_umap'], index = mdata.mod['RNA'].obs.index)
umap_df.to_csv(prefix+"/umap_python.csv", index=True)

print("Save mde.csv")
# mdata.obsm["X_mde"] = scvi.model.utils.mde(mdata.obsm[TOTALVI_LATENT_KEY])
# mde_df = pd.DataFrame(mdata.obsm["X_mde"], index = mdata.mod['RNA'].obs.index) 
# mde_df.to_csv(prefix+"/mde.csv", index=True)

mde = pymde.preserve_neighbors(
    mdata.obsm[TOTALVI_LATENT_KEY], #mdata.obsm[TOTALVI_LATENT_KEY].values
    embedding_dim=2,
    constraint=pymde.Standardized(),
    repulsive_fraction=0.7,
    n_neighbors = 15,
    verbose=True,
    device = 'cuda'
)

import pickle
model_pkl_file = prefix+"/mde_model.pkl"  
with open(model_pkl_file, 'wb') as file:  
    pickle.dump(mde, file)

# with open(model_pkl_file, 'rb') as file:  
#     model = pickle.load(file)

embedding = mde.embed(verbose=True)
mde_df = pd.DataFrame(embedding.cpu(), index = mdata.mod['RNA'].obs.index) 
mde_df.to_csv(prefix+"/mde2.csv", index=True)

pairs, distorsions = mde.high_distortion_pairs()
pairs_df = pd.DataFrame(pairs.cpu(), columns=['cell1', 'cell2'])
pairs_df['cell1'] = mdata.mod['RNA'].obs.index[pairs_df['cell1']]
pairs_df['cell2'] = mdata.mod['RNA'].obs.index[pairs_df['cell2']]
distorsions_df = pd.DataFrame(distorsions.cpu(),columns=['distortion'])
pairs_df['distortion'] = distorsions_df['distortion']
#pairs_df.sample(n=1000000, replace=False).to_csv(prefix+"/mde_distorsions_sample.csv")
pairs_df.to_csv(prefix+"/mde_distorsions.csv")
#most_distorted_pairs = pairs[:500]
mde.distortions_cdf()
fig = plt.gcf() 
fig.savefig(prefix+"/mde_distorsion_cdf.png")

umap_mde_aligned = pymde.align(source=embedding.cpu(), target=mdata.obsm['X_umap'])
umap_mde_aligned_df = pd.DataFrame(umap_mde_aligned, index = mdata.mod['RNA'].obs.index)
umap_mde_aligned_df.to_csv(prefix+"/umap_mde_aligned.csv", index=True)
    
# load model from pickle file
# with open(model_pkl_file, 'rb') as file:  
#     mde2 = pickle.load(file)



#mdata.write(prefix+"/adata.h5mu")

# if corrected_counts:
#     print("Save corrected counts")
#     import numpy as np
#     from scipy.sparse import csr_matrix
#     from scipy.io import mmwrite
#     newcounts = model.posterior_predictive_sample() #RNA AND Proteins!
#     newcounts = newcounts.astype(int)
#     newcounts_sparse = csr_matrix(newcounts)
#     mmwrite(prefix+"/counts_corrected.mtx", newcounts_sparse)
#     genes_prot = pd.DataFrame(mdata.var_names)
#     genes_prot.to_csv(prefix+"/genes.csv", index = False,header=False )
#     cells = pd.DataFrame(mdata.obs_names)
#     cells.to_csv(prefix+"/cells.csv", index = False,header=False )
# 
# if denoised_data:
#     print("Calculate denoised data")
#     denoised = model.get_normalized_expression()
#     mdata.mod['protein'].layers["counts"] = mdata.mod["protein"].X.copy()
#     mdata.mod['protein'].layers["protein_denoised"] = denoised[1]
#     mdata.mod['RNA'].layers["rna_denoised"] = denoised[0]
#     print("Save denoised data")
#     mdata.write(prefix+"/adata_withdenoised.h5mu")
#     #export mu data in separate AnnData for R
#     print("Export data rna_data.h5ad and protein_data.h5ad in AnnData for R")
#     mdata.mod['RNA'].write_h5ad(prefix+"/rna_data.h5ad")
#     mdata.mod['protein'].write_h5ad(prefix+"/protein_data.h5ad")
#     print("Done")

