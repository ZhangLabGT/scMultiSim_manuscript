# In[]

import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
pandas2ri.activate()

import numpy as np 
import anndata
import scanpy as sc
import celloracle as co
import matplotlib.pyplot as plt
import pandas as pd
import os, sys


# In[]
# --------------------------------------------------------------------
#
# Load the data
#
# --------------------------------------------------------------------
# path = "/project/hli691/scm_paper/1/tree1_500_cells110_genes_sigma0.1_1/"
path = sys.argv[1]
print(path)
try:
    assert os.path.exists(path)
except:
    raise ValueError(f"path: {path} does not exist.")

print("Read count matrix and base GRN...")
readRDS = robjects.r['readRDS']
res = readRDS(path + "res.rds")
region_tf = readRDS(path + "region_to_tf.rds")
# count matrix: gene x cell
counts_rna = res["counts"]
# count matrix observed
counts_rna_obs = res["counts_obs"]
# region x cell
counts_atac = res["atacseq_data"]
# gene x tf
grn = res[".grn"].rx2("geff")
regulators = ["gene_" + str(int(x)) for x in res[".grn"].rx2("regulators")] 
# region x tf
region2tf = region_tf.rx2(2)
region2gene = res["region_to_gene"]
cell_meta = robjects.conversion.rpy2py(res["cell_meta"])

# In[]
# --------------------------------------------------------------------
#
# Preprocessing
#
# --------------------------------------------------------------------
base_grn = (region2gene.T @ region2tf > 0).astype(int)
adata = anndata.AnnData(X = counts_rna.T, obs = cell_meta)
adata.var.index = ["gene_" + str(x + 1) for x in range(adata.shape[1])]
base_grn = pd.DataFrame(data = base_grn, index = adata.var.index.values, columns = regulators)

print("Preprocess the count matrix...")
# preprocessing, counts_per_cell_after cannot be too small, or import_anndata_as_raw_count will fail
sc.pp.normalize_per_cell(adata, counts_per_cell_after = 100)
adata.raw = adata
adata.layers["raw_count"] = adata.raw.X.copy()
# log transformation and scaling
sc.pp.log1p(adata)
sc.pp.scale(adata)
# dimension reduction
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=30, n_pcs=30)


# sc.tl.diffmap(adata)
# # Calculate neihbors again based on diffusionmap
# sc.pp.neighbors(adata, n_neighbors=30, use_rep='X_diffmap')
# sc.tl.paga(adata, groups='leiden')
# sc.tl.draw_graph(adata, init_pos='paga', random_state=123)

sc.tl.umap(adata)
# clustering
# sc.tl.leiden(adata, resolution=0.2)
adata.obs["leiden"] = "cluster0"
# visualization
sc.pl.umap(adata, color = "leiden")

# In[]
# --------------------------------------------------------------------
#
# Load count matrix and base GRN into CellOracle
#
# --------------------------------------------------------------------
print("Running cell oracle....")
# Instantiate Oracle object
oracle = co.Oracle()
# In this notebook, we use the unscaled mRNA count for the nput of Oracle object.
adata.X = adata.layers["raw_count"].copy()
# Instantiate Oracle object.
oracle.import_anndata_as_raw_count(adata=adata, cluster_column_name="leiden", embedding_name="X_pca", transform = "log2")

# Creat TFinfo_dictionary, each key is a target gene, the value is the list of TF connecting to the target gene in baseGRN
TF_to_TG_dictionary = {}
for tf_id, tf in enumerate(base_grn.columns.values):
    TF_to_TG_dictionary[tf] = [x for x in base_grn.index.values[base_grn.loc[:, tf] != 0]] 

TG_to_TF_dictionary = co.utility.inverse_dictionary(TF_to_TG_dictionary)

# You can load TF info dataframe with the following code.
oracle.import_TF_data(TFdict=TG_to_TF_dictionary)


# In[]
# --------------------------------------------------------------------
#
# Data imputation
#
# --------------------------------------------------------------------
# knn imputation
# Perform PCA
oracle.perform_PCA()

# Select important PCs
plt.plot(np.cumsum(oracle.pca.explained_variance_ratio_)[:100])
n_comps = np.where(np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_))>0.002))[0][0]
plt.axvline(n_comps, c="k")
plt.show()
print(n_comps)
n_comps = min(n_comps, 50)

n_cell = oracle.adata.shape[0]
k = int(0.025*n_cell)
oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k*8, b_maxl=k*4, n_jobs=4)

# In[]
# --------------------------------------------------------------------
#
# Infer GRNs
#
# --------------------------------------------------------------------
links = oracle.get_links(cluster_name_for_GRN_unit="leiden", alpha=10, verbose_level=10)

# In[]
# --------------------------------------------------------------------
#
# Save results
#
# --------------------------------------------------------------------
print("Save results...")
if not os.path.exists(path + "celloracle_results"):
    os.makedirs(path + "celloracle_results")

# save cluster assignment result
cluster_assign = adata.obs[["leiden"]]
cluster_assign.to_csv(path + "celloracle_results/cluster_assignment.csv")

# save grn inference result
for cluster_id in links.links_dict.keys():
    # direct output of celloracle
    grn_cluster = links.links_dict[cluster_id]
    grn_cluster.to_csv(path + "celloracle_results/grn_df_" + str(cluster_id) + ".csv")

    # the coeff mean can be treated as the edge weight of the grn
    grn_coef_mean = pd.DataFrame(data = 0, index = adata.var.index.values, columns = adata.var.index.values)
    for i in range(grn_cluster.shape[0]):
        source = grn_cluster.loc[i, "source"]
        target = grn_cluster.loc[i, "target"]
        # row source, column target
        grn_coef_mean.loc[source, target] = grn_cluster.loc[i, "coef_mean"]
    # save results
    np.savetxt(fname = path + "celloracle_results/grn_coef_mean_" + str(cluster_id) + ".txt", X = grn_coef_mean.values)


# %%
