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

import warnings
warnings.filterwarnings("ignore")
# In[]
# --------------------------------------------------------------------
#
# Load the data
#
# --------------------------------------------------------------------
# dataset = "tree3_500_cells110_genes_sigma0.1_1"
dataset = sys.argv[1]
path = f"/project/hli691/scm_paper/1/{dataset}/"
print(path)
try:
    assert os.path.exists(path)
except:
    raise ValueError(f"path: {path} does not exist.")

result_dir = f"{dataset}/"
if not os.path.exists(result_dir):
    os.makedirs(result_dir)

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
regulators = ["gene" + str(int(x)) for x in res[".grn"].rx2("regulators")] 
# region x tf
region2tf = region_tf.rx2(2)
region2gene = res["region_to_gene"]
cell_meta = robjects.conversion.rpy2py(res["cell_meta"])

adata = anndata.AnnData(X = counts_rna.T, obs = cell_meta)
adata.var.index = ["gene" + str(x + 1) for x in range(adata.shape[1])]
base_grn = (region2gene.T @ region2tf > 0).astype(int)
base_grn = pd.DataFrame(data = base_grn, index = adata.var.index.values, columns = regulators)

# In[]
# --------------------------------------------------------------------
#
# Preprocessing
#
# --------------------------------------------------------------------
print("Preprocessing...")
sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)

# In[]
# --------------------------------------------------------------------
#
# Save files
#
# --------------------------------------------------------------------
# 1. run with single cluster 
# counts = pd.DataFrame(data = adata.X, index = adata.obs.index.values, columns = np.array(["gene" + str(x + 1) for x in range(adata.shape[1])])).T
# counts.to_csv(result_dir + "cluster1.table", sep = "\t")

# np.savetxt(result_dir + "regulators.txt", np.array(regulators)[:,None], fmt = "%s")

# base_grn = (region2gene.T @ region2tf > 0).astype(int)
# base_grn = pd.DataFrame(data = base_grn, index = adata.var.index.values, columns = regulators)

# prior = pd.DataFrame(columns = ["regulator", "target", "strength"])
# for target in base_grn.index.values.squeeze():
#     for regulator in base_grn.columns.values.squeeze():
#         if base_grn.loc[target,regulator] != 0:
#             prior = pd.concat([prior, pd.DataFrame(data = [[regulator + "_cluster1", target + "_cluster1", base_grn.loc[target,regulator]]], columns = ["regulator", "target", "strength"])], axis = 0, ignore_index = True)

# prior.to_csv(result_dir + "cluster1_network.txt", sep = "\t", header = False, index = False)

# filelist = pd.DataFrame(data = [["cluster1", result_dir + f"cluster1.table"]])

# filelist.to_csv(result_dir + "filelist.txt", sep = "\t", header = False, index = False)

# In[]
# 2. run with ground truth cluster
counts = pd.DataFrame(data = adata.X, index = adata.obs.index.values, columns = np.array(["gene" + str(x + 1) for x in range(adata.shape[1])])).T
counts_list = []
if dataset[:5] == "tree1":
    counts_cluster1 = counts.loc[:, (adata.obs["depth"] < 0.5).values]
    counts_cluster1.to_csv(result_dir + f"cluster1.table", sep = "\t")
    counts_cluster2 = counts.loc[:, (adata.obs["depth"] >= 0.5).values]
    counts_cluster2.to_csv(result_dir + f"cluster2.table", sep = "\t")
else:
    for id, cluster in enumerate(np.unique(adata.obs["pop"])):
        counts_list.append(counts.loc[:, (adata.obs["pop"] == cluster).values])
        counts_list[id].to_csv(result_dir + f"cluster{id+1}.table", sep = "\t")

np.savetxt(result_dir + "regulators.txt", np.array(regulators)[:,None], fmt = "%s")

base_grn = (region2gene.T @ region2tf > 0).astype(int)
base_grn = pd.DataFrame(data = base_grn, index = adata.var.index.values, columns = regulators)

prior = pd.DataFrame(columns = ["regulator", "target", "strength"])
for target in base_grn.index.values.squeeze():
    for regulator in base_grn.columns.values.squeeze():
        if base_grn.loc[target,regulator] != 0:
            prior = pd.concat([prior, pd.DataFrame(data = [[regulator + "_cluster1", target + "_cluster1", base_grn.loc[target,regulator]]], columns = ["regulator", "target", "strength"])], axis = 0, ignore_index = True)

if dataset[:5] == "tree1":
    prior.to_csv(result_dir + f"cluster1_network.txt", sep = "\t", header = False, index = False)
    prior.to_csv(result_dir + f"cluster2_network.txt", sep = "\t", header = False, index = False)
    filelist = pd.DataFrame(data = [["cluster1", result_dir + f"cluster1.table"], ["cluster2", result_dir + f"cluster2.table"]])
    celltype_tree_ancestor = pd.DataFrame([["cluster2", "cluster1", "0.2", "0.2"]])
    celltype_tree_ancestor.to_csv(result_dir + "celltype_tree_ancestor.txt", sep = "\t", header = False, index = False)
    filelist.to_csv(result_dir + "filelist.txt", sep = "\t", header = False, index = False)

else:
    filelist = pd.DataFrame(columns = ["cluster", "directory"])
    for id, cluster in enumerate(np.unique(adata.obs["pop"])):
        prior.to_csv(result_dir + f"cluster{id+1}_network.txt", sep = "\t", header = False, index = False)
        filelist = pd.concat([filelist, pd.DataFrame(data = [[f"cluster{id+1}", result_dir + f"cluster{id+1}.table"]], columns = ["cluster", "directory"])], axis = 0)
    filelist.to_csv(result_dir + "filelist.txt", sep = "\t", header = False, index = False)
    
    if dataset[:5] == "tree3":
        celltype_tree_ancestor = pd.DataFrame([["cluster1", "cluster2", "0.2", "0.2"], ["cluster3", "cluster2", "0.2", "0.2"], ["cluster4", "cluster2", "0.2", "0.2"]])
        celltype_tree_ancestor.to_csv(result_dir + "celltype_tree_ancestor.txt", sep = "\t", header = False, index = False)
    elif dataset[:5] == "tree5":
        # give probability 0.5, to make sure that the clusters are independent
        celltype_tree_ancestor = pd.DataFrame([["cluster2", "cluster1", "0.5", "0.5"], ["cluster3", "cluster2", "0.5", "0.5"], ["cluster4", "cluster3", "0.5", "0.5"], ["cluster5", "cluster4", "0.5", "0.5"]])
        celltype_tree_ancestor.to_csv(result_dir + "celltype_tree_ancestor.txt", sep = "\t", header = False, index = False)
        




# %%
