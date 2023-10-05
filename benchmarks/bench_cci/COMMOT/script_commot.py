# In[]
import commot as ct
import scanpy as sc
import pandas as pd
import numpy as np
from anndata import AnnData
import scanpy as sc
from sklearn.metrics import pairwise_distances
import os, sys

# In[]
dataset = sys.argv[1]
n_spurious = eval(sys.argv[2])
K = eval(sys.argv[3])

data_dir = f"/localscratch/hechen/cci_sc/{dataset}/"
result_dir = data_dir + "result_commot/"
if not os.path.exists(result_dir):
    os.makedirs(result_dir)

# Cell-type-level CCI: ligand-receptor interaction strength for each pair of cell types
cci = pd.read_csv(data_dir + "cci.csv")
# Cell-level CCI: cell * cell * LR pair matrix, binary, indicating cci ground truth
cci_gt = np.load(data_dir + "cci_gt.npy")
# cell type annotation
cell_info = pd.read_csv(data_dir + "cell_info.csv", index_col = 0)
# the gene by cell count matrix
counts = pd.read_csv(data_dir + "counts.csv", index_col = 0)
# the GRN ground truth
grn = pd.read_csv(data_dir + "grn.csv")
# the spatial location of cells
locs = pd.read_csv(data_dir + "locs.csv", index_col = 0)


# In[]
adata = AnnData(X = counts.T.values)
adata.var.index = ["gene" + str(x + 1) for x in range(counts.shape[0])]
adata.obs = cell_info
adata.obsm["spatial"] = locs.values

# preprocessing
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)

# using ground truth cci set: 6 in total
ligrec_db = cci.groupby(["ligand", "receptor"]).size().reset_index().rename(columns={0:'count'})
ligrec_db = ligrec_db[["ligand", "receptor"]]
for idx in ligrec_db.index:
    ligrec_db.loc[idx, "ligand"] = "gene" + str(ligrec_db.loc[idx, "ligand"])
    ligrec_db.loc[idx, "receptor"] = "gene" + str(ligrec_db.loc[idx, "receptor"])
    ligrec_db["ligand_receptor"] = str(ligrec_db.loc[idx, "ligand"]) + "_" + str(ligrec_db.loc[idx, "receptor"])
ligrec_db["pathway"] = None


# TODO: add spurious LR pairs

np.random.seed(0)

if n_spurious > 0:
    for iter in range(n_spurious):
        lig, rec = np.random.randint(low = 1, high = counts.shape[0], size = 2)
        lig = "gene" + str(lig)
        rec = "gene" + str(rec)
        if not lig + "_" + rec in ligrec_db["ligand_receptor"]:
            ligrec_db = pd.concat([ligrec_db, pd.DataFrame(data = [[lig, rec,  lig+"_"+rec, None]], columns = ["ligand", "receptor", "ligand_receptor", "pathway"])], axis = 0, ignore_index = True)


# In[]

D = pairwise_distances(locs.values)
knn_index = np.argpartition(D, kth = K-1, axis = 1)[:, (K-1)]
kth_dist = np.take_along_axis(D, knn_index[:,None], axis = 1)
# maximum knn distance as dist cut-off, change with K selection
dis_thr = np.max(kth_dist)

# In[]
# heterometric: boolean, whether the protein complexes can form ligands/receptors. 
# If heterometric is True, commot will include ligand/receptor complexes of the name "ligand1_ligand2" in df_ligrec, where "_" is the heteromeric_deliminator
# If heterometric is False, commot will only consider the ligand/receptor exists in df_ligrec.
# dis_thr: maximum distance of ligand-receptor interaction, ligand-receptor with distance larger can dis_thr will not be considered as interacting
# pathway_sum: If True, calculate the cell*cell*pathway score matrix.
ct.tl.spatial_communication(adata, database_name='simulated', df_ligrec=ligrec_db, dis_thr=dis_thr, heteromeric=False, pathway_sum = False)

LR_pairs = []
cci_predict = []
for LR_pair in adata.obsp.keys():
    LR_pairs.append(LR_pair)
    cci_predict.append(adata.obsp[LR_pair].toarray()[None,:,:])
LR_pairs = np.array(LR_pairs)
cci_predict = np.concatenate(cci_predict, axis = 0)
print(cci_predict.shape)
print(LR_pairs)
np.savetxt(result_dir + f"LR_index_{n_spurious}_{K}.txt", LR_pairs, fmt = "%s")
np.save(result_dir + f"cci_predict_{n_spurious}_{K}.npy", cci_predict)



# %%
