# In[]
import numpy as np 
import pandas as pd
import os, sys

# In[]
# dataset = "tree3_500_cells110_genes_sigma0.1_1"
dataset = sys.argv[1]
path = f"/project/hli691/scm_paper/1/{dataset}/"
try:
    assert os.path.exists(path)
except:
    raise ValueError(f"path: {path} does not exist.")

result_dir = path + "scmtni_results/"
if not os.path.exists(result_dir):
    os.makedirs(result_dir)

# In[]
# ----------------------------------------
#
# After running the method
#
# ----------------------------------------
print("Postprocessing...")
ngenes = eval(dataset.split("_")[2][5:])

grn_clusters = []
grn_clusters.append(pd.read_table(f"results_{dataset}/cluster1/fold0/var_mb_pw_k50.txt", header = None))
grn_clusters.append(pd.read_table(f"results_{dataset}/cluster2/fold0/var_mb_pw_k50.txt", header = None))
    
if (dataset[:5] == "tree3") or (dataset[:5] == "tree5"):
    grn_clusters.append(pd.read_table(f"results_{dataset}/cluster3/fold0/var_mb_pw_k50.txt", header = None))
    grn_clusters.append(pd.read_table(f"results_{dataset}/cluster4/fold0/var_mb_pw_k50.txt", header = None))

    if dataset[:5] == "tree5":
        grn_clusters.append(pd.read_table(f"results_{dataset}/cluster5/fold0/var_mb_pw_k50.txt", header = None))

G_clusters = []

for grn_cluster in grn_clusters:
    # source by target
    G_clusters.append(np.zeros((1, ngenes, ngenes)))
    for row in range(grn_cluster.shape[0]):
        reg = grn_cluster.loc[row,0].split("_")[0]
        targ = grn_cluster.loc[row,1].split("_")[0]
        score = grn_cluster.loc[row,2]
        reg_idx = eval(reg[4:])-1
        targ_idx = eval(targ[4:])-1
        G_clusters[-1][0, reg_idx, targ_idx] += score

G_clusters = np.concatenate(G_clusters, axis = 0)
G = np.mean(G_clusters, axis = 0)
np.save(result_dir + "GRN_scmtni.npy", G)
print("Done.")
# %%
