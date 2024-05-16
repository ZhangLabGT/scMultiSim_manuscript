import pandas as pd
from sklearn.metrics import adjusted_rand_score

ground_truth = pd.read_csv("./stagate/domains.csv", index_col=0).to_numpy().flatten()
pred1 = pd.read_csv("./scHybridNMF/result.csv", header=None).to_numpy().flatten()
pred2 = pd.read_csv("./stagate/stagate_out.csv", index_col=0).to_numpy().flatten()

# calculate ARI
ari1 = adjusted_rand_score(ground_truth, pred1)
ari2 = adjusted_rand_score(ground_truth, pred2)

print(ari1)
print(ari2)
