from collections import namedtuple
from typing import List
from typing_extensions import assert_type
import scvelo as scv
import pandas as pd
import anndata as ad
import numpy as np
import scanpy
from scipy.sparse import coo_matrix
import loompy
import velocyto

from session import BenchSession

scv.set_figure_params()

N_NEIGHBORS = 30
N_PCS = 30

VeloResult = namedtuple('VeloResult', ['obj', 'velocity', 'expression_future', 'transition_matrix'])


def cosine_sim(velo_pred, velo_true):
    norm_velo_pred = velo_pred/np.linalg.norm(velo_pred, ord = 2, axis = 1)[:,None]
    norm_velo_true = velo_true/np.linalg.norm(velo_true, ord = 2, axis = 1)[:,None]

    cos_score = np.sum(norm_velo_pred * norm_velo_true, axis = 1)
    aver_cos_score = 1/cos_score.shape[0] * np.sum(cos_score)
    return aver_cos_score


def bench_velo(method: str, session: BenchSession,
    n_nbs=N_NEIGHBORS, n_pcs=N_PCS, log_scale=False,
    scv_mode='stochastic', var_names='velocity_genes', layer='spliced', assumption='constant_velocity'):

    if method not in ['scvelo', 'velocyto']:
        raise RuntimeError(f'Method {method} is not available.')
    
    if method == 'scvelo':
        if scv_mode not in ['stochastic', 'deterministic', 'dynamical', 'dynamical_residuals']:
            raise RuntimeError(f'scvelo mode {scv_mode} is not available.')
        if var_names not in ['velocity_genes', 'all']:
            raise RuntimeError(f'scvelo var_names {var_names} is not available.')
        if layer not in ['spliced', 'imputed']:
            raise RuntimeError(f'scvelo {layer} is not available.')
    elif method == 'velocyto':
        if assumption not in ['constant_velocity', 'constant_unspliced']:
            raise RuntimeError(f'velocyto assumption {assumption} is not available.')

    # S
    S = pd.read_csv(session.data_file('counts_s.csv'), index_col=0)
    # U
    U = pd.read_csv(session.data_file('counts_u.csv'), index_col=0)
    # pop
    pop = pd.read_csv(session.data_file('meta.csv')).loc[:, 'pop'].to_numpy()

    if log_scale:
        S = np.log(S + 1)
        U = np.log(U + 1)

    if method == 'scvelo':
        adata = ad.AnnData(S)
        adata.layers['spliced'] = S
        # U
        adata.layers['unspliced'] = U
        # pop
        adata.obs['clusters'] = pd.Categorical(pop, categories=np.unique(pop), ordered=False)

        # scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
        scv.pp.normalize_per_cell(adata)
        scv.pp.log1p(adata)
        scv.pp.moments(adata, n_pcs=n_pcs, n_neighbors=n_nbs)

        # visualization
        scv.tl.umap(adata)
        scv.tl.louvain(adata)
        scanpy.tl.tsne(adata)

        if scv_mode == 'dynamical':
            # scv.tl.velocity(adata, mode="deterministic")
            # scv.tl.velocity_graph(adata)
            scv.tl.recover_dynamics(adata, var_names=var_names)
        scv.tl.velocity(adata, mode=scv_mode)
        scv.tl.velocity_graph(adata)

        velocity_vector = np.nan_to_num(adata.layers['velocity'], nan=0)

        if (layer == 'imputed'):
            imputed = adata.layers['Ms']
            expression_future = imputed + velocity_vector
        else:
            expression_future = adata.layers['spliced'] + velocity_vector

        transition_matrix = scv.tl.transition_matrix(adata)

        return VeloResult(
            adata, velocity_vector, expression_future, transition_matrix
        )
    elif method == 'velocyto':
        loom_fn = session.res_file(f'velocyto_{assumption}.loom')
        S_t = S.to_numpy().T
        U_t = U.to_numpy().T

        layers = {
            "": S_t,
            "spliced": S_t,
            "unspliced": U_t,
            "ambiguous": coo_matrix(S_t.shape, dtype=np.int)
        }
        row_attrs = {
            "Gene": S.columns.to_numpy()
        }
        col_attrs = {
            "CellID": S.index.values,
            "Cluster": pop
        }
        loompy.create(loom_fn, layers=layers, row_attrs=row_attrs, col_attrs=col_attrs)
        vlm = velocyto.VelocytoLoom(loom_fn)
        vlm._normalize_S(relative_size=vlm.initial_cell_size,
                target_size=np.mean(vlm.initial_cell_size))
        vlm._normalize_U(relative_size=vlm.initial_Ucell_size,
                target_size=np.mean(vlm.initial_Ucell_size))
        vlm.perform_PCA()
        vlm.knn_imputation(k = n_nbs, n_pca_dims=n_pcs)
        vlm.fit_gammas()
        vlm.predict_U()
        vlm.calculate_velocity()
        vlm.calculate_shift(assumption=assumption)
        vlm.extrapolate_cell_at_t(delta_t=1.)

        velocity_vector = vlm.Sx_sz_t - vlm.Sx

        return VeloResult(
            vlm, velocity_vector.T, None, None
        )


def scv_inspect(adata):
    scv.pl.proportions(adata)
    scanpy.pl.umap(adata)
    scanpy.pl.tsne(adata)
    scv.pl.velocity_embedding_stream(adata, basis='tsne')


def eval_velo(session: BenchSession, res: VeloResult):
    velo = res.velocity
    gt = pd.read_csv(session.data_file('velo.csv'), index_col=0).to_numpy()

    _genes_subset = ~np.isnan(velo).any(axis=0)
    velo = velo[:, _genes_subset]
    gt = gt[:, _genes_subset]

    return cosine_sim(velo, gt)


def eval_all(datasets: List[str]):
    configs = [
        (1, "velocyto", { "assumption": "constant_velocity" }),
        (2, "velocyto", { "assumption": "constant_unspliced" }),
        (3, "scvelo", { "scv_mode": "stochastic" }),
        (4, "scvelo", { "scv_mode": "deterministic" }),
        (5, "scvelo", { "scv_mode": "dynamical", "var_names": "all" }),
        # (6, "scvelo", { "scv_mode": "dynamical_residuals", "var_names": "all" }),
    ]

    config_names = [
        "velocyto_constant_velocity",
        "velocyto_constant_unspliced",
        "scvelo_stochastic",
        "scvelo_determinstic",
        "scvelo_dynamical",
        # "scvelo_dynamical_residuals"
    ]

    scores = { k: [] for k in config_names}
    sessions = map(lambda x: BenchSession(x), datasets)

    for session in sessions:
        for conf_id, method, args in configs:
            res = bench_velo(method, session, log_scale=True, **args)
            conf_name = config_names[conf_id - 1]
            scores[conf_name].append(eval_velo(session, res))

    return scores