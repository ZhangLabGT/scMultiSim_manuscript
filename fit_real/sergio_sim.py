import os
import numpy as np
import pandas as pd
from SERGIO.sergio import sergio


def run_sim(ds, ngenes, ncells):
    sim = sergio(
        number_genes=ngenes,
        number_bins=4,
        number_sc=ncells,
        noise_params=1,
        decays=0.8,
        sampling_state=15,
        noise_type="dpd",
    )
    sim.build_graph(
        input_file_taregts=f"grn_sergio_{ds}.txt",
        input_file_regs=f"tf_sergio_{ds}.txt",
        shared_coop_state=2,
    )
    sim.simulate()
    expr = sim.getExpressions()

    def filter_nan_rows_cols(arr):
        # Check if any element in a row/column is not NaN
        mask_rows = ~np.isnan(arr).all(axis=(1, 2))
        mask_cols = ~np.isnan(arr).all(axis=(0, 2))
        mask_depth = ~np.isnan(arr).all(axis=(0, 1))

        def filter(a):
            return a[mask_rows][:, mask_cols][:, :, mask_depth]
        return filter

    expr = np.nan_to_num(expr)
    filter = filter_nan_rows_cols(expr)
    expr_clean = np.concatenate(expr, axis=1)

    print("clean")

    """
    Add outlier genes
    """
    expr_O = sim.outlier_effect(expr, outlier_prob=0.01, mean=2, scale=2)

    """
    Add Library Size Effect
    """
    libFactor, expr_O_L = sim.lib_size_effect(expr_O, mean=5.2, scale=0.2)

    """
    Add Dropouts
    """
    binary_ind = sim.dropout_indicator(expr_O_L, shape=3.5, percentile=62)
    expr_O_L_D = np.multiply(binary_ind, expr_O_L)

    """
    Convert to UMI count
    """
    count_matrix = sim.convert_to_UMIcounts(expr_O_L_D)

    """
    Make a 2d gene expression matrix
    """
    count_matrix = np.concatenate(count_matrix, axis=1)

    # write csv
    np.savetxt(f'expr_{ds}.csv', count_matrix, delimiter=',')
    np.savetxt(f'expr_clean_{ds}.csv', expr_clean, delimiter=',')

if __name__ == "__main__":
    run_sim("10x", 145, 735)
    run_sim("issac", 140, 750)
    run_sim("merfish", 92, round(2853 / 4))
    run_sim("seqfish", 75, round(523 / 4))
