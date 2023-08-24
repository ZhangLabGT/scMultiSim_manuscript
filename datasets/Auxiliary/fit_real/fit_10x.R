METHOD <- "scMultiSim"
# METHOD <- "dyngen"


library("devtools")
devtools::load_all(".")
library("scMultiSim")

path <- "~/gt/scMultiSim/tmp/datasets/pbmc_unsorted_3k"
library(Seurat)

# Load the PBMC dataset
set.seed(0)
pbmc.data <- Read10X_h5(file.path(path, "filtered_feature_bc_matrix.h5"))
pbmc.atac <- pbmc.data[["Peaks"]]
pbmc.rna <- pbmc.data[["Gene Expression"]]
rm(pbmc.data)
pbmc <- CreateSeuratObject(counts = pbmc.rna, project = "10x_multi_pbmc_3k", min.cells = 3, min.features = 200)
pbmc <- FindVariableFeatures(pbmc, nfeatures = 2000)
pbmc <- FindVariableFeatures(pbmc, nfeatures = 1000)
pbmc <- NormalizeData(pbmc)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc, npcs = 30)
pbmc <- FindNeighbors(object = pbmc)
pbmc <- FindClusters(object = pbmc, resolution = 0.02)
pbmc <- RunUMAP(pbmc, reduction = "pca", dims = 1:30)
DimPlot(pbmc, reduction = "umap")
rm(pbmc.rna)

# write.table(t(GetAssayData(pbmc[VariableFeatures(pbmc)], slot = "scale.data")),
#             file = "~/vm_shared/pbmc_rna.csv",
#             quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

tenx_data <- GetAssayData(pbmc[VariableFeatures(pbmc)], slot = "count")
# stat_res <- get_stats(GetAssayData(pbmc[VariableFeatures(pbmc)], slot = "count"))
#

grn_out <- read.csv("~/vm_shared/out.csv", sep = "\t", header = TRUE)[1:200,]
uniq_genes <- unique(c(grn_out[,"TF"], grn_out[,"target"]))
remaining_genes <- setdiff(rownames(pbmc), uniq_genes)
tenx_data <- GetAssayData(pbmc[
                            c(uniq_genes, sample(remaining_genes, 1000 - length(uniq_genes)))
                          ], slot = "count")

idx_sel <- Idents(pbmc)
pop_sizes <- sapply(0:3, \(i) sum(idx_sel == i))
# c(1495, 1111, 305, 28)

# ===========================
# ATAC
tenx_atac <- as.matrix(pbmc.atac[sample(1:nrow(pbmc.atac), 3000), colnames(tenx_data)])
atac.dens <- atac_dens_nonzero(tenx_atac)


# dyngen
if (METHOD == "dyngen") {

  library(tidyverse)
  library(dyngen)
  library(tibble)
  set.seed(1)

  grn_out$TF <- gsub("[.-]", "_", grn_out$TF)
  grn_out$target <- gsub("[.-]", "_", grn_out$target)

  .dyn_module_info <-
    rbind(
      data.frame(
        module_id = c("Burn1", "Burn2"),
        basal = 1,
        burn = T, #sc_grn_genes %in% sc_grn$TF,
        independence = 1
      ),
      data.frame(
        module_id = unique(c(grn_out$TF, grn_out$target)),
        basal = 1,
        burn = F, #sc_grn_genes %in% sc_grn$TF,
        independence = 1
      )
    )

  .dyn_network <-
    rbind(
      data.frame(
        from = "Burn1",
        to = "Burn2",
        effect = as.integer(1),
        strength = 1,
        hill = 2
      ),
      data.frame(
        from = "Burn2",
        to = unique(setdiff(grn_out$TF, grn_out$target)),
        effect = as.integer(1),
        strength = 1,
        hill = 2
      ),
      data.frame(
        from = grn_out$TF,
        to = grn_out$target,
        effect = as.integer(1),
        strength = grn_out$importance,
        hill = 2
      )
    )

  backbone <- backbone(
    module_info = as_tibble(.dyn_module_info),
    module_network = as_tibble(.dyn_network),
    expression_patterns = tribble(
      ~from, ~to,  ~module_progression, ~start, ~burn, ~time,
      "s0",  "s1", "+Burn1,+Burn2",               TRUE,   TRUE,  30,
      "s1",  "s2", "+PLXDC2",               FALSE,   FALSE,  300,
    )
  )

  .dyn_config <-
    initialise_model(
      backbone = backbone,
      num_tfs = nrow(backbone$module_info),
      num_cells = ncol(tenx_data),
      num_targets = length(uniq_genes),
      num_hks = 900, # length(remaining_genes),
      verbose = FALSE,
      tf_network_params = tf_network_default(
        min_tfs_per_module = 1L,
        sample_num_regulators = function() 1,
        weighted_sampling = FALSE
      ),
      feature_network_params = feature_network_default(
        realnet = NULL,
        damping = 0.01,
        target_resampling = Inf,
        max_in_degree = 5
      ),
      simulation_params = simulation_default(
        total_time = 1000,
        census_interval = 2,
        ssa_algorithm = ssa_etl(tau = 300/3600),
        experiment_params = simulation_type_wild_type(num_simulations = 10)
      ),
      experiment_params = experiment_snapshot(
        realcount = Matrix(data = tenx_data, sparse = T)
      )
    )

  .dyn_model <- generate_tf_network(.dyn_config)
  plot_feature_network(.dyn_model, show_targets = FALSE)

  dg_out <- generate_dataset(
    .dyn_config,
    format = "sce",
    make_plots = FALSE
  )


  plot_dist(counts(dg_out$dataset), tenx_data)
} else {

  set.seed(0)
  options_ <- list(
    rand.seed = 0,
    GRN = data.frame(
      target = grn_out[,"target"],
      tf = grn_out[,"TF"],
      importance = grn_out[,"importance"] / 50
    ),
    num.genes = 1000,
    num.cells = ncol(pbmc),
    num.cifs = 20,
    tree = rtree(4),
    discrete.cif = T,
    discrete.pop.size = as.integer(pop_sizes),
    diff.cif.fraction = 0.6,
    do.velocity = F,
    atac.p_zero = sum(tenx_atac == 0) / length(tenx_atac),
    atac.density = atac.dens,
    riv.mean = 2,
    riv.prob = 0.15,
    riv.sd = 3,
    hge.prop = 0.1,
    hge.mean = 2,
    hge.sd = 10,
    hge.max.var = 1000000,
    cif.sigma = 5,
    giv.sd = 5,
    scale.s = 0.1
  )

  res <- sim_true_counts(options_)
  set.seed(0)
  add_expr_noise(res, alpha_mean = 0.001, alpha_sd = 0.08,
                 depth_mean = 50, depth_sd = 200,
                 nPCR1 = 16)

  # stat_res0 <- get_stats(res$counts_obs)
  # compare_res <- .compare_stats(stat_res0, stat_res)

  # plot_dist_3(tenx_data, res$counts_obs, counts(dg_out$dataset))

  plot_atac(tenx_atac, res$atacseq_data)
}

