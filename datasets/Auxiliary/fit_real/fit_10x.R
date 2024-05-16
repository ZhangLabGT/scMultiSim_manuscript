library("devtools")
devtools::load_all(".")

path <- "figures/real_data/pbmc_unsorted_3k"
library(Seurat)
library(scDesign3)

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

tenx_sce <- as.SingleCellExperiment(pbmc[VariableFeatures(pbmc)])
tenx_data <- GetAssayData(pbmc[VariableFeatures(pbmc)], slot = "count")
tenx_ct <- Idents(pbmc)

grn_out <- read.csv("~/vm_shared/out.csv", sep = "\t", header = TRUE)[1:400,]

adj_list <- list()
for (i in seq(nrow(grn_out))) {
  tf <- grn_out[i, "TF"]
  tg <- grn_out[i, "target"]
  if (is.null(adj_list[[tf]])) adj_list[[tf]] <- c()
  adj_list[[tf]] <- c(adj_list[[tf]], tg)
}
.remove_cycle <- function(node, visited) {
  if (is.null(adj_list[[node]])) return()
  if (any(visited %in% adj_list[[node]])) {
    adj_list[[node]] <<- setdiff(adj_list[[node]], visited)
    return()
  }
  for (i in adj_list[[node]]) {
    .remove_cycle(i, c(visited, i))
  }
}
for (i in names(adj_list)) {
  .remove_cycle(i, i)
}
grn_out0 <- grn_out
grn_out <- NULL
for (tf in names(adj_list)) {
  for (tg in adj_list[[tf]]) {
    df1 <- grn_out0[grn_out0$TF == tf & grn_out0$target == tg, "importance"]
    df2 <- grn_out0[grn_out0$TF == tg & grn_out0$target == tf, "importance"]
    if (length(df1) > 0 || length(df2) > 0)
      grn_out <- rbind(grn_out, data.frame(TF = tf, target = tg, importance = max(df1, df2)))
  }
}
grn_out <- grn_out[order(grn_out$importance, decreasing = TRUE),][1:200,]

uniq_genes <- unique(c(grn_out$TF, grn_out$target))

remaining_genes <- setdiff(rownames(pbmc), uniq_genes)
tenx_data <- GetAssayData(pbmc[
                            c(uniq_genes, sample(remaining_genes, 1000 - length(uniq_genes)))
                          ], slot = "count")

idx_sel <- Idents(pbmc)
pop_sizes <- sapply(0:3, \(i) sum(idx_sel == i))
# c(1495, 1111, 305, 28)

# ATAC
tenx_atac <- as.matrix(pbmc.atac[sample(1:nrow(pbmc.atac), 3000), colnames(tenx_data)])
atac.dens <- atac_dens_nonzero(tenx_atac)

# ############################
# Simulating
# ############################


sim_dyngen <- function() {
  library(tidyverse)
  library(dyngen)
  library(tibble)
  library(Matrix)
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
}

sim_scmultisim <- function () {
  grn_out_scm <- data.frame(TF = character(), target = character(), importance = numeric())
  .find_parent <- function(tg) {
    tf <- grn_out[grn_out$target == tg, "TF"]
    if (length(tf) == 0) {
      return(tg)
    } else if (length(tf) == 1) {
      return(.find_parent(tf))
    } else {
      unique(unlist(sapply(tf, .find_parent, USE.NAMES = F)))
    }
  }

  for (i in 1:nrow(grn_out)) {
    tf <- grn_out[i, "TF"]
    target <- grn_out[i, "target"]
    importance <- grn_out[i, "importance"]
    p <- .find_parent(target)
    for (j in p) {
      grn_out_scm <- rbind(grn_out_scm,
                           data.frame(TF = j, target = target, importance = importance))
    }
  }
  grn_out_scm <- unique(grn_out_scm)
  remain_rows <- c()
  for (i in 1:nrow(grn_out_scm)) {
    tf <- grn_out_scm[i, "TF"]
    target <- grn_out_scm[i, "target"]
    rows <- grn_out_scm$target == target & grn_out_scm$TF == tf
    if (sum(rows) > 1) {
      remain_rows <- c(remain_rows, which(rows)[1])
    } else {
      remain_rows <- c(remain_rows, i)
    }
  }
  grn_out_scm <- grn_out_scm[unique(remain_rows),]

  set.seed(0)
  options_ <- list(
    rand.seed = 1,
    GRN = data.frame(
      target = grn_out_scm[,"target"],
      tf = grn_out_scm[,"TF"],
      importance = grn_out_scm[,"importance"] / 5
    ),
    # grn.effect = 5,
    num.genes = 1000,
    num.cells = ncol(pbmc),
    num.cifs = 50,
    tree = rtree(4),
    discrete.cif = T,
    discrete.pop.size = as.integer(pop_sizes),
    diff.cif.fraction = 0.8,
    do.velocity = F,
    atac.effect = 0.2,
    atac.p_zero = sum(tenx_atac == 0) / length(tenx_atac),
    atac.density = atac.dens,
    riv.mean = 2,
    riv.prob = 0.15,
    riv.sd = 3,
    cif.mean = 0,
    cif.sigma = 0.1,
    giv.sd = 1,
    scale.s = 0.1
  )

  res <- sim_true_counts(options_)
  set.seed(0)
  add_expr_noise(res, alpha_mean = 0.001, alpha_sd = 0.05, depth_mean = 50, depth_sd = 300)
  set.seed(0)
  add_outliers(res, 0.01, 1000, 200, 4, max.var = 5000)
  res
}


sim_sergio1 <- function () {
  # target gene id, number of target’s regulators, regulator ID_1,…, regulator ID_n, K_1,…,K_n, hill_coeff_1, …, hill_coeff_n
  out_data <- c()
  tf <- grn_out$TF
  tg <- grn_out$target
  all_genes <- unique(c(tf, tg))
  id_map <- setNames(seq_along(all_genes), all_genes)
  tg_uniq <- unique(tg)
  for (i in tg_uniq) {
    tfs <- grn_out[grn_out$target == i, "TF"]
    out_str <- paste(id_map[i], length(tfs), sep = ",")
    for (j in id_map[tfs]) {
      out_str <- paste(out_str, j, sep = ",")
    }
    for (j in tfs) {
      out_str <- paste(out_str, grn_out[grn_out$target == i & grn_out$TF == j, "importance"], sep = ",")
    }
    for (j in tfs) {
      out_str <- paste(out_str, 2, sep = ",")
    }
    out_data <- c(out_data, out_str)
  }
  out_data <- paste(out_data, collapse = "\n")
  write(out_data, file = "figures/real_data/grn_sergio_10x.txt")

  tf_uniq <- setdiff(unique(tf), tg_uniq)
  tf_data <- c()
  for (i in tf_uniq) {
    tf_expr <- rnorm(4, 2, 1)
    tf_data <- c(tf_data, paste(id_map[i], tf_expr[1], tf_expr[2], tf_expr[3], tf_expr[4], sep = ","))
  }
  tf_data <- paste(tf_data, collapse = "\n")
  write(tf_data, file = "figures/real_data/tf_sergio_10x.txt")
  id_map
}

sergio_id_map <- sim_sergio1()
sergio_id_map <- setNames(names(sergio_id_map), sergio_id_map)

sim_scdesign3 <- function () {
  set.seed(0)
  scd3_res <- scdesign3(
    sce = tenx_sce,
    assay_use = "counts",
    celltype = "ident",
    pseudotime = NULL,
    spatial = NULL,
    other_covariates = NULL,
    mu_formula = "ident",
    sigma_formula = "1",
    family_use = "nb",
    n_cores = 1,
    usebam = FALSE,
    corr_formula = "1",
    copula = "gaussian",
    DT = TRUE,
    pseudo_obs = FALSE,
    return_model = FALSE,
    nonzerovar = FALSE
  )
}

# dg_out <- sim_dyngen()
res <- sim_scmultisim()
stop()
scd3_out <- sim_scdesign3()
sg_out <- read.csv("~/gt/SERGIO/expr_clean_10x.csv", header = FALSE)
rownames(sg_out) <- sergio_id_map[rownames(sg_out)]
sg_out <- as.matrix(sg_out)

# ############################
# Plotting
# ############################

ggmap_real <- plot_gg_heatmap(as.matrix(tenx_data)[uniq_genes,])

# plot_tsne(log2(res$counts + 1), res$cell_meta$pop, runPCA = T)

p_dist <- compare_real_dist(
  as.matrix(tenx_data),
  res$counts_obs,
  as.matrix(scd3_out$new_count),
  as.matrix(counts(dg_out$dataset)),
  labels = c("true", "scMultiSim", "scDesign3", "dyngen"))

ggsave("figures/real_data/dist_10x.pdf", p_dist, width = 20, height = 3.5)

plot_gg_heatmap(as.matrix(tenx_data)[uniq_genes,], row_order = ggmap_real$rowDendrogram,
                filename = "figures/real_data/gg_10x_true.png")

plot_gg_heatmap(res$counts[uniq_genes,], row_order = ggmap_real$rowDendrogram,
                filename = "figures/real_data/gg_10x_scmultisim.png")

plot_gg_heatmap(sg_out, row_order = ggmap_real$rowDendrogram,
                filename = "figures/real_data/gg_10x_sergio.png")

plot_gg_heatmap(
  as.matrix(counts(dg_out$dataset))[paste0(gsub( "[.-]", "_", uniq_genes), "_TF1"),],
  row_order = ggmap_real$rowDendrogram,
  filename = "figures/real_data/gg_10x_dyngen.png")

plot_gg_heatmap(as.matrix(scd3_out$new_count)[uniq_genes,], row_order = ggmap_real$rowDendrogram,
                filename = "figures/real_data/gg_10x_scdesign3.png")

plot_cc_heatmap(res$counts_obs, res$cell_meta$pop,
                filename = "figures/real_data/cc_10x_scmultisim.png")

# plot_cc_heatmap(sg_out)

plot_cc_heatmap(as.matrix(tenx_data), tenx_ct,
                filename = "figures/real_data/cc_10x_true.png")

# plot_cc_heatmap(as.matrix(counts(dg_out$dataset))[paste0(gsub( "[.-]", "_", uniq_genes), "_TF1"),],
#                 filename = "figures/real_data/cc_10x_dyngen.png")

plot_cc_heatmap(as.matrix(scd3_out$new_count), tenx_ct,
                filename = "figures/real_data/cc_10x_scdesign3.png")
