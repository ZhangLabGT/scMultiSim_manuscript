library("devtools")
devtools::load_all(".")
library("scMultiSim")


path <- "~/gt/scMultiSim/tmp/datasets/issaac-seq"
library(Seurat)

# Load the PBMC dataset
set.seed(0)
iss.data <- read.csv(file.path(path, "mCortex_all_gene_expression_matrix.csv.gz"), row.names = 1)

iss <- CreateSeuratObject(counts = iss.data, project = "issaac", min.cells = 3, min.features = 200)
iss <- FindVariableFeatures(iss, nfeatures = 2000)
iss <- FindVariableFeatures(iss, nfeatures = 1000)
iss <- NormalizeData(iss)
iss <- ScaleData(iss)
iss <- RunPCA(iss, npcs = 30)
iss <- FindNeighbors(object = iss)
iss <- FindClusters(object = iss, resolution = 0.5)
iss <- RunUMAP(iss, reduction = "pca", dims = 1:30)
DimPlot(iss, reduction = "umap")

sel_cells <- sample(colnames(iss), 3000)
# write.table(t(GetAssayData(iss[VariableFeatures(iss), sel_cells], slot = "scale.data")),
#             file = "~/vm_shared/iss_rna.csv",
#             quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

# tenx_data <- GetAssayData(iss[VariableFeatures(iss), sel_cells], slot = "count")
# stat_res <- get_stats(GetAssayData(iss[VariableFeatures(iss)], slot = "count"))

grn_out0 <- read.csv("~/vm_shared/out_iss.csv", sep = "\t", header = TRUE)[1:400,]
adj_list <- list()
for (i in seq(nrow(grn_out0))) {
  tf <- grn_out0[i, "TF"]
  tg <- grn_out0[i, "target"]
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

grn_out <- NULL
for (tf in names(adj_list)) {
  for (tg in adj_list[[tf]]) {
    df1 <- grn_out0[grn_out0$TF == tf & grn_out0$target == tg, "importance"]
    df2 <- grn_out0[grn_out0$TF == tg & grn_out0$target == tf, "importance"]
    print(c(tf, tg, max(df1, df2)))
    if (length(df1) > 0 || length(df2) > 0)
      grn_out <- rbind(grn_out, data.frame(TF = tf, target = tg, importance = max(df1, df2)))
  }
}
grn_out <- distinct(grn_out)
grn_out <- grn_out[order(grn_out$importance, decreasing = TRUE),][1:200,]

uniq_genes <- unique(c(grn_out$TF, grn_out$target))

remaining_genes <- setdiff(rownames(iss), uniq_genes)
tenx_data <- GetAssayData(iss[
                            c(uniq_genes, sample(remaining_genes, 1000 - length(uniq_genes))),
                            sel_cells
                          ], slot = "count")
tenx_sce <- as.SingleCellExperiment(iss[
                                      c(uniq_genes, sample(remaining_genes, 1000 - length(uniq_genes))),
                                      sel_cells
                                    ])

# grn_out
#
# idx <- Idents(iss)
# pop_sizes <- sapply(1:10, \(i) sum(idx == i))
idx_sel <- Idents(iss[,sel_cells])
pop_sizes <- sapply(0:10, \(i) sum(idx_sel == i))
c(1495, 1111, 305, 28)

# ===========================
# ATAC
iss.atac <- read.csv(file.path(path, "mCortex_all_ATAC_peak_matrix.csv.gz"), row.names = 1)
tenx_atac <- as.matrix(iss.atac[sample(1:nrow(iss.atac), 3000), colnames(tenx_data)])
atac.dens <- atac_dens_nonzero(tenx_atac)

sim_dyngen <- function () {
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
      "s1",  "s2", "+Rasgrf2",               FALSE,   FALSE,  300,
    )
  )

  .dyn_config <-
    initialise_model(
      backbone = backbone,
      num_tfs = nrow(backbone$module_info),
      num_cells = ncol(tenx_data),
      num_targets = length(uniq_genes),
      num_hks = 800, # length(remaining_genes),
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

  generate_dataset(
    .dyn_config,
    format = "sce",
    make_plots = FALSE
  )
}

sim_scmultisim <- function () {
  set.seed(0)
  options_ <- list(
    rand.seed = 0,
    GRN = data.frame(
      target = grn_out[,"target"],
      tf = grn_out[,"TF"],
      importance = grn_out[,"importance"] / 50
    ),
    num.genes = 1000,
    num.cells = 3000,
    num.cifs = 30,
    tree = ape::read.tree(text = '(A:1.0,B:1.0,C:1.0,D:1.0:,E:1.0,F:1.0,G:1.0,H:1.0,I:1.0,J:1.0,K:1.0);'),
    discrete.cif = T,
    discrete.pop.size = as.integer(pop_sizes),
    diff.cif.fraction = 0.6,
    do.velocity = F,
    atac.p_zero = sum(tenx_atac == 0) / length(tenx_atac),
    atac.density = atac.dens,
    riv.mean = 2,
    riv.prob = 0.25,
    riv.sd = 3,
    hge.prop = 0.15,
    hge.mean = 2,
    hge.sd = 100,
    hge.max.var = 1000000,
    cif.sigma = 10,
    giv.sd = 10,
    scale.s = 0.1
  )

  res <- sim_true_counts(options_)
  set.seed(0)
  add_expr_noise(res, alpha_mean = 0.01, alpha_sd = 0.5,
                 depth_mean = 50, depth_sd = 500,
                 nPCR1 = 16)

  # stat_res0 <- get_stats(res$counts_obs)
  # compare_res <- .compare_stats(stat_res0, stat_res)

  # plot_dist(res$counts_obs, tenx_data)
  # plot_dist_3(tenx_data, res$counts_obs, counts(dg_out$dataset))

  # plot_atac(tenx_atac, res$atacseq_data)
  res
}


sim_scdesign3 <- function () {
  library(scDesign3)
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

sim_sergio1 <- function () {
  # target gene id, number of target’s regulators, regulator ID_1,…, regulator ID_n, K_1,…,K_n, hill_coeff_1, …, hill_coeff_n
  out_data <- c()
  tf <- grn_out$TF
  tg <- grn_out$target
  all_genes <- unique(c(tf, tg))
  id_map <- setNames(seq_along(all_genes) - 1, all_genes)
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
  write(out_data, file = "figures/real_data/grn_sergio_issac.txt")

  tf_uniq <- setdiff(unique(tf), tg_uniq)
  tf_data <- c()
  for (i in tf_uniq) {
    tf_expr <- rnorm(4, 2, 1)
    tf_data <- c(tf_data, paste(id_map[i], tf_expr[1], tf_expr[2], tf_expr[3], tf_expr[4], sep = ","))
  }
  tf_data <- paste(tf_data, collapse = "\n")
  write(tf_data, file = "figures/real_data/tf_sergio_issac.txt")
  id_map
}

sergio_id_map <- sim_sergio1()
sergio_id_map <- setNames(names(sergio_id_map), sergio_id_map)

dg_out <- sim_dyngen()
res <- sim_scmultisim()
scd3_out <- sim_scdesign3()

sg_out_expr <- read.csv("~/gt/SERGIO/expr_issac.csv", header = FALSE)
rownames(sg_out_expr) <- sergio_id_map[as.character(as.integer(rownames(sg_out_expr)) - 1)]
sg_out_expr <- as.matrix(sg_out_expr)

p_dist <- compare_real_dist(
  as.matrix(tenx_data),
  res$counts_obs,
  as.matrix(scd3_out$new_count),
  sg_out_expr,
  as.matrix(counts(dg_out$dataset)),
  labels = c("true", "scMultiSim", "scDesign3", "SERGIO", "dyngen")
)

ggsave("figures/real_data/dist_issac.pdf", p_dist, width = 20, height = 3.5)
