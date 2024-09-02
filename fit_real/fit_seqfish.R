library("devtools")
devtools::load_all(".")
library("scMultiSim")

library(Giotto)
library(SingleCellExperiment)

SS_seqfish <- readRDS("figures/real_data/SS_seqfish.rds")
spatial_all_scores <- readRDS("figures/real_data/spatial_all_scores.rds")
sc_grn <- read.delim("~/vm_shared/out_seqfish.txt", header = T)[1:200,]

selected_spat = spatial_all_scores[p.adj <= 0.01 & abs(log2fc) > 0.25 & lig_nr >= 4 & rec_nr >= 4]
data.table::setorder(selected_spat, -PI)
# ====
top_LR_ints = unique(selected_spat[order(-abs(PI))]$LR_comb)[1:33]
top_LR_cell_ints = unique(selected_spat[order(-abs(PI))]$LR_cell_comb)[1:33]

lrpair_used <- selected_spat[lig_nr > 10 & rec_nr > 10 & lig_expr > 0 & rec_expr > 0]
lrpair_used <- lrpair_used[log2fc > 0.7 | log2fc < -1]

gene_metadata = fDataDT(SS_seqfish)
mtx <- as.matrix(SS_seqfish@expression$cell$rna$custom@exprMat)
hvf <- gene_metadata$hvf == 'yes' & gene_metadata$perc_cells > 4 & gene_metadata$mean_expr_det > 0.5

sc_grn_genes <- unique(c(sc_grn$TF, sc_grn$target))

.remove <- numeric()
for (i in 1:nrow(sc_grn)) {
  if ((sc_grn[i, 1] %in% sc_grn[1:i, 2]) && (sc_grn[i, 2] %in% sc_grn[1:i, 1])) {
    .remove <- c(.remove, i)
  }
}
sc_grn <- sc_grn[-.remove, ]

in_grn <- rownames(mtx) %in% sc_grn_genes
table(in_grn)

other_g <- sample(which(!in_grn), 200 - sum(in_grn))
sce <- as.matrix(SS_seqfish@expression$cell$rna$raw@exprMat)
sce <- sce[sort(c(which(in_grn), other_g)),]
sce <- SingleCellExperiment(assays = list(counts = as.matrix(sce)),
                            colData = data.frame(ident = SS_seqfish@cell_metadata$cell$rna$cell_types))

sc_grn$importance <- sc_grn$importance * 100


lrpair_used <- selected_spat[lig_nr > 10 & rec_nr > 10 & lig_expr > 0 & rec_expr > 0]
lrpair_used <- lrpair_used[log2fc > 0.7 | log2fc < -1]


#########

grn_out0 <- sc_grn
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
grn_out <- grn_out[order(grn_out$importance, decreasing = TRUE),]

#########

uniq_genes <- unique(c(grn_out[,"TF"], grn_out[,"target"]))
remaining_genes <- setdiff(rownames(mtx), uniq_genes)

sim_dyngen <- function () {
  library(tidyverse)
  library(dyngen)
  library(tibble)
  library(Matrix)
  set.seed(1)
  sc_grn_1 <- sc_grn

  .dyn_module_info <-
    rbind(
      data.frame(
        module_id = c("Burn1", "Burn2"),
        basal = 1,
        burn = T, #sc_grn_genes %in% sc_grn$TF,
        independence = 1
      ),
      data.frame(
        module_id = unique(c(sc_grn_1$TF, sc_grn_1$target)),
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
        to = unique(setdiff(sc_grn_1$TF, sc_grn_1$target)),
        effect = as.integer(1),
        strength = 1,
        hill = 2
      ),
      data.frame(
        from = sc_grn_1$TF,
        to = sc_grn_1$target,
        effect = as.integer(1),
        strength = sc_grn_1$importance,
        hill = 2
      )
    )

  backbone <- backbone(
    module_info = as_tibble(.dyn_module_info),
    module_network = as_tibble(.dyn_network),
    expression_patterns = tribble(
      ~from, ~to,  ~module_progression, ~start, ~burn, ~time,
      "s0",  "s1", "+Burn1,+Burn2",               TRUE,   TRUE,  30,
      "s1",  "s2", "+Slc7a5",               FALSE,   FALSE,  300,
    )
  )

  .dyn_config <-
    initialise_model(
      backbone = backbone,
      num_tfs = nrow(backbone$module_info),
      num_cells = 523,
      num_targets = 110,
      num_hks = 10,
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
        realcount = Matrix(data = sce, sparse = T)
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
  tree10 <- ape::read.tree(text = '(A:1.0,B:1.0,C:1.0,D:1.0:,E:1.0,F:1.0,G:1.0,H:1.0,I:1.0,J:1.0);')

  sc_grn_scm <- sc_grn[,c(2, 1, 3)]

  cell_types_pop <- table(SS_seqfish@cell_metadata$cell$rna$cell_types)
  cell_types <- setNames(seq_along(cell_types_pop), names(cell_types_pop))
  n_lr_pair <- nrow(lrpair_used)
  ctp <- cci_cell_type_params(tree10, total.lr = n_lr_pair, step.size = 1, rand = F, discrete = T)
  for (k in seq(n_lr_pair)) {
    i <- cell_types[unlist(lrpair_used[k, 'lig_cell_type'])]
    j <- cell_types[unlist(lrpair_used[k, 'rec_cell_type'])]
    print(c(i, j, k))
    ctp$params[i, j, k] <- ctp$params[j, i, k] <- 1
  }

  lig_params <- data.frame(
    receptor = lrpair_used$receptor,
    ligand = lrpair_used$ligand,
    effect = lrpair_used$log2fc * 2
  )

  res <- sim_true_counts(list2(
    rand.seed = 1,
    threads = 1,
    GRN = sc_grn_scm,
    cif.sigma = 3,
    diff.cif.fraction = 0.8,
    num.genes = 200,
    num.cells = 523,
    num.cifs = 50,
    tree = tree10,
    discrete.cif = T,
    discrete.pop.size = cell_types_pop,
    do.velocity = F,
    intrinsic.noise = 1,
    giv.sd = 0.03,
    hge.prop = 0.05,
    hge.mean = 1.6,
    hge.sd = 0.5,
    hge.max.var = 80,
    scale.s = 0.8,
    cci = list(
      params = lig_params,
      cell.type.interaction = ctp,
      max.neighbors = 4,
      cell.type.interaction = "random",
      step.size = 1
    ),
    debug = T,
  ))

  set.seed(0)
  add_expr_noise(res, protocol = "UMI",
                 alpha_mean = 0.15, alpha_sd = 0.2,
                 alpha_gene_mean = 0.1, alpha_gene_sd = 0,
                 depth_mean = 40000)
  res
}

sim_scdesign3 <- function () {
  library(scDesign3)
  set.seed(0)
  scd3_res <- scdesign3(
    sce = sce,
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
  write(out_data, file = "figures/real_data/grn_sergio_seqfish.txt")

  tf_uniq <- setdiff(unique(tf), tg_uniq)
  tf_data <- c()
  for (i in tf_uniq) {
    tf_expr <- rnorm(4, 2, 1)
    tf_data <- c(tf_data, paste(id_map[i], tf_expr[1], tf_expr[2], tf_expr[3], tf_expr[4], sep = ","))
  }
  tf_data <- paste(tf_data, collapse = "\n")
  write(tf_data, file = "figures/real_data/tf_sergio_seqfish.txt")
  id_map
}

sergio_id_map <- sim_sergio1()
sergio_id_map <- setNames(names(sergio_id_map), sergio_id_map)

dg_out <- sim_dyngen()
res <- sim_scmultisim()
scd3_out <- sim_scdesign3()

stop()
sg_out_expr <- read.csv("~/gt/SERGIO/expr_seqfish.csv", header = FALSE)
rownames(sg_out_expr) <- sergio_id_map[as.character(as.integer(rownames(sg_out_expr)) - 1)]
sg_out_expr <- as.matrix(sg_out_expr)

p_dist <- compare_real_dist(
  counts(sce),
  res$counts_obs$counts,
  as.matrix(scd3_out$new_count),
  sg_out_expr,
  as.matrix(counts(dg_out$dataset)),
  labels = c("true", "scMultiSim", "scDesign3", "SERGIO", "dyngen")
)

ggsave("figures/real_data/dist_seqfish.pdf", p_dist, width = 20, height = 3.5)
