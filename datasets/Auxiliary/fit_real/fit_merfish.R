METHOD <- "scMultiSim"
# METHOD <- "dyngen"


library(MerfishData)

spe <- MouseHypothalamusMoffitt2018()
colnames(spe) <- paste0("cell_", 1:ncol(spe))

library(Seurat)
set.seed(0)
sel_cells <- sample(seq(ncol(spe)), 3000)
spe <- spe[,sel_cells]
counts <- matrix(0, nrow = nrow(spe), ncol = ncol(spe))
rownames(counts) <- rownames(spe)
for (g in rownames(spe)) {
  counts[g,] <- sapply(assay(spe, "molecules")[g,], \(x) nrow(x))
}

colnames(counts) <- colData(spe)$cell_id
mf <- CreateSeuratObject(counts = counts, project = "mf", min.cells = 3, min.features = 20)
mf <- FindVariableFeatures(mf, nfeatures = 2000)
mf <- NormalizeData(mf)
mf <- ScaleData(mf)
mf <- RunPCA(mf, npcs = 30)
mf <- FindNeighbors(object = mf)
mf <- FindClusters(object = mf, resolution = 0.5)
mf <- RunUMAP(mf, reduction = "pca", dims = 1:30)
DimPlot(mf, reduction = "umap")

# write.table(t(GetAssayData(mf, slot = "scale.data")),
#             file = "~/vm_shared/mf_rna.csv",
#             quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

grn_out <- read.csv("~/vm_shared/out_mf.csv", sep = "\t", header = TRUE)[1:200,]
uniq_genes <- unique(c(grn_out[,"TF"], grn_out[,"target"]))
remaining_genes <- setdiff(rownames(mf), uniq_genes)

idx_sel <- Idents(mf)
pop_sizes <- sapply(0:8, \(i) sum(idx_sel == i))
tenx_data <- GetAssayData(mf, slot = "count")

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
      "s1",  "s2", "+Tmem108",               FALSE,   FALSE,  300,
    )
  )

  .dyn_config <-
    initialise_model(
      backbone = backbone,
      num_tfs = nrow(backbone$module_info),
      num_cells = ncol(tenx_data),
      num_targets = length(uniq_genes),
      num_hks = length(remaining_genes),
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

} else {

  set.seed(0)
  options_ <- list(
    rand.seed = 0,
    GRN = data.frame(
      target = grn_out[,"target"],
      tf = grn_out[,"TF"],
      importance = grn_out[,"importance"] / 50
    ),
    num.genes = nrow(mf),
    num.cells = ncol(mf),
    num.cifs = 20,
    tree = rtree(9),
    discrete.cif = T,
    discrete.pop.size = as.integer(pop_sizes),
    diff.cif.fraction = 0.6,
    do.velocity = F,
    hge.prop = 0.1,
    hge.mean = 2,
    hge.sd = 5,
    hge.max.var = 100000,
    cif.sigma = 2,
    giv.sd = 2,
    scale.s = 1
  )

  res <- sim_true_counts(options_)
  set.seed(0)
  add_expr_noise(res, alpha_mean = 0.05, alpha_sd = 0.2,
                 depth_mean = 50, depth_sd = 200,
                 nPCR1 = 16)

  plot_dist_3(tenx_data, res$counts_obs, counts(dg_out$dataset))
}