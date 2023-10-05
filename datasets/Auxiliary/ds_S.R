# devtools::load_all(".")
library("scMultiSim")


(function() {
  .output_cci_res <- function(res, root) {
    message(file.path(root, "res.rds"))
    dir.create(root, recursive = T)
    write.csv(res$counts, file = file.path(root, "counts.csv"), quote = F)
    write.csv(res$cell_meta, file = file.path(root, "cell_info.csv"), quote = F)
    write.csv(res$grn_params, file = file.path(root, "grn.csv"), row.names = F, quote = F)
    write.csv(res$cci_cell_type_param, file = file.path(root, "cci.csv"), row.names = F, quote = F)
    write.csv(res$cci_locs, file = file.path(root, "locs.csv"), row.names = T, quote = F)
    saveRDS(res, file.path(root, "res.rds"))
  }
  lig_params <- data.frame(
    target    = c(2,   6,   10,   8,  20,  30),
    regulator = c(101, 102, 103, 104, 105, 106),
    effect    = 5
  )
  
  ctp <- cci_cell_type_params(Phyla1(), 6, 3:6, 0.3, rand = F)
  for (i in 2:4) {
    for (j in 2:4) {
      if (i == j) next
      n_pairs <- sample(4:6, 1)
      pairs <- sample(1:6, n_pairs)
      for (p in pairs)
        ctp$params[i, j, p] <- ctp$params[j, i, p] <- 1
    }
  }
  
  # foreach (seed = 1:8) %dopar% {
  for (seed in 1:8) {
    options_ <- list(
      rand.seed = seed,
      GRN = GRN_params_100,
      num.genes = 120,
      num.cells = 400,
      num.cifs = 50,
      tree = Phyla1(),
      intrinsic.noise = 1,
      cci = list(
        params = lig_params,
        max.neighbors = 4,
        cell.type.interaction = ctp,
        cell.type.lr.pairs = 3:6,
        step.size = 0.3
      ),
      debug = F
    )
    
    results <- sim_true_counts(options_)
    .output_cci_res(results, sprintf("/home/[USER]/scMultiSim/bench/cci/cci_sc_%g", seed))
  }
})()