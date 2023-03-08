# setwd("/home/[USER]/scMultiSim")
# devtools::load_all(".")

library("scMultiSim")

# Set your root directory to save data
ROOT_DIR <- "./scMultiSim/unif/"
# ROOT_DIR <- "/home/[USER]/scMultiSim/bench/unif/1"

.unif_fd <- function(conf) {
  fn <- sprintf("tree%g_%g_cells%g_genes_sigma%g_%g",
                conf$tree, conf$ncells, conf$ngenes, conf$sigma, conf$seed)
  fd <- file.path(ROOT_DIR, fn)
  dir.create(fd, showWarnings = F, recursive = T)
  fd
}

.unif_fn <- function(conf) {
  file.path(.unif_fd(conf), "res.rds")
}

configs1 <- expand.grid(
  tree = c(1),
  ngenes = c(200, 500),
  ncells = c(500),
  sigma = c(0.1, 0.5),
  seed = 1:4
) %>% split(., seq(nrow(.)))

configs3 <- expand.grid(
  tree = c(3),
  ngenes = c(500),
  ncells = c(800),
  sigma = c(0.1, 0.5),
  seed = 4
) %>% split(., seq(nrow(.)))

configs5 <- expand.grid(
  tree = c(5),
  ngenes = c(110, 200, 500),
  ncells = c(500),
  sigma = c(0.1, 0.5),
  seed = 1:4
) %>% split(., seq(nrow(.)))

# Set config set; configs1=ML, configs3=MT, configs5=MD
CURR_CONFIG <- configs1

# ==============================================================================

ctp3 <- (function(){
  set.seed(1)
  ctp <- cci_cell_type_params(Phyla3(), 6, 3:6, 1, rand = F)
  a <- c(1, 2, 1, 4, 1)
  b <- c(2, 3, 3, 5, 4)
  for (m in seq_along(a)) {
    i <- a[m]; j <- b[m]
    n_pairs <- sample(4:6, 1)
    pairs <- sample(1:6, n_pairs)
    for (p in pairs)
      ctp$params[i, j, p] <- ctp$params[j, i, p] <- 1
  }
  ctp
})()

ctp1 <- (function(){
  set.seed(1)
  ctp <- cci_cell_type_params(Phyla1(), 6, 3:6, 0.2, rand = F)
  for (i in 3:5) {
    for (j in 3:5) {
      if (i == j) next
      n_pairs <- sample(5:6, 1)
      pairs <- sample(1:6, n_pairs)
      for (p in pairs)
        ctp$params[i, j, p] <- ctp$params[j, i, p] <- 1
    }
  }
  ctp
})()


.sim_unif_dataset <- function(conf) {
  # reg: 2, 6, 10, 19, 80, 91
  grn_params <- GRN_params_100
  
  lig_params <- data.frame(
    target    = c(2,   6,   10,  19,  80,  91),
    regulator = c(101, 102, 103, 104, 105, 99),
    effect    = 5
  )
  
  tree <- if (conf$tree == 5) Phyla5() else if (conf$tree == 1) Phyla1() else Phyla3()
  cci <- if(conf$tree == 5) {
    list(
      params = lig_params,
      max.neighbors = 4,
      cell.type.lr.pairs = 3:6,
      step.size = 1,
      grid.size = ceiling(sqrt(conf$ncells) * 2),
      cell.type.interaction = "random",
      same.type.prob = 0.8,
      layout = "enhanced"
    )
  } else if(conf$tree == 1) {
    list(
      params = lig_params,
      max.neighbors = 4,
      cell.type.interaction = ctp1,
      cell.type.lr.pairs = 3:6,
      step.size = 0.2,
      grid.size = ceiling(sqrt(conf$ncells) * 2),
      same.type.prob = 0.8,
      layout = "enhanced"
    )
  } else {
    list(
      params = lig_params,
      max.neighbors = 4,
      cell.type.interaction = ctp3,
      cell.type.lr.pairs = 3:6,
      step.size = 1,
      grid.size = ceiling(sqrt(conf$ncells) * 1.5),
      same.type.prob = 0.8,
      layout = "enhanced2"
    )
  }
  
  giv_mean <- if (conf$ngene > 110) -1 else 0
  
  options_ <- list(
    rand.seed = conf$seed,
    threads = 1,
    GRN = grn_params,
    num.genes = conf$ngene,
    num.cells = conf$ncell,
    num.cifs = 50,
    tree = tree,
    cif.sigma = conf$sigma,
    intrinsic.noise = 1,
    cif.sigma = conf$sigma,
    giv.mean = giv_mean,
    giv.sd = 0.2,
    diff.cif.fraction = 0.8,
    cci = cci,
    debug = F
  )
  
  if (conf$tree == 5) {
    options$discrete.cif = T
    options$discrete.min.pop.size = ceiling(conf$ncells / 6)
  }
  
  cat("start",
      file= file.path(.unif_fd(conf), "sig.txt"), append = T)
  
  results <- sim_true_counts(options_)
  cat(format(Sys.time(), "%a %b %d %X %Y"),
      file= file.path(.unif_fd(conf), "sig.txt"), append = T)
  add_expr_noise(results)
  divide_batches(results)
  results
}

.write_unif_dataset <- function(conf, res) {
  fd <- .unif_fd(conf)
  saveRDS(res, file.path(fd, "res.rds"))
}

cat("===========BEGIN==========\n")

(function() {
  args <- commandArgs(trailingOnly=TRUE)
  conf <- configs5[[as.integer(args[1])]]
  print(.unif_fd(conf))
  res <- .sim_unif_dataset(conf)
  .write_unif_dataset(conf, res)
})()

