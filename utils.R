# setwd("/home/[USER]/scMultiSim")
# devtools::load_all(".")
library("scMultiSim")

# Dataset M helpers

configs1 <- expand.grid(
  tree = c(1),
  ngenes = c(110, 200, 500),
  ncells = c(500, 800),
  sigma = c(0.1, 0.5),
  seed = 1:4
) %>% split(., seq(nrow(.)))
configs3 <- expand.grid(
  tree = c(3),
  ngenes = c(110, 200, 500),
  ncells = c(500, 800),
  sigma = c(0.1, 0.5),
  seed = 1:4
) %>% split(., seq(nrow(.)))
configs5 <- expand.grid(
  tree = c(5),
  ngenes = c(110, 200, 500),
  ncells = c(500, 800),
  sigma = c(0.1, 0.5),
  seed = 1:4
) %>% split(., seq(nrow(.)))

configs_all <- c(configs1, configs3, configs5)


.conf_path <- function(conf, fn = NULL, grn = F) {
  ROOT_DIR <- "/home/[USER]/scMultiSim/bench/unif/0"
  GRN_DIR <- "/home/[USER]/scMultiSim/bench/unif"
  fd <- sprintf("tree%g_%g_cells%g_genes_sigma%g_%g",
                conf$tree, conf$ncells, conf$ngenes, conf$sigma, conf$seed)
  fd <- if (grn) {
    file.path(GRN_DIR, "grn_all", fd)
  } else {
    file.path(ROOT_DIR, fd)
  }
  if (is.null(fn)) {
    return(fd)
  } else {
    return(file.path(fd, fn))
  }
}

read_res <- function(conf) {
  readRDS(.conf_path(conf, "res.rds"))
}

# Dataset A helpers

af_list <- c(0.2, 0.5, 0.9)
configs_a <- expand.grid(
  velo = F,
  tree = 5,
  af = af_list,
  intr = c(0.3, 1),
  diff = c(0.8, 0.2),
  dis = c(T, F),
  sigma = c(0.1, 0.5, 1),
  seed = 1:2
) %>% split(., seq(nrow(.)))

.ds_A_name <- function(conf, fn = NULL) {
  dis_str <- if (conf$dis) "dis_" else "cnt_"
  out_dir <- sprintf("%stree_%d_velo_%g_af_%g_diff_%g_intr_%g_sig_%g_seed_%d",
                     dis_str, conf$tree, conf$velo, conf$af, conf$diff, conf$intr, conf$sigma, conf$seed)
  path <- file.path("/home/[USER]/scMultiSim/bench/a", out_dir)
  dir.create(path, showWarnings = F, recursive = T)
  if (is.null(fn)) {
    path
  } else {
    file.path(path, fn)    
  }
}