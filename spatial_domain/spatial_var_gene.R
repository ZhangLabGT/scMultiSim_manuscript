devtools::load_all(".")

ROOT <- "bench/figures/spatial_domain"
fp <- \(x) file.path(ROOT, x)

dir.create(fp("scHybridNMF"), showWarnings = FALSE)
dir.create(fp("stagate"), showWarnings = FALSE)

set.seed(0)

.simSpatial <- function(seed) {
  set.seed(seed)

  lig_params <- data.frame(
    target    = c(101, 102),
    regulator = c(103, 104),
    effect    = c(5.2, 5.9)
  )

  zone_gt <- NULL
  new_center <- NULL
  locs <- NULL

  layout0 <- function(grid_size, final_types) {
    message("layout0")
    ncells <- length(final_types)
    grid_center <- c(round(grid_size / 2), round(grid_size / 2))
    all_locs <- gen_clutter(ncells, grid_size, grid_center)
    # center is bottom-left
    left_ones <- which(all_locs[,1] == min(all_locs[,1]))
    new_center <<- all_locs[left_ones[which.min(all_locs[left_ones, 2])],]
    dist_to_center <- sqrt(colSums((t(all_locs) - new_center)^2))
    new_locs <- all_locs[order(dist_to_center),]
    # prob of a cell type being in a zone (cell_type x zone)
    ct_matrix <- matrix(c(
      0.9, 0.1, 0.0,
      0.1, 0.8, 0.1,
      0.1, 0.7, 0.2,
      0.0, 0.1, 0.9
    ), nrow = 4, byrow = TRUE)
    # number of cells per type
    ct_pop <- c(160, 80, 100, 140)
    pop_mtx <- round(ct_matrix * ct_pop)
    if (sum(pop_mtx) != ncells) {
      diffrence <- ncells - sum(pop_mtx)
      pop_mtx[1, 1] <- pop_mtx[1, 1] + diffrence
    }
    # number of cells per zone
    zone_pop <- colSums(pop_mtx)
    # assign cells to zones
    cs <- cumsum(zone_pop)
    # sample cells
    cell_idx <- lapply(1:3, function(izone) {
      sample(rep(1:4, pop_mtx[,izone]), zone_pop[izone])
    }) %>% unlist()
    locs <<- new_locs[order(cell_idx),]
    zone_gt <<- rep(1:3, zone_pop)[order(cell_idx)]
    return(locs)
  }

  mod_cif <- function(i, cell_i, path_i, cif, giv) {
    list(cif, giv * 0.1)
  }

  ext_cif <- function(i) {
    message("ext_cif")
    spatial_genes <- 290:300
    dist_to_center <- colSums((t(locs) - new_center)^2)
    dist_to_center <- dist_to_center / max(dist_to_center)
    if (i == 3) {
      idx <- list(1:2, 3:4)
      ex_cif <- cbind(
        rnorm(500, dist_to_center, 0.02),
        rnorm(500, dist_to_center, 0.05),
        rnorm(500, 0.5 * dnorm(abs(dist_to_center - 0.5), 0, 0.04), 0.02),
        rnorm(500, 0.5 * dnorm(abs(dist_to_center - 0.5), 0, 0.04), 0.02)
      )
      ex_giv <- matrix(0, nrow = 300, ncol = 4)
      for (i in spatial_genes) {
        ex_giv[i, idx[[i %% 2 + 1]]] <- rnorm(2, 1, 0.5)
      }
      list(ex_cif, ex_giv * 2)
    } else {
      NULL
    }
  }

  res0 <- sim_true_counts(list2(
    num.cells = 500,
    num.genes = 300,
    num.cifs = 40,
    GRN = NA,
    speed.up = T,
    cif.sigma = 0.8,
    tree = ape::read.tree(text = "(A:1,B:1,C:1,D:1);"),
    diff.cif.fraction = 0.8,
    discrete.cif = T,
    discrete.pop.size = as.integer(c(120,150,100,130)),
    do.velocity = F,
    mod.cif.giv = mod_cif,
    ext.cif.giv = ext_cif,
    debug = T,
    cci = list(
      params = lig_params,
      max.neighbors = 4,
      cell.type.interaction = "random",
      layout = layout0,
      step.size = 1
    )
  ))

  res0$zone_gt <- as.character(zone_gt)
  res0
}


res <- .simSpatial(0)
stop()

saveRDS(res, fp("sim.rds"))

p <- plot_cell_loc(res, show.arrows = FALSE)
ggsave(fp("cell_type.pdf"), p, width = 8, height = 8)
p

p <- plot_cell_loc(res, show.arrows = FALSE, .cell.pop = log(res$counts[300,] + 1)) + scale_colour_viridis_c(direction = -1)
ggsave(fp("gene_300.pdf"), p, width = 8, height = 8)
p
p <- plot_cell_loc(res, show.arrows = FALSE, .cell.pop = log(res$counts[299,] + 1)) +
  scale_colour_viridis_c(direction = -1, begin = 0.3, end = 1)
ggsave(fp("gene_299.pdf"), p, width = 8, height = 8)
p
p <- plot_cell_loc(res, show.arrows = FALSE, .cell.pop = log(res$counts[298,] + 1)) + scale_colour_viridis_c(direction = -1)
ggsave(fp("gene_298.pdf"), p, width = 8, height = 8)
p
p <- plot_cell_loc(res, show.arrows = FALSE, .cell.pop = log(res$counts[297,] + 1)) +
scale_colour_viridis_c(direction = -1, begin = 0.2, end = 1) + theme_bw()
ggsave(fp("gene_297.pdf"), p, width = 7, height = 7)
p

p <- plot_cell_loc(res, show.arrows = FALSE, .cell.pop = log(res$counts[293,] + 1)) +
  scale_colour_viridis_c(direction = -1, begin = 0.5, end = 1) + theme_bw()
ggsave(fp("gene_293.pdf"), p, width = 8, height = 8)

p <- plot_cell_loc(res, .cell.pop = res$zone_gt, show.arrows = FALSE)
ggsave(fp("domain_gt.pdf"), p, width = 8, height = 8)
p

write.table(t(log2(res$counts + 1)), file = fp("scHybridNMF/expression_1.txt"),
            row.names = F, col.names = F, quote = F, sep = " ")

write.table(rownames(res$counts), file = fp("scHybridNMF/genes_1.txt"),
            row.names = F, col.names = F, quote = F, sep = " ")

# celltypes_
write.table(res$cell_meta$pop, file = fp("scHybridNMF/celltypes_1.txt"),
            row.names = F, col.names = F, quote = F, sep = " ")

write.csv(
  data.frame(
    "cell type" = res$cell_meta$pop,
    "coor X" = res$cci_locs[,1],
    "coor Y" = res$cci_locs[,2]
  ),
  file = fp("scHybridNMF/meta_1.txt"),
  row.names = F, col.names = T, quote = F, sep = ",")


.getEdges <- function(locs) {
  str <- ""
  ncells <- nrow(locs)
  for (i in seq_len(ncells)) {
    for (j in seq_len(ncells)) {
      if (i == j) {
        next
      }
      d <- sqrt(sum((locs[i,] - locs[j,])^2))
      if (d < 1.1) {
        str <- paste0(str, i - 1, " ", j - 1, "\n")
      }
    }
  }
  str
}

write(.getEdges(res$cci_locs), file = fp("scHybridNMF/neighborhood_1.txt"))

write.csv(log2(res$counts + 1), file = fp("stagate/expression.txt"),
          row.names = T, quote = F, sep = "\t")
write.csv(
  data.frame(
    "xcoord" = res$cci_locs[,1],
    "ycoord" = res$cci_locs[,2]
  ),
  file = fp("stagate/locations.csv"),
  row.names = T, col.names = T, quote = F, sep = ",")

write.csv(res$zone_gt, file = fp("stagate/domains.csv"),
          row.names = T, col.names = T, quote = F, sep = ",")

stop()

meta <- read.csv(fp("scHybridNMF/meta_1.txt"), header = T)
p <- plot_cell_loc(NULL, .cell.pop = as.character(meta$cell.type), show.arrows = FALSE, .locs = t(meta[,2:3]))
ggsave(fp("cell_type_gt.pdf"), p, width = 8, height = 8)

res_schybridnmf <- read.csv(fp("scHybridNMF/C.csv"), header = F)
p <- plot_cell_loc(NULL, .cell.pop = factor(as.character(res_schybridnmf[,1]), levels = c("3", "1", "2")),
                   show.arrows = FALSE, .locs = t(meta[,2:3]))
ggsave(fp("cell_type_schybridnmf.pdf"), p, width = 8, height = 8)

res_stagate <- read.csv(fp("stagate/stagate_out.csv"), header = T)
p <- plot_cell_loc(NULL, .cell.pop = factor(as.character(res_stagate$louvain + 1), levels = c("3", "1", "2")),
                   show.arrows = FALSE, .locs = t(meta[,2:3]))
ggsave(fp("cell_type_stagate.pdf"), p, width = 8, height = 8)

generateSpatialLoc()
