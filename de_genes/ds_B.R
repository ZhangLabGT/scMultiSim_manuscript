devtools::load_all(".")
library(Seurat)

# pop_size <- c(1053, 1220, 578, 171)
pop_size <- c(1550, 810, 237, 120)
num_cells <- sum(pop_size)
shift_mean <- c(3.4, 1.8, 3.1, 5.0)

.simCondEffect <- function(seed) {
  set.seed(seed)
  res0 <- sim_true_counts(list2(
    num.cells = num_cells,
    num.genes = 2000,
    # num.cifs = 40,
    num.cifs = 50,
    GRN = NA,
    speed.up = T,
    cif.sigma = 0.5,
    tree = ape::read.tree(text = "(A:1,B:1,C:1,D:1);"),
    diff.cif.fraction = 0.8,
    discrete.cif = T,
    discrete.pop.size = as.integer(pop_size),
    do.velocity = F,
  ))

  gc()
  modified_cif <- list()

  set.seed(seed)
  res_norm_2 <- sim_true_counts(list2(
    num.cells = num_cells,
    num.genes = 2000,
    num.cifs = 50,
    GRN = NA,
    speed.up = T,
    cif.sigma = 0.5,
    tree = ape::read.tree(text = "(A:1,B:1,C:1,D:1);"),
    diff.cif.fraction = 0.8,
    discrete.cif = T,
    discrete.pop.size = as.integer(pop_size),
    do.velocity = F,
    # debug = T,
    mod.cif.giv = function(i, cif, giv, meta) {
      message(sprintf("mod.cif %g", i))
      res <- matrix(cif, nrow = nrow(cif), ncol = ncol(cif))
      all_pop <- unique(meta$pop)
      cs <- colSums(giv > 0) - colSums(giv < 0)
      for (p in all_pop) {
        idx <- meta$pop == p
        # change_idx <- sample(which(cs > 0), 4, replace = FALSE)
        change_idx <- sample(11:50, 4, replace = FALSE)
        # ci1 <- change_idx[1:2]
        ci1 <- change_idx
        if (i == 1) {
          # res[idx, ci1] <- pmax(cif[idx, ci1] + runif(sum(idx) * ci1, 2, 4), 0)
          # res[idx, ci2] <- pmax(cif[idx, ci2] + runif(sum(idx) * ci2, 2, 6), 0)
        } else if (i == 2) {
          # res[idx, ci1] <- pmax(cif[idx, ci1] + runif(sum(idx) * ci1, -2, -1), 0)
          # res[idx, ci2] <- pmax(cif[idx, ci2] + runif(sum(idx) * ci2, -4, -2), 0)
        } else if (i == 3) {
          message(sprintf("-- pop %g", p))
          print(ci1)
          noise <- rnorm(sum(idx) * ci1, shift_mean[p], 0.1)
          print(noise)
          modified_cif[[p]] <<- list(ci1, noise)
          res[idx, ci1] <- cif[idx, ci1] + noise
          # res[idx, ci2] <- pmax(cif[idx, ci2] + runif(sum(idx) * ci2, -2, -1), 0)
        }
      }
      list(res, giv)
    },
  ))

  list(res0, res_norm_2, modified_cif)
}

.simCondEffect2 <- function(seed) {
  set.seed(seed)
  cif_shift <- runif(num_cells * 4, -1, 1) * 2.5

  set.seed(seed)
  res0 <- sim_true_counts(list2(
    num.cells = num_cells,
    num.genes = 2000,
    # num.cifs = 40,
    num.cifs = 50,
    GRN = NA,
    speed.up = T,
    cif.sigma = 0.5,
    tree = ape::read.tree(text = "(A:1,B:1,C:1,D:1);"),
    diff.cif.fraction = 0.8,
    discrete.cif = T,
    discrete.pop.size = as.integer(pop_size),
    do.velocity = F,
    mod.cif.giv = function(i, cif, giv, meta) {
      message(sprintf("mod.cif %g", i))
      res <- matrix(cif, nrow = nrow(cif), ncol = ncol(cif))
      res[,c(2,4,6,8)] <- res[,c(2,4,6,8)] + matrix(cif_shift, nrow = nrow(cif), ncol = 4)
      list(res, giv)
    },
  ))

  gc()
  modified_cif <- list()

  set.seed(seed)
  res_norm_2 <- sim_true_counts(list2(
    num.cells = num_cells,
    num.genes = 2000,
    num.cifs = 50,
    GRN = NA,
    speed.up = T,
    cif.sigma = 0.5,
    tree = ape::read.tree(text = "(A:1,B:1,C:1,D:1);"),
    diff.cif.fraction = 0.8,
    discrete.cif = T,
    discrete.pop.size = as.integer(pop_size),
    do.velocity = F,
    # debug = T,
    mod.cif.giv = function(i, cif, giv, meta) {
      message(sprintf("mod.cif %g", i))
      res <- matrix(cif, nrow = nrow(cif), ncol = ncol(cif))
      res[,c(2,4,6,8)] <- res[,c(2,4,6,8)] + matrix(cif_shift, nrow = nrow(cif), ncol = 4)
      all_pop <- unique(meta$pop)
      cs <- colSums(giv > 0) - colSums(giv < 0)
      for (p in all_pop) {
        idx <- meta$pop == p
        # change_idx <- sample(which(cs > 0), 4, replace = FALSE)
        change_idx <- sample(11:50, 4, replace = FALSE)
        # ci1 <- change_idx[1:2]
        # ci2 <- change_idx[3:4]
        ci1 <- change_idx
        if (i == 1) {
          # res[idx, ci1] <- pmax(cif[idx, ci1] + runif(sum(idx) * ci1, 2, 4), 0)
          # res[idx, ci2] <- pmax(cif[idx, ci2] + runif(sum(idx) * ci2, 2, 6), 0)
        } else if (i == 2) {
          # res[idx, ci1] <- pmax(cif[idx, ci1] + runif(sum(idx) * ci1, -2, -1), 0)
          # res[idx, ci2] <- pmax(cif[idx, ci2] + runif(sum(idx) * ci2, -4, -2), 0)
        } else if (i == 3) {
          message(sprintf("-- pop %g", p))
          print(ci1)
          noise <- rnorm(sum(idx) * ci1, shift_mean[p], 0.1)
          print(noise)
          modified_cif[[p]] <<- list(ci1, noise)
          res[idx, ci1] <- cif[idx, ci1] + noise
          # res[idx, ci2] <- pmax(cif[idx, ci2] + runif(sum(idx) * ci2, -2, -1), 0)
        }
      }
      list(res, giv)
    },
  ))

  list(res0, res_norm_2, modified_cif)
}


deDist <- function(x, sample_label, id1, id2, min.diff.pct = 0.01) {
  # if x is seurat object
  if (is(x, "Seurat")) {
    obj <- x[VariableFeatures(x),]
  } else {
    obj <- CreateSeuratObject(x, project = "de")
    print(obj)
    obj <- NormalizeData(obj)
    obj <- ScaleData(obj)
    obj <- FindVariableFeatures(obj, nfeatures = 2000)
    obj <- RunPCA(obj, verbose = FALSE)
  }
  Idents(obj) <- sample_label
  mono.de <- FindMarkers(obj, ident.1 = id1, ident.2 = id2, verbose = FALSE, min.diff.pct = min.diff.pct)
  plt <- hist(mono.de$avg_log2FC, breaks = 100, main = "Distribution of Log Fold Changes", xlab = "Log Fold Change", xlim = c(-4, 3))
  list(plt, mono.de)
}

.getDRGroundTruth <- function(res, cif_idx, noise, thres = 1) {
  paste0("gene", which(
    # rowSums(res_norm_2$giv$s[,c(1,33)]) > 0.5
    abs(res$giv$s[,cif_idx] %*% noise) > thres
  ))
}

.testDE <- function(res0, res_norm_2, pop, cif_idx, noise, thres = 1, plot.name = NULL) {
  ncells <- sum(res0$cell_meta$pop == pop)
  sample_label <- c(rep("1", ncells), rep("0", ncells))
  x <- res0$counts[,res0$cell_meta$pop == pop]
  colnames(x) <- paste0(colnames(x), "_0")

  if (!is.null(plot.name))
    pdf(plot.name)
  plt <- deDist(cbind(
    res_norm_2$counts[,res_norm_2$cell_meta$pop == pop], x
  ), sample_label, "1", "0")
  if (!is.null(plot.name))
    dev.off()

  mono.de <- plt[[2]]

  de.det <- paste0("gene", which(
    # rowSums(res_norm_2$giv$s[,c(1,33)]) > 0.5
    abs(res_norm_2$giv$s[,cif_idx] %*% noise) > thres
  ))

  int <- intersect(de.det, rownames(mono.de))

  TP <- length(int)
  FP <- nrow(mono.de) - TP
  FN <- length(de.det) - TP
  TN <- 2000 - TP - FP - FN
  message(sprintf("TP: %g, FP: %g, FN: %g, TN: %g", TP, FP, FN, TN))
  message(sprintf("Accuracy: %g", sum(TP + TN) / 2000))
  message(sprintf("Precision: %g", sum(TP) / (sum(TP) + sum(FP))))
  message(sprintf("F1: %g", 2 * sum(TP) / (2 * sum(TP) + sum(FP) + sum(FN))))

  list(TP = TP, FP = FP, FN = FN, TN = TN,
       Accuracy = sum(TP + TN) / 2000,
       Precision = sum(TP) / (sum(TP) + sum(FP)),
       F1 = 2 * sum(TP) / (2 * sum(TP) + sum(FP) + sum(FN)))
}

res0 <- .simCondEffect(2)
saveRDS(res0, "figures/de_genes/cond_res0_n0.5.rds")
res1 <- .simCondEffect2(2)
saveRDS(res1, "figures/de_genes/cond_res1_n0.5.rds")

# res0 <- readRDS("cond_res0_n0.5.rds")
# res1 <- readRDS("cond_res1_n0.5.rds")

precision_list <- numeric()
accuracy_list <- numeric()
f1_list <- numeric()

# DE histogram
lapply(1:4, function (i) {
  gt <- res0[[3]][[i]]
  ds_res <- .testDE(res0[[1]], res0[[2]], as.character(i),
                    gt[[1]], gt[[2]], 1,
                    sprintf("figures/de_genes/cond1_de_%g.pdf", i))
  precision_list <<- c(precision_list, ds_res$Precision)
  accuracy_list <<- c(accuracy_list, ds_res$Accuracy)
  f1_list <<- c(f1_list, ds_res$F1)
})

lapply(1:4, function (i) {
  gt <- res0[[3]][[i]]
  ds_res <- .testDE(res1[[1]], res1[[2]], as.character(i),
                    gt[[1]], gt[[2]], 1,
                    sprintf("figures/de_genes/cond2_de_%g.pdf", i))
  precision_list <<- c(precision_list, ds_res$Precision)
  accuracy_list <<- c(accuracy_list, ds_res$Accuracy)
  f1_list <<- c(f1_list, ds_res$F1)
})

# box plot
plt <- ggplot(data.frame(
  score = c(precision_list, accuracy_list, f1_list),
  metric = factor(c(rep("Precision", 8), rep("Accuracy", 8), rep("F1", 8)), levels = c("Precision", "Accuracy", "F1"))
), aes(x = metric, y = score)) +
  geom_boxplot() +
  geom_jitter(width = 0.2)
ggsave("figures/de_genes/cond_de_boxplot.pdf", plt, width = 4, height = 4)

ncells <- ncol(res0[[1]]$counts)

random_2ncells <- sample(1:(2*ncells), 2*ncells, replace = FALSE)
p1 <- plot_tsne(
  log2(cbind(res0[[1]]$counts, res0[[2]]$counts) + 1)[,random_2ncells],
  c(rep("Condition A", ncells), rep("Condition B", ncells))[random_2ncells],
  runPCA = TRUE,
  alpha = 0.25
)
p1 <- p1 + scale_alpha_manual(values = c(0.5, 0.5)) + theme_bw()
ggsave("figures/de_genes/cond_de_tsne.pdf", p1, width = 6, height = 6)

stop()
p <- plot_tsne(
  log2(cbind(res0[[1]]$counts, res0[[2]]$counts, res1[[1]]$counts, res1[[2]]$counts) + 1),
  c(rep("Sample 1", 2*ncells), rep("Sample 2", 2*ncells)),
  runPCA = TRUE,
  alpha = 0.25
)
p <- p + scale_alpha_manual(values = c(0.5, 0.5)) + theme_bw()
ggsave("figures/de_genes/cond_by_sample_tsne.pdf", p, width = 6, height = 6)

p <- plot_tsne(
  log2(cbind(res0[[1]]$counts, res0[[2]]$counts, res1[[1]]$counts, res1[[2]]$counts) + 1),
  c(rep("Condition A", ncells), rep("Condition B", ncells), rep("Condition A", ncells), rep("Condition B", ncells)),
  runPCA = TRUE,
  alpha = 0.25
)
p <- p + scale_alpha_manual(values = c(0.5, 0.5)) + theme_bw()
ggsave("figures/de_genes/cond_by_cond_tsne.pdf", p, width = 6, height = 6)


#
# plot_tsne(log2(cbind(res0[[1]]$counts, res0[[2]]$counts, res1[[1]]$counts, res1[[2]]$counts) + 1),
#           c(rep("A1", ncells), rep("B1", ncells), rep("A2", ncells), rep("B2", ncells)))


