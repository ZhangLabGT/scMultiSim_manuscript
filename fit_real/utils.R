library(patchwork)

compare_real_dist <- \(..., labels = NULL){
  datasets <- list(...)
  xorder <- labels

  lib.size <- \(ct) apply(ct, 2, function(pc) sum(pc))
  zero.cnt <- \(ct) apply(ct, 2, function(pc) sum(pc == 0) / length(pc))
  zero.cnt.G <- \(ct) apply(ct, 1, function(pc) sum(pc == 0) / length(pc))
  mean.cnt.G <- \(ct) apply(ct, 1, function(pc) mean(pc))
  var.cnt.G <- \(ct) apply(ct, 1, function(pc) var(pc))

  subplot_box <- function (plt_title, fn, ylim = NULL) {
    df <- NULL
    for (i in seq_along(datasets)) {
      df <- rbind(
        df,
        data.frame(x = fn(datasets[[i]]), group = labels[i])
      )
    }
    df$group <- factor(df$group, levels = xorder)
    # boxplot(df$x ~ df$group, ylab="", xlab="", ylim = ylim)
    # title(plt_title)
    p <- ggplot(df, aes(group, x)) + geom_boxplot() + ggtitle(plt_title) +
        theme(axis.title.x=element_blank(), axis.title.y = element_blank())
    p
  }

  subplot_zero_vs_mean <- function (plt_title) {
    df <- NULL
    for (i in seq_along(datasets)) {
      df <- rbind(
        df,
        data.frame(
          y = zero.cnt.G(datasets[[i]]),
          x = mean.cnt.G(datasets[[i]]),
          group = labels[i]
        )
      )
    }
    df$group <- factor(df$group, levels = xorder)
    p <- ggplot(df, aes(x, y, colour = group)) + geom_point() + ggtitle(plt_title) +
      theme(axis.title.x=element_blank(), axis.title.y = element_blank())
    p
  }

  subplot_box("Lib size (per cell)", lib.size) +
    subplot_box("Zero count (per cell)", zero.cnt, ylim = c(0, 1)) +
    subplot_box("Zero count (per gene)", zero.cnt.G, ylim = c(0, 1)) +
    subplot_box("Mean count (per gene)", mean.cnt.G) +
    subplot_box("Var count (per gene)", var.cnt.G) +
    subplot_zero_vs_mean("Zero count vs Mean count (per gene)") +
    plot_layout(nrow = 1, ncol = 6)
}

gg_cor_hist <- function(counts) {
  counts <- log(counts + 1)
  cor_list <- numeric()
  for (i in seq(nrow(counts))) {
    for (j in seq_len(i-1)) {
      cor_list <- c(cor_list, cor(counts[i,], counts[j,]))
    }
  }
  hist(cor_list, breaks = 100)
}


plot_gg_heatmap <- function (counts, row_order = TRUE, filename = NULL) {
  ngenes <- nrow(counts)
  m <- cor(t(counts), method = "pearson")
  assert_that(nrow(m) == ngenes)
  colors <- unique(c(seq(-1,0.2,length=100),seq(-0.2,0.2,length=100),seq(0.2,1,length=100)))
  my_palette <- colorRampPalette(c("blue", "white", "red"))(n = length(colors) - 1)

  if (!is.null(filename)) png(filename, width = 800, height = 800)
  p <- heatmap.2(m,
                 col=my_palette, breaks=colors,
                 dendrogram = "none",
                 distfun = dist, hclustfun = hclust, key = TRUE, trace = "none",
                 cexRow = 1, cexCol = 1,
                 Rowv=row_order, Colv= if (is.null(row_order)) FALSE else "Rowv",
                 symm=F,symkey=F,symbreaks=T, scale="none",
                 # RowSideColors = gene_module_color_vector[1:num_GRN_genes],
                 # ColSideColors = gene_module_color_vector[1:num_GRN_genes],
                 main = 'Gene-gene correlation')
  if (!is.null(filename)) dev.off()
  p
}

plot_cc_heatmap <- function (counts, cell_types, filename = NULL) {
  ncells <- ncol(counts)
  color_vec <- c("red", "green", "blue", "yellow")
  m <- cor(counts, method = "pearson")
  assert_that(nrow(m) == ncells)
  m[is.na(m)] <- 0

  png(filename, width = 800, height = 800)
  heatmap.2(m,
            scale = "none", Rowv = TRUE, Colv = TRUE, dendrogram = "none",
            distfun = dist, hclustfun = hclust, key = TRUE, trace = "none",
            cexRow = 1, cexCol = 1,
            RowSideColors = color_vec[cell_types],
            ColSideColors = color_vec[cell_types],
            col = bluered(75), main = 'Cell-cell correlation')
  dev.off()
}
