paired_simil <- function(x, y, method = "spearman", margin = 1) {
    assertthat::assert_that(
        nrow(x) == nrow(y),
        ncol(x) == ncol(y)
    )
    n <- if (margin == 1) nrow(x) else ncol(x)

    out <- sapply(seq_len(n), function(i) {
        subx <- if (margin == 1) x[i, , drop = FALSE] else x[, i, drop = FALSE]
        suby <- if (margin == 1) y[i, , drop = FALSE] else y[, i, drop = FALSE]
        dynutils::calculate_similarity(subx, suby, method = method, margin = margin)[[1]]
    })

    names(out) <- if (margin == 1) rownames(x) else colnames(x)

    out
}

rna_velo_knn <- function(results, velocity, perplexity = 70, randseed = 0, raw = FALSE) {
  set.seed(randseed)
  counts_s <- results$counts
  pop <- results$cell_meta$pop
  depth <- results$cell_meta$depth
  
  counts_s_lg <- t(log2(counts_s + 1))
  
  if (is.null(results$velocity)) {
    stop("The result object is not produced in velocity mode.")
  }
  
  process_velocity <- function(v) {
    assertthat::assert_that(
      nrow(counts_s) == nrow(v),
      ncol(counts_s) == ncol(v)
    )
    
    future_counts_s <- counts_s + v
    future_counts_s[future_counts_s < 0] <- 0
    future_counts_s_lg <- t(log2(future_counts_s + 1))
    future_counts_s_lg - counts_s_lg
  }
  
  
  normalize_velocity <- function(v) {
    v_normalizer <- apply(v, 2, \(vi) vi^2) %>% rowSums() %>% sqrt()
    t(t(v) / v_normalizer)
  }
  
  if (raw) {
    return(
      paired_simil(velocity, results$velocity, method = "cosine")
    )
  }
  
  dist_obj <- dist(counts_s_lg)
  dist_mat <- as.matrix(dist_obj)
  n_cells <- nrow(dist_mat)
  k <- ceiling(n_cells / 50)
  
  v_knn <- process_velocity(velocity) %>%
    apply(2, \(vi)
      distMat.KernelKnn(dist_mat, TEST_indices = NULL,
                        weights_function = 'gaussian',
                        y = vi, k = k, regression = TRUE)
    ) %>%
    normalize_velocity()
  
  v_true_knn <- process_velocity(results$velocity) %>%
    apply(2, \(vi)
      distMat.KernelKnn(dist_mat, TEST_indices = NULL,
                        weights_function = 'gaussian',
                        y = vi, k = k, regression = TRUE)
    ) %>%
    normalize_velocity()
  
  sim <- paired_simil(v_knn, v_true_knn, method = "cosine")
  
  mean(sim)
}


p_grn <- (\(){
  configs <- list(
    "1" = "velocyto_constant_velocity",
    # "2" = "velocyto_constant_unspliced",
    "3" = "scvelo_stochastic",
    "4" = "scvelo_determinstic",
    "5" = "scvelo_dynamical"
    # "6" = "scvelo_dynamical_residuals"
  )
  
  datasets <- expand.grid(
    ngenes = c(100, 200, 500),
    ncells = c(500, 750, 1000),
    seed = c(1, 2, 3, 4)
  ) %>% split(., seq(nrow(.)))
  
  plot_data <- data.frame()
  for (ds in datasets) {
    print(ds)
    c(ngenes, ncells, seed) %<-% ds
    ds_name <- sprintf("grn_%dcells_%dgenes_%d", ncells, ngenes, seed)
    res <- readRDS(file.path("../bench/velo/data", ds_name, "res.rds"))
    
    for (conf_id in names(configs)) {
      conf_name <- configs[[conf_id]]
      res_fn <- file.path("../bench/velo/result_0", ds_name, paste0(conf_id, ".csv"))
      res_data <- read.csv(res_fn, header = T, row.names = 1) %>% t()
      res_data[is.na(res_data)] = 0
      siml <- rna_velo_knn(res, res_data, raw=F) 
      plot_data <- rbind(plot_data, data.frame(
        name = conf_name, ds = ds_name, value = siml
      ))
    }
  }
  
  plot_data
})()

p_grn_raw <- (\(){
  configs <- list(
    "1" = "velocyto_constant_velocity",
    # "2" = "velocyto_constant_unspliced",
    "3" = "scvelo_stochastic",
    "4" = "scvelo_determinstic",
    "5" = "scvelo_dynamical"
    # "6" = "scvelo_dynamical_residuals"
  )
  
  datasets <- expand.grid(
    ngenes = c(100, 200, 500),
    ncells = c(500, 750, 1000),
    seed = c(1, 2, 3, 4)
  ) %>% split(., seq(nrow(.)))
  
  plot_data <- data.frame()
  for (ds in datasets) {
    print(ds)
    c(ngenes, ncells, seed) %<-% ds
    ds_name <- sprintf("grn_%dcells_%dgenes_%d", ncells, ngenes, seed)
    res <- readRDS(file.path("../bench/velo/data", ds_name, "res.rds"))
    
    for (conf_id in names(configs)) {
      conf_name <- configs[[conf_id]]
      res_fn <- file.path("../bench/velo/result_0", ds_name, paste0(conf_id, ".csv"))
      res_data <- read.csv(res_fn, header = T, row.names = 1) %>% t()
      res_data[is.na(res_data)] = 0
      siml <- rna_velo_knn(res, res_data, raw=T) 
      plot_data <- rbind(plot_data, data.frame(
        name = conf_name, ds = ds_name, value = siml
      ))
    }
  }
  
  plot_data
})()

p_nogrn <- (\(){
  configs <- list(
    "1" = "velocyto_constant_velocity",
    # "2" = "velocyto_constant_unspliced",
    "3" = "scvelo_stochastic",
    "4" = "scvelo_determinstic",
    "5" = "scvelo_dynamical"
    # "6" = "scvelo_dynamical_residuals"
  )
  
  datasets <- expand.grid(
    ngenes = c(100, 200, 500),
    ncells = c(500, 750, 1000),
    seed = c(1, 2, 3, 4)
  ) %>% split(., seq(nrow(.)))
  
  plot_data <- data.frame()
  for (ds in datasets) {
    print(ds)
    c(ngenes, ncells, seed) %<-% ds
    ds_name <- sprintf("nogrn_%dcells_%dgenes_%d", ncells, ngenes, seed)
    res <- readRDS(file.path("../bench/velo/data", ds_name, "res.rds"))
    
    for (conf_id in names(configs)) {
      conf_name <- configs[[conf_id]]
      res_fn <- file.path("../bench/velo/result_0", ds_name, paste0(conf_id, ".csv"))
      res_data <- read.csv(res_fn, header = T, row.names = 1) %>% t()
      res_data[is.na(res_data)] = 0
      siml <- rna_velo_knn(res, res_data, raw=F) 
      plot_data <- rbind(plot_data, data.frame(
        name = conf_name, ds = ds_name, value = siml
      ))
    }
  }
  
  plot_data
})()

p_nogrn_raw <- (\(){
  configs <- list(
    "1" = "velocyto_constant_velocity",
    # "2" = "velocyto_constant_unspliced",
    "3" = "scvelo_stochastic",
    "4" = "scvelo_determinstic",
    "5" = "scvelo_dynamical"
    # "6" = "scvelo_dynamical_residuals"
  )
  
  datasets <- expand.grid(
    ngenes = c(100, 200, 500),
    ncells = c(500, 750, 1000),
    seed = c(1, 2, 3, 4)
  ) %>% split(., seq(nrow(.)))
  
  plot_data <- data.frame()
  for (ds in datasets) {
    print(ds)
    c(ngenes, ncells, seed) %<-% ds
    ds_name <- sprintf("nogrn_%dcells_%dgenes_%d", ncells, ngenes, seed)
    res <- readRDS(file.path("../bench/velo/data", ds_name, "res.rds"))
    
    for (conf_id in names(configs)) {
      conf_name <- configs[[conf_id]]
      res_fn <- file.path("../bench/velo/result_0", ds_name, paste0(conf_id, ".csv"))
      res_data <- read.csv(res_fn, header = T, row.names = 1) %>% t()
      res_data[is.na(res_data)] = 0
      siml <- rna_velo_knn(res, res_data, raw=T) 
      plot_data <- rbind(plot_data, data.frame(
        name = conf_name, ds = ds_name, value = siml
      ))
    }
  }
  
  plot_data
})()

write.csv(p_grn, "grn.csv", quote = FALSE)
write.csv(p_nogrn, "nogrn.csv", quote = FALSE)
write.csv(p_grn_raw, "grn_i.csv", quote = FALSE)
write.csv(p_nogrn_raw, "nogrn_i.csv", quote = FALSE)
