devtools::load_all(".")
library(peakRAM)

set.seed(1)

lig_params <- data.frame(
  target    = c(101, 102),
  regulator = c(103, 104),
  effect    = c(10, 10)
)

data(GRN_params_100)

n <- 5000
g <- 500

res <- NULL

time_res <- peakRAM::peakRAM(
  res <<- sim_true_counts(list2(
    num.cells = n,
    num.genes = g,
    num.cifs = 20,
    threads = 1,
    GRN = GRN_params_100,
    speed.up = T,
    cif.sigma = 0.8,
    tree = Phyla5(),
    diff.cif.fraction = 0.8,
    discrete.cif = F,
    do.velocity = F,
    cci = list(
      params = lig_params,
      max.neighbors = 4,
      start.layer = n,
      cell.type.interaction = "random",
      step.size = 1
    )
  ))
)

res$time <- time_res

res_out <- list()
for (i in names(res)) {
  if (i == "cif" || i == "giv" || i == "kinetic_params") next
  res_out[[i]] <- res[[i]]
}

saveRDS(res_out, sprintf("./sim_spatial_%dcells_%dgenes.rds", n, g))

