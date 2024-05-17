devtools::load_all(".")
library(peakRAM)

set.seed(1)

lig_params <- data.frame(
  target    = c(101, 102),
  regulator = c(103, 104),
  effect    = c(10, 10)
)

data(GRN_params_100)

ncells <- c(500, 1000, 2000, 5000)
ngenes <- c(500, 1000, 2000)

res <- data.frame(
  ncells = integer(),
  ngenes = integer(),
  ram = numeric(),
  time = numeric()
)


for (n in ncells) {
  for (g in ngenes) {
    res0 <- peakRAM::peakRAM(
      sim_true_counts(list2(
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
    res <- rbind(res, data.frame(
      ncells = n, ngenes = g, ram = res0$Total_RAM_Used_MiB, time = res0$Elapsed_Time_sec))
    gc()
  }
}

stop()
write.csv(res, "./time_curve_sp.csv", row.names = F)

curve_data <- read.csv("./time_curve_sp.csv")
p <- ggplot(curve_data, aes(x = ncells, y = time, color = factor(ngenes))) +
  geom_point() +
  geom_line() +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "Number of cells", y = "Time (s)", color = "Number of genes")
ggsave("./timing_spatial.pdf", p, width = 6, height = 3)

p <- ggplot(curve_data, aes(x = ncells, y = ram, color = factor(ngenes))) +
  geom_point() +
  geom_line() +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "Number of cells", y = "RAM (MiB)", color = "Number of genes")
ggsave("./ram_spatial.pdf", p, width = 6, height = 3)

res_out <- list()
for (n in names(rs0)) {
  res_out[[n]] <- rs0[[n]]
}
saveRDS(res_out, sprintf("./sim_spatial_%dcells_%dgenes.rds", n, g))
