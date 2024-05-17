devtools::load_all(".")
library(peakRAM)

set.seed(1)

data(GRN_params_100)

ncells <- c(500, 1000, 2000, 5000, 10000)
ngenes <- c(1000, 2000, 5000)

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
        GRN = GRN_params_100,
        speed.up = T,
        cif.sigma = 0.8,
        tree = Phyla3(),
        diff.cif.fraction = 0.8,
        discrete.cif = F,
        do.velocity = F,
      ))
    )
    res <- rbind(res, data.frame(
      ncells = n, ngenes = g, ram = res0$Total_RAM_Used_MiB, time = res0$Elapsed_Time_sec))
    gc()
  }
}

write.csv(res, "./time_curve.csv", row.names = F)

curve_data <- read.csv("figures/timings/time_curve.csv")
p <- ggplot(curve_data, aes(x = ncells, y = time, color = factor(ngenes))) +
  geom_point() +
  geom_line() +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "Number of cells", y = "Time (s)", color = "Number of genes")
ggsave("figures/timings/timing_non_spatial.pdf", p, width = 6, height = 3)

p <- ggplot(curve_data, aes(x = ncells, y = ram, color = factor(ngenes))) +
  geom_point() +
  geom_line() +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "Number of cells", y = "RAM (MiB)", color = "Number of genes")
ggsave("figures/timings/ram_non_spatial.pdf", p, width = 6, height = 3)