library(ggplot2)
library(patchwork)

### Load the functions
source("Exhaustive_Procedures.R")
source("plot_generating_function.R")

# Generates Figure 3 of the paper "An exhaustive ADDIS principle for online FWER control"
ggsave("results/Figure3.pdf",
  plot = plot_generator(4, 1000, 6 / (pi^2 * ((1:1000)^2))), width = 7.64,
  height = 3.59
)

# Generates Figure 4 of the paper "An exhaustive ADDIS principle for online FWER control"
ggsave("results/Figure4.pdf",
  plot = plot_generator(4, 1000, (1 / (((1:1000) + 1) * log((1:1000) + 1)^2)) / 2.06227),
  width = 7.64, height = 3.59
)

# Generates Figure 5 of the paper "An exhaustive ADDIS principle for online FWER control"
ggsave("results/Figure5.pdf", plot = plot_generator(2, 10, 6 / (pi^2 * ((1:10)^2))), width = 7.64, height = 3.59)
