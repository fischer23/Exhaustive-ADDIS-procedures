# Generates Figure 6 of the paper "An exhaustive ADDIS principle for online FWER control"

library(ggplot2)
library(patchwork)

### load the procedures
source("Exhaustive_Procedures.R")

### load data set
mydata <- read.csv("IMPC_data.csv")

# Number of hypotheses
n <- 5000

### Initialise Hyperparameters
gamma <- (1 / (((1:n) + 1) * log((1:n) + 1)^1.5)) / 2.47167 # 2.47167 is the approximated value such
# that the series equals 1
tau <- 0.8
lambda <- 0.16

w <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - (1:n - 1))
w[w == 0] <- 1
w <- matrix(gamma[w], n, n)
w[upper.tri(w) == 0] <- 0

# p-values
p <- mydata$Sex.P.Val

# Choosing different overall significance levels
alpha_vec <- c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4)

### Predefine vectors for number of rejections of the different procedures
n_rej_Alpha <- matrix(0, length(alpha_vec))
n_rej_Graph <- matrix(0, length(alpha_vec))
n_rej_EIGraph <- matrix(0, length(alpha_vec))
n_rej_Spending <- matrix(0, length(alpha_vec))
n_rej_ESpending <- matrix(0, length(alpha_vec))
n_rej_EGraph <- matrix(0, length(alpha_vec))

alpha_ind_Graph <- matrix(0, n, length(alpha_vec))
alpha_ind_EIGraph <- matrix(0, n, length(alpha_vec))

b <- 1 # Counter

for (alpha in alpha_vec) {
  ## Alpha-Spending
  alpha_ind <- alpha_spending(alpha, gamma, p, n)
  hypo_est <- alpha_ind >= p
  n_rej_Alpha[b] <- sum(hypo_est)

  ## ADDIS-Graph
  alpha_ind_Graph[, b] <- ADDIS_Graph(alpha, gamma, w, tau, lambda, 0, p, n)
  hypo_est <- alpha_ind_Graph[, b] >= p
  n_rej_Graph[b] <- sum(hypo_est)

  ## ADDIS-Spending
  alpha_ind <- ADDIS_Spending(alpha, gamma, tau, lambda, 0, p, n)
  hypo_est <- alpha_ind >= p
  n_rej_Spending[b] <- sum(hypo_est)

  ## E-ADDIS-Graph
  alpha_ind <- E_ADDIS_Graph(alpha, gamma, w, tau, lambda, p, n)
  hypo_est <- alpha_ind >= p
  n_rej_EGraph[b] <- sum(hypo_est)

  ## E-ADDIS-Spending
  alpha_ind <- E_ADDIS_Spending(alpha, gamma, tau, lambda, p, n)
  hypo_est <- alpha_ind >= p
  n_rej_ESpending[b] <- sum(hypo_est)

  ## EI-ADDIS-Graph
  alpha_ind_EIGraph[, b] <- EI_ADDIS_Graph(alpha, gamma, w, tau, lambda, p, n)
  hypo_est <- alpha_ind_EIGraph[, b] >= p
  n_rej_EIGraph[b] <- sum(hypo_est)

  b <- b + 1
}

### Create Plot
cols <- c("limegreen", "#f84f4f", "#e52c10", "#3a14f5", "#3dbeff", "cornflowerblue")

lab <- c("Alpha-Spending", "ADDIS-Graph", "ADDIS-Spending", "E-ADDIS-Graph", "E-ADDIS-Spending", "EI-ADDIS-Graph")

results_df <- data.frame(
  idx = seq(0.05, 0.4, 0.05), n_rej_Alpha, n_rej_Graph, n_rej_Spending, n_rej_EGraph,
  n_rej_ESpending, n_rej_EIGraph
)

plot_real <- ggplot(results_df, aes(idx)) +
  geom_line(aes(y = n_rej_Alpha, colour = "Alpha-Spending")) +
  geom_point(aes(y = n_rej_Alpha, colour = "Alpha-Spending", shape = "Alpha-Spending")) +
  geom_line(aes(y = n_rej_Graph, colour = "ADDIS-Graph")) +
  geom_point(aes(y = n_rej_Graph, colour = "ADDIS-Graph", shape = "ADDIS-Graph")) +
  geom_line(aes(y = n_rej_Spending, colour = "ADDIS-Spending")) +
  geom_point(aes(y = n_rej_Spending, colour = "ADDIS-Spending", shape = "ADDIS-Spending")) +
  geom_line(aes(y = n_rej_EGraph, colour = "E-ADDIS-Graph")) +
  geom_point(aes(y = n_rej_EGraph, colour = "E-ADDIS-Graph", shape = "E-ADDIS-Graph")) +
  geom_line(aes(y = n_rej_ESpending, colour = "E-ADDIS-Spending")) +
  geom_point(aes(y = n_rej_ESpending, colour = "E-ADDIS-Spending", shape = "E-ADDIS-Spending")) +
  geom_line(aes(y = n_rej_EIGraph, colour = "EI-ADDIS-Graph")) +
  geom_point(aes(y = n_rej_EIGraph, colour = "EI-ADDIS-Graph", shape = "EI-ADDIS-Graph")) +
  scale_shape_manual(
    name = "Procedure",
    values = c(
      "Alpha-Spending" = 8, "ADDIS-Graph" = 16, "ADDIS-Spending" = 4, "E-ADDIS-Graph" = 5,
      "E-ADDIS-Spending" = 6, "EI-ADDIS-Graph" = 17
    ),
    labels = c(
      "Alpha-Spending" = lab[1], "ADDIS-Graph" = lab[2], "ADDIS-Spending" = lab[3],
      "E-ADDIS-Graph" = lab[4], "E-ADDIS-Spending" = lab[5], "EI-ADDIS-Graph" = lab[6]
    )
  ) +
  scale_color_manual(
    name = "Procedure",
    values = c(
      "Alpha-Spending" = cols[1], "ADDIS-Graph" = cols[2], "ADDIS-Spending" = cols[3],
      "E-ADDIS-Graph" = cols[4], "E-ADDIS-Spending" = cols[5], "EI-ADDIS-Graph" = cols[6]
    ),
    labels = c(
      "Alpha-Spending" = lab[1], "ADDIS-Graph" = lab[2], "ADDIS-Spending" = lab[3],
      "E-ADDIS-Graph" = lab[4], "E-ADDIS-Spending" = lab[5], "EI-ADDIS-Graph" = lab[6]
    )
  ) +
  xlab("FWER level") +
  ylab("Number of rejections") +
  scale_x_continuous(breaks = seq(0.1, 0.35, 0.05), limits = c(0.05, 0.4), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(1000, 1600, 100), limits = c(1100, 1550), expand = c(0, 0)) +
  theme(
    panel.grid.major = element_line(color = "grey", size = 0.5, linetype = 1),
    panel.border = element_rect(colour = "black", fill = NA, size = 1), legend.text.align = 0
  )


ggsave("results/Figure6.pdf", plot = plot_real, width = 5.75, height = 3.65)
