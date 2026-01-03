# ==============================================================================
# Comprehensive Simulation Study (Extended) - B = 10000
# Reproducibility script used for the paper
# Outputs:
#   Tables  -> analysis/outputs/
#   Figures -> analysis/figures/
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# Load the *local* package code so agb() matches the repo version
if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("Install devtools first: install.packages('devtools')")
}
devtools::load_all(".")

# Output folders
out_tbl_dir <- file.path("analysis", "outputs")
out_fig_dir <- file.path("analysis", "figures")
dir.create(out_tbl_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_fig_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Global Configuration
# -----------------------------
sample_sizes <- c(10, 15, 20, 30, 50, 75, 100, 200)

M <- 500
B <- 10000
alpha <- 0.05
set.seed(123)

message("Running EXTENDED simulation with M = ", M, ", B = ", B, ", alpha = ", alpha, ", seed = 123")

# -----------------------------
# Standard bootstrap means
# -----------------------------
std_boot_means <- function(x, B = 10000L) {
  n <- length(x)
  idx <- sample.int(n, size = B * n, replace = TRUE)
  mat <- matrix(x[idx], nrow = B, ncol = n)
  rowMeans(mat)
}

ci_percentile <- function(boot_means, alpha = 0.05) {
  q <- quantile(boot_means, probs = c(alpha/2, 1 - alpha/2), names = FALSE)
  c(L = q[1], U = q[2])
}

# -----------------------------
# Data-generating scenarios
# -----------------------------
scenarios <- list(
  Normal = list(
    gen = function(n) rnorm(n, 0, 1),
    true_mean = function(n) 0
  ),
  MixtureSymmetric = list(
    gen = function(n) {
      u <- runif(n)
      ifelse(u < 0.5, rnorm(n, -2, 1), rnorm(n,  2, 1))
    },
    true_mean = function(n) 0
  ),
  SkewedLognormalCentered = list(
    gen = function(n) {
      mu <- 0; sigma <- 1
      y <- rlnorm(n, meanlog = mu, sdlog = sigma)
      y - exp(mu + 0.5 * sigma^2)  # centered to mean 0
    },
    true_mean = function(n) 0
  ),
  HeavyTail_t3 = list(
    gen = function(n) {
      rt(n, df = 3) / sqrt(3)  # scaled so Var approx 1
    },
    true_mean = function(n) 0
  ),
  Contaminated = list(
    gen = function(n) {
      u <- runif(n)
      ifelse(u < 0.9, rnorm(n, 0, 1), rnorm(n, 0, 10))
    },
    true_mean = function(n) 0
  )
)

# -----------------------------
# Run simulation
# -----------------------------
results <- data.frame()

for (sc_name in names(scenarios)) {
  message("Scenario: ", sc_name)

  gen_fun <- scenarios[[sc_name]]$gen
  theta_fun <- scenarios[[sc_name]]$true_mean

  for (n in sample_sizes) {
    message("  n = ", n)

    cp_std <- logical(M)
    cp_agb <- logical(M)

    rmse_std <- numeric(M)
    rmse_agb <- numeric(M)

    len_std <- numeric(M)
    len_agb <- numeric(M)

    theta <- theta_fun(n)

    for (m in seq_len(M)) {
      x <- gen_fun(n)

      # Standard bootstrap percentile CI
      bm_std <- std_boot_means(x, B = B)
      ci_std <- ci_percentile(bm_std, alpha = alpha)

      # AGB percentile CI (from package)
      agb_out <- agb(x, n_boot = B, alpha = alpha, seed = 123 + m)
      ci_agb <- c(L = unname(agb_out$ci[1]), U = unname(agb_out$ci[2]))

      # Coverage
      cp_std[m] <- (theta >= ci_std["L"] && theta <= ci_std["U"])
      cp_agb[m] <- (theta >= ci_agb["L"] && theta <= ci_agb["U"])

      # CI length
      len_std[m] <- ci_std["U"] - ci_std["L"]
      len_agb[m] <- ci_agb["U"] - ci_agb["L"]

      # RMSE of bootstrap-mean estimator
      rmse_std[m] <- (mean(bm_std) - theta)^2
      rmse_agb[m] <- (agb_out$mean_boot - theta)^2
    }

    results <- rbind(
      results,
      data.frame(
        Scenario = sc_name,
        SampleSize = n,
        Method = "Standard Bootstrap",
        CP = mean(cp_std),
        RMSE = sqrt(mean(rmse_std)),
        AvgLength = mean(len_std)
      ),
      data.frame(
        Scenario = sc_name,
        SampleSize = n,
        Method = "Adaptive Bootstrap (AGB)",
        CP = mean(cp_agb),
        RMSE = sqrt(mean(rmse_agb)),
        AvgLength = mean(len_agb)
      )
    )
  }
}

# -----------------------------
# Save tables (same filenames as before)
# -----------------------------
write.csv(results,
          file.path(out_tbl_dir, "Simulation_Results_Table_B10000_Extended.csv"),
          row.names = FALSE)

wide <- results %>%
  pivot_wider(names_from = Method, values_from = c(CP, RMSE, AvgLength)) %>%
  arrange(Scenario, SampleSize)

write.csv(wide,
          file.path(out_tbl_dir, "Simulation_Results_Wide_B10000_Extended.csv"),
          row.names = FALSE)

# -----------------------------
# Figures (Normal: main, others: supplementary)
# -----------------------------
normal_df <- results %>% filter(Scenario == "Normal")

p_cp <- ggplot(normal_df, aes(x = SampleSize, y = CP, color = Method, shape = Method)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "black", alpha = 0.6) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_x_continuous(breaks = sample_sizes) +
  ylim(0.85, 1.0) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Coverage Probability vs. Sample Size (Normal benchmark)",
    subtitle = paste0("Standard vs. AGB (Target CP = 0.95), M=", M, ", B=", B),
    x = "Sample Size (n)",
    y = "Coverage Probability"
  )

ggsave(file.path(out_fig_dir, "Figure_3_Coverage_Convergence.png"),
       plot = p_cp, width = 8, height = 6, dpi = 300)

p_len <- ggplot(normal_df, aes(x = factor(SampleSize), y = AvgLength, fill = Method)) +
  geom_col(position = "dodge", alpha = 0.85) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "top",
    legend.title = element_blank()
  ) +
  labs(
    title = "Average Confidence Interval Length (Normal benchmark)",
    subtitle = paste0("M=", M, ", B=", B, " (95% CI)"),
    x = "Sample Size (n)",
    y = "Average Length of 95% CI"
  )

ggsave(file.path(out_fig_dir, "Figure_4_Interval_Length.png"),
       plot = p_len, width = 8, height = 6, dpi = 300)

for (sc_name in setdiff(unique(results$Scenario), "Normal")) {
  df_sc <- results %>% filter(Scenario == sc_name)

  p_sc <- ggplot(df_sc, aes(x = SampleSize, y = CP, color = Method, shape = Method)) +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "black", alpha = 0.6) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    scale_x_continuous(breaks = sample_sizes) +
    ylim(0.80, 1.0) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = paste0("Coverage Probability vs. Sample Size (", sc_name, ")"),
      subtitle = paste0("Standard vs. AGB (Target CP = 0.95), M=", M, ", B=", B),
      x = "Sample Size (n)",
      y = "Coverage Probability"
    )

  out_name <- paste0("Figure_S_", sc_name, "_Coverage.png")
  ggsave(file.path(out_fig_dir, out_name),
         plot = p_sc, width = 8, height = 6, dpi = 300)
}

message("----------------------------------------------------------")
message("Success! Outputs saved in:")
message("Tables:  ", normalizePath(out_tbl_dir, winslash = "/", mustWork = FALSE))
message("Figures: ", normalizePath(out_fig_dir, winslash = "/", mustWork = FALSE))
message("----------------------------------------------------------")
