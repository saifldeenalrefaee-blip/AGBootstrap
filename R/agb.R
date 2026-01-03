#' Adaptive Generative Bootstrap (AGB)
#'
#' This function calculates confidence intervals for small or multimodal samples
#' using the Adaptive Generative Bootstrap method.
#'
#' @param data A numeric vector of the observed data.
#' @param n_boot Integer. Number of bootstrap iterations (default: 10000).
#' @param alpha Numeric. Significance level (default: 0.05 for 95% CI).
#' @param seed Integer. Random seed for reproducibility (optional).
#'
#' @return A list containing:
#' \item{ci}{The calculated confidence interval (Lower, Upper).}
#' \item{width}{The width of the confidence interval.}
#' \item{mean_boot}{The bootstrap mean estimate.}
#' \item{sample_size}{The size of the input sample.}
#'
#' @examples
#' data <- c(576, 635, 558, 578, 666, 580, 555, 661, 651, 605, 653, 575, 545, 572, 594)
#' result <- agb(data)
#' print(result$ci)
#'
#' @importFrom stats density approx rnorm quantile
#' @export
agb <- function(data, n_boot = 10000, alpha = 0.05, seed = NULL) {

  if (!is.null(seed)) set.seed(seed)
  n <- length(data)

  # 1. Pilot Density Estimation (Standard KDE)
  d_pilot <- density(data, bw = "nrd0")

  # Interpolate pilot density at data points
  f_pilot <- approx(d_pilot$x, d_pilot$y, xout = data)$y
  # Handle numerical zeros
  f_pilot[is.na(f_pilot) | f_pilot < 1e-10] <- 1e-10

  # Geometric mean of pilot density
  g <- exp(mean(log(f_pilot)))

  # 2. Abramson's Adaptive Bandwidth Factors
  # Lambda is proportional to 1/sqrt(f)
  lambda_factors <- (f_pilot / g) ^ -0.5

  # Clip extreme factors for stability (between 0.2 and 3.0)
  lambda_factors <- pmin(pmax(lambda_factors, 0.2), 3.0)

  # Calculate local bandwidths
  h_base <- d_pilot$bw
  h_local <- h_base * lambda_factors

  # 3. Generative Bootstrap Sampling
  indices <- sample(1:n, n_boot * n, replace = TRUE)
  centers <- data[indices]
  scales  <- h_local[indices]
  noise   <- rnorm(n_boot * n, mean = 0, sd = 1)

  # Generate new points: X* = X_i + h_i * epsilon
  resamples_vec <- centers + scales * noise

  # Reshape into matrix (rows = boot iterations, cols = sample size)
  resamples_matrix <- matrix(resamples_vec, nrow = n_boot, ncol = n)

  # Compute Statistics
  boot_means <- rowMeans(resamples_matrix)

  # Confidence Interval (percentile)
  ci <- quantile(boot_means, probs = c(alpha / 2, 1 - alpha / 2))

  list(
    ci = ci,
    width = as.numeric(ci[2] - ci[1]),
    mean_boot = mean(boot_means),
    sample_size = n
  )
}
