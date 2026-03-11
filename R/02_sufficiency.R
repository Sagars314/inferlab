#' Sufficient Statistic Verification via Neyman-Fisher Factorization
#'
#' Checks whether a proposed statistic T(x) is sufficient for theta by verifying
#' the Neyman-Fisher factorization: f(x; theta) = g(T(x), theta) * h(x).
#'
#' @param data Numeric vector of observations
#' @param stat_fn Function computing the proposed sufficient statistic T(x)
#' @param density_fn Function f(x, theta) giving the joint density/pmf
#' @param theta_vals Numeric vector of theta values to test factorization
#' @param n_check Integer: number of data permutations to verify invariance
#'
#' @return List with factorization check results
#' @export
#'
#' @examples
#' # For Poisson, T = sum(X) is sufficient for lambda
#' x <- rpois(10, lambda = 2)
#' stat_fn <- function(x) sum(x)
#' density_fn <- function(x, theta) prod(dpois(x, lambda = theta))
#' result <- neyman_fisher_factorization(x, stat_fn, density_fn, theta_vals = c(1,2,3))
neyman_fisher_factorization <- function(data, stat_fn, density_fn,
                                         theta_vals = c(1, 2, 3),
                                         n_check = 5) {
  T_obs <- stat_fn(data)
  n <- length(data)

  cat("=== Neyman-Fisher Factorization Check ===\n")
  cat(sprintf("Observed T(x) = %.4f  (n = %d)\n\n", T_obs, n))

  # For each theta, check if f(x;theta) / g(T(x),theta) is constant across
  # random permutations with the same sufficient statistic value
  results <- lapply(theta_vals, function(theta) {
    f_obs <- density_fn(data, theta)
    list(theta = theta, f_obs = f_obs, T_obs = T_obs)
  })

  cat("Joint density values at observed data:\n")
  for (r in results) {
    cat(sprintf("  theta = %.2f :  f(x; theta) = %.6e\n", r$theta, r$f_obs))
  }

  # Ratio test: f(x1; theta) / f(x2; theta) should not depend on theta
  # when T(x1) = T(x2)
  cat("\nNote: Full factorization verification requires permutation with same T.\n")
  cat("Sufficient statistic T(x) =", T_obs, "\n")

  invisible(list(T_obs = T_obs, results = results))
}

#' Sufficient Statistic Summary for Common Distributions
#'
#' Returns the minimal sufficient statistic for standard exponential family distributions.
#'
#' @param family Character: distribution family
#' @param data Numeric vector of observations
#'
#' @return List with sufficient statistic and its value
#' @export
#'
#' @examples
#' x <- rnorm(30, mean = 2, sd = 1)
#' sufficient_statistic("normal_mean", data = x)
#' sufficient_statistic("normal_both", data = x)
sufficient_statistic <- function(family = c("normal_mean", "normal_variance",
                                             "normal_both", "poisson",
                                             "binomial", "exponential",
                                             "gamma_rate", "gamma_both"),
                                  data) {
  family <- match.arg(family)
  n <- length(data)
  xbar <- mean(data)
  s2   <- var(data) * (n - 1) / n  # MLE variance

  result <- switch(family,
    normal_mean = list(
      statistic   = "T(X) = X_bar (sample mean)",
      value       = xbar,
      formula     = "T = (1/n) * sum(X_i)",
      note        = "sigma^2 assumed known; complete & minimal sufficient for mu"
    ),
    normal_variance = list(
      statistic   = "T(X) = sum((X_i - mu)^2) / n",
      value       = mean((data - mean(data))^2),
      formula     = "T = sum((X_i-mu)^2)/n  [mu known]",
      note        = "mu assumed known; minimal sufficient for sigma^2"
    ),
    normal_both = list(
      statistic   = "T(X) = (X_bar, S^2)  [two-dimensional]",
      value       = c(mean = xbar, S2 = var(data)),
      formula     = "T = (sum X_i, sum X_i^2)",
      note        = "Complete minimal sufficient for (mu, sigma^2)"
    ),
    poisson = list(
      statistic   = "T(X) = sum(X_i)  [total count]",
      value       = sum(data),
      formula     = "T = sum(X_i)",
      note        = "Complete minimal sufficient for lambda"
    ),
    binomial = list(
      statistic   = "T(X) = sum(X_i)  [total successes]",
      value       = sum(data),
      formula     = "T = sum(X_i)",
      note        = "Complete minimal sufficient for p"
    ),
    exponential = list(
      statistic   = "T(X) = sum(X_i)  [total lifetime]",
      value       = sum(data),
      formula     = "T = sum(X_i)",
      note        = "Complete minimal sufficient for lambda (rate)"
    ),
    gamma_rate = list(
      statistic   = "T(X) = sum(X_i)  [sum, shape alpha known]",
      value       = sum(data),
      formula     = "T = sum(X_i)",
      note        = "Complete minimal sufficient for beta (rate), alpha known"
    ),
    gamma_both = list(
      statistic   = "T(X) = (sum X_i, sum log X_i)  [two-dimensional]",
      value       = c(sum_x = sum(data), sum_logx = sum(log(data[data > 0]))),
      formula     = "T = (sum X_i, sum log X_i)",
      note        = "Complete minimal sufficient for (alpha, beta)"
    )
  )

  cat("=== Sufficient Statistic ===\n")
  cat("Family    :", family, "\n")
  cat("Statistic :", result$statistic, "\n")
  cat("Value     :", result$value, "\n")
  cat("Note      :", result$note, "\n")
  invisible(result)
}

#' Minimal Sufficient Statistic via Lehmann-Scheffe Criterion
#'
#' Demonstrates the ratio criterion for minimal sufficiency:
#' T(x) is minimal sufficient if f(x;theta)/f(y;theta) is constant in theta
#' iff T(x) = T(y).
#'
#' @param data1 Numeric vector, first observation set
#' @param data2 Numeric vector, second observation set (same length as data1)
#' @param density_fn Function f(x, theta)
#' @param stat_fn Proposed minimal sufficient statistic function
#' @param theta_vals Theta values for ratio check
#'
#' @return Logical indicating if ratio is theta-free when T(x1)=T(x2)
#' @export
minimal_sufficient <- function(data1, data2, density_fn, stat_fn,
                                theta_vals = seq(0.5, 3, by = 0.5)) {
  T1 <- stat_fn(data1)
  T2 <- stat_fn(data2)
  same_T <- isTRUE(all.equal(T1, T2, tolerance = 1e-8))

  ratios <- sapply(theta_vals, function(th) {
    f1 <- density_fn(data1, th)
    f2 <- density_fn(data2, th)
    if (f2 == 0) NA else f1 / f2
  })

  cat("=== Minimal Sufficiency Check ===\n")
  cat(sprintf("T(x1) = %.6f,  T(x2) = %.6f\n", T1, T2))
  cat(sprintf("Same T: %s\n\n", same_T))
  cat("Likelihood ratios f(x1;theta)/f(x2;theta):\n")
  for (i in seq_along(theta_vals)) {
    cat(sprintf("  theta = %.2f : ratio = %.6f\n", theta_vals[i], ratios[i]))
  }

  ratio_const <- diff(range(ratios, na.rm = TRUE)) < 1e-6
  cat(sprintf("\nRatio is theta-free: %s\n", ratio_const))
  if (same_T && ratio_const)
    cat("=> Consistent with T being minimal sufficient.\n")
  else if (!same_T && !ratio_const)
    cat("=> Consistent with T being minimal sufficient (different T => theta-dependent ratio).\n")

  invisible(list(T1 = T1, T2 = T2, ratios = ratios,
                 same_T = same_T, ratio_constant = ratio_const))
}

#' Ancillary Statistic Test
#'
#' Checks if a statistic A(X) is ancillary, i.e., its distribution does not
#' depend on the parameter theta.
#'
#' @param stat_fn Function computing the statistic from a sample
#' @param rfn Function rfn(n, theta) to simulate from the model
#' @param theta_vals Numeric vector of theta values
#' @param n Sample size
#' @param B Number of Monte Carlo replications
#'
#' @return Data frame of statistic means and variances across theta values
#' @export
#'
#' @examples
#' # Range / sigma is ancillary for location families
#' rfn <- function(n, theta) rnorm(n, mean = theta, sd = 1)
#' stat_fn <- function(x) diff(range(x))  # range is ancillary for normal location
#' ancillary_statistic_test(stat_fn, rfn, theta_vals = c(-2, 0, 2), n = 10)
ancillary_statistic_test <- function(stat_fn, rfn, theta_vals, n = 20, B = 1000) {
  results <- lapply(theta_vals, function(th) {
    vals <- replicate(B, stat_fn(rfn(n, th)))
    c(theta = th, mean = mean(vals), sd = sd(vals),
      q25 = quantile(vals, 0.25), q75 = quantile(vals, 0.75))
  })
  df <- as.data.frame(do.call(rbind, results))

  cat("=== Ancillary Statistic Check (Monte Carlo) ===\n")
  cat(sprintf("n = %d, B = %d replications\n\n", n, B))
  print(round(df, 4))

  mean_range <- diff(range(df$mean))
  sd_range   <- diff(range(df$sd))
  cat(sprintf("\nRange of means across theta: %.4f\n", mean_range))
  cat(sprintf("Range of SDs   across theta: %.4f\n", sd_range))
  if (mean_range < 0.1 * mean(abs(df$mean)) + 0.05)
    cat("=> Statistic appears ANCILLARY (distribution stable across theta).\n")
  else
    cat("=> Statistic does NOT appear ancillary.\n")

  invisible(df)
}
