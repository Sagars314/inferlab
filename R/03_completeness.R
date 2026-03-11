#' Completeness Check via Monte Carlo
#'
#' Checks approximate completeness of a family: a statistic T is complete if
#' E_theta[g(T)] = 0 for all theta implies g = 0 a.e.
#' This function tests whether no non-trivial zero-mean function of T exists
#' by examining whether the distribution of T varies with theta.
#'
#' @param stat_fn Function: computes T from a sample
#' @param rfn Function rfn(n, theta): generates a sample from the model
#' @param theta_vals Numeric vector of theta values
#' @param n Sample size
#' @param B Monte Carlo replications
#'
#' @return List with distributional summaries and completeness assessment
#' @export
#'
#' @examples
#' rfn <- function(n, theta) rpois(n, lambda = theta)
#' stat_fn <- function(x) sum(x)
#' completeness_check(stat_fn, rfn, theta_vals = c(1, 2, 3), n = 20)
completeness_check <- function(stat_fn, rfn, theta_vals, n = 20, B = 2000) {
  dist_list <- lapply(theta_vals, function(th) {
    vals <- replicate(B, stat_fn(rfn(n, th)))
    list(theta = th, mean = mean(vals), var = var(vals),
         vals = vals)
  })

  cat("=== Completeness Assessment (Monte Carlo) ===\n\n")
  cat("Distribution of T(X) across theta values:\n")
  summary_df <- data.frame(
    theta = sapply(dist_list, `[[`, "theta"),
    E_T   = sapply(dist_list, `[[`, "mean"),
    Var_T = sapply(dist_list, `[[`, "var")
  )
  print(round(summary_df, 4))

  # A complete statistic's distribution must vary with theta
  mean_variation <- var(summary_df$E_T)
  cat(sprintf("\nVariation in E[T] across theta: %.4f\n", mean_variation))
  if (mean_variation > 1e-4)
    cat("=> T(X) is likely COMPLETE: its distribution changes with theta,\n",
        "   so no non-trivial zero-expectation function exists.\n")
  else
    cat("=> WARNING: E[T] does not vary with theta — T may NOT be complete.\n")

  invisible(list(summary = summary_df, dist_list = dist_list))
}

#' Basu's Theorem Demonstration
#'
#' Demonstrates Basu's theorem: if T is a complete sufficient statistic and
#' A is an ancillary statistic, then T and A are independent.
#'
#' Tests empirical independence between a complete sufficient statistic and
#' a proposed ancillary statistic via Monte Carlo simulation.
#'
#' @param rfn Function rfn(n, theta): generates a sample
#' @param suff_fn Function: computes complete sufficient statistic
#' @param anc_fn Function: computes proposed ancillary statistic
#' @param theta Reference theta value
#' @param n Sample size
#' @param B Monte Carlo replications
#'
#' @return List with correlation and independence test
#' @export
#'
#' @examples
#' # Normal(mu, 1): X_bar is complete sufficient, Range is ancillary
#' rfn <- function(n, theta) rnorm(n, mean = theta, sd = 1)
#' suff_fn <- function(x) mean(x)              # complete sufficient
#' anc_fn  <- function(x) diff(range(x))        # ancillary
#' basu_theorem_demo(rfn, suff_fn, anc_fn, theta = 2, n = 10, B = 2000)
basu_theorem_demo <- function(rfn, suff_fn, anc_fn, theta = 1, n = 20, B = 2000) {
  samples <- replicate(B, {
    x <- rfn(n, theta)
    c(T = suff_fn(x), A = anc_fn(x))
  })

  T_vals <- samples["T", ]
  A_vals <- samples["A", ]

  cor_TA <- cor(T_vals, A_vals)
  # Test independence via chi-squared on binned data
  T_bin <- cut(T_vals, breaks = 5)
  A_bin <- cut(A_vals, breaks = 5)
  chi_test <- tryCatch(chisq.test(table(T_bin, A_bin)), error = function(e) NULL)

  cat("=== Basu's Theorem Demonstration ===\n")
  cat(sprintf("theta = %.2f, n = %d, B = %d replications\n\n", theta, n, B))
  cat(sprintf("Correlation between T and A : %.4f\n", cor_TA))
  if (!is.null(chi_test)) {
    cat(sprintf("Chi-sq independence test    : X^2 = %.3f, df = %d, p = %.4f\n",
                chi_test$statistic, chi_test$parameter, chi_test$p.value))
    if (chi_test$p.value > 0.05)
      cat("=> Fail to reject independence (p > 0.05): consistent with Basu's theorem.\n")
    else
      cat("=> Evidence against independence: verify if T is complete sufficient.\n")
  }
  cat("\nBasu's Theorem: If T is complete sufficient and A is ancillary,\n")
  cat("then T _|_ A for all theta.\n")

  invisible(list(correlation = cor_TA, T_vals = T_vals, A_vals = A_vals,
                 chi_test = chi_test))
}

#' Rao-Blackwell Improvement
#'
#' Applies the Rao-Blackwell theorem to improve an unbiased estimator
#' by conditioning on a sufficient statistic, via Monte Carlo approximation.
#'
#' The Rao-Blackwell estimator: delta*(T) = E[delta(X) | T(X) = t]
#'
#' @param data Numeric vector of observations
#' @param estimator_fn Function: initial unbiased estimator of theta from data
#' @param suff_stat_fn Function: computes the sufficient statistic T(x)
#' @param rfn Function rfn(n, theta_hat): simulates from the model
#' @param n_cond Integer: Monte Carlo size for conditional expectation
#' @param tol Tolerance for matching T values
#'
#' @return List with original estimate, RB estimate, and variance comparison
#' @export
#'
#' @examples
#' # Improving indicator estimator for Poisson P(X=0)=exp(-lambda)
#' set.seed(42)
#' x <- rpois(15, lambda = 1.5)
#' # Estimator: I(X_1 = 0) is unbiased for e^{-lambda}
#' estimator_fn <- function(x) as.numeric(x[1] == 0)
#' suff_fn <- function(x) sum(x)
#' rfn <- function(n, lam) rpois(n, lambda = lam)
#' rao_blackwell(x, estimator_fn, suff_fn, rfn)
rao_blackwell <- function(data, estimator_fn, suff_stat_fn, rfn,
                           n_cond = 5000, tol = 0.5) {
  n <- length(data)
  delta_obs  <- estimator_fn(data)
  T_obs      <- suff_stat_fn(data)
  lambda_hat <- mean(data)   # crude estimate for simulation

  # Approximate E[delta(X) | T(X) = T_obs] by simulation
  delta_vals <- T_vals <- numeric(n_cond)
  for (i in seq_len(n_cond)) {
    x_sim <- rfn(n, lambda_hat)
    delta_vals[i] <- estimator_fn(x_sim)
    T_vals[i]     <- suff_stat_fn(x_sim)
  }

  # Condition on T ≈ T_obs
  idx_match  <- which(abs(T_vals - T_obs) <= tol)
  if (length(idx_match) < 10) {
    message("Too few matches in conditioning. Widening tolerance.")
    idx_match <- which(abs(T_vals - T_obs) <= 2 * tol)
  }

  delta_rb   <- mean(delta_vals[idx_match])

  cat("=== Rao-Blackwell Improvement ===\n")
  cat(sprintf("n = %d,  T(x) = %.4f\n\n", n, T_obs))
  cat(sprintf("Original estimator  delta(x)   = %.4f\n", delta_obs))
  cat(sprintf("Rao-Blackwell est.  delta*(T)  = %.4f\n", delta_rb))
  cat(sprintf("(Based on %d conditioning samples)\n", length(idx_match)))
  cat("\nRao-Blackwell guarantees: Var(delta*(T)) <= Var(delta(X))\n")
  cat("The RB estimator has equal or smaller variance than the original.\n")

  invisible(list(
    T_obs       = T_obs,
    delta_orig  = delta_obs,
    delta_rb    = delta_rb,
    n_matched   = length(idx_match)
  ))
}
