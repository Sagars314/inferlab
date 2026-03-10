# ============================================================
#  MTH209 :: MODULE 4 — UNBIASEDNESS & UMVUE
# ============================================================

#' Check basic unbiasedness of an estimator
#'
#' @param estimator    Character string: estimator T(X)
#' @param parameter    Character string: parameter being estimated (e.g. "mu", "sigma^2")
#' @param distribution Character string naming the distribution
#' @return A message describing what will be computed
#' @export
unbiased_check <- function(estimator, parameter, distribution) {
  cat("=== UNBIASEDNESS CHECK ===\n")
  cat(sprintf("Estimator    : T(X) = %s\n", estimator))
  cat(sprintf("Parameter    : %s\n", parameter))
  cat(sprintf("Distribution : %s\n", distribution))
  cat("I will:\n")
  cat("  1. Compute E_theta[T(X)] analytically\n")
  cat("  2. Check whether E_theta[T(X)] = theta (or g(theta)) for ALL theta\n")
  cat("  3. If not equal: compute Bias = E[T(X)] - theta\n")
  cat("  4. Return: Unbiased / Biased (with bias expression)\n")
  invisible(NULL)
}

#' Find a Locally Minimum Variance Unbiased Estimator (LMVUE)
#'
#' @param estimator    Character string: candidate estimator
#' @param parameter    Character string: parameter of interest
#' @param theta0       Numeric: the specific point theta_0 for local optimality
#' @param distribution Character string naming the distribution
#' @return A message describing what will be computed
#' @export
umvue_local <- function(estimator, parameter, theta0, distribution) {
  cat("=== LOCALLY MINIMUM VARIANCE UNBIASED ESTIMATOR (LMVUE) ===\n")
  cat(sprintf("Estimator    : %s\n", estimator))
  cat(sprintf("Parameter    : %s\n", parameter))
  cat(sprintf("Local point  : theta_0 = %s\n", theta0))
  cat(sprintf("Distribution : %s\n", distribution))
  cat("I will:\n")
  cat("  1. Confirm unbiasedness of the estimator\n")
  cat("  2. Find all unbiased estimators of zero (the 'U_0' class)\n")
  cat("  3. Minimise Var_theta0[T(X)] subject to unbiasedness at theta_0\n")
  cat("  4. Return the LMVUE and its variance at theta_0\n")
  invisible(NULL)
}

#' Find the Uniformly Minimum Variance Unbiased Estimator (UMVUE)
#'
#' @param parameter    Character string: parameter of interest (e.g. "theta", "theta^2")
#' @param distribution Character string naming the distribution
#' @return A message describing what will be computed
#' @export
umvue_find <- function(parameter, distribution) {
  cat("=== UNIFORMLY MINIMUM VARIANCE UNBIASED ESTIMATOR (UMVUE) ===\n")
  cat(sprintf("Estimating   : g(theta) = %s\n", parameter))
  cat(sprintf("Distribution : %s\n", distribution))
  cat("I will:\n")
  cat("  1. Find a complete sufficient statistic T(X) for theta\n")
  cat("  2. Strategy A — Direct: Express g(theta) as E[h(T)] and identify h\n")
  cat("  3. Strategy B — Rao-Blackwell: Start with ANY unbiased estimator,\n")
  cat("     condition on T(X) to get delta*(T) = E[estimator | T]\n")
  cat("  4. By Lehmann-Scheffe: delta*(T) is the unique UMVUE\n")
  cat("  5. Return the UMVUE formula and its variance\n")
  invisible(NULL)
}

#' Apply Lehmann-Scheffe's Theorem
#'
#' @param estimator    Character string: unbiased estimator that is a function of T
#' @param suff_stat    Character string: complete sufficient statistic T(X)
#' @param distribution Character string naming the distribution
#' @return A message describing what will be computed
#' @export
umvue_lehmann_scheffe <- function(estimator, suff_stat, distribution) {
  cat("=== LEHMANN-SCHEFFE THEOREM ===\n")
  cat(sprintf("Estimator (fn of T)  : %s\n", estimator))
  cat(sprintf("Complete suff. stat. : T(X) = %s\n", suff_stat))
  cat(sprintf("Distribution         : %s\n", distribution))
  cat("I will:\n")
  cat("  1. Verify T(X) is complete and sufficient\n")
  cat("  2. Verify the estimator is a function of T(X) only\n")
  cat("  3. Verify the estimator is unbiased for g(theta)\n")
  cat("  4. Conclude: this estimator IS the UNIQUE UMVUE for g(theta)\n")
  cat("  5. Note: uniqueness holds a.s. — no other unbiased fn of T can exist\n")
  invisible(NULL)
}
