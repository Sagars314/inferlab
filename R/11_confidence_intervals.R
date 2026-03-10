# ============================================================
#  MTH209 :: MODULE 11 — CONFIDENCE INTERVALS
# ============================================================

#' Construct a confidence interval using a pivotal quantity
#'
#' @param distribution Character string naming the distribution
#' @param parameter    Character string: parameter of interest theta
#' @param pivot        Character string: the pivotal quantity Q(X, theta)
#' @param conf_level   Numeric: confidence level (e.g. 0.95)
#' @return A message describing what will be computed
#' @export
ci_pivot <- function(distribution, parameter, pivot, conf_level = 0.95) {
  cat("=== CONFIDENCE INTERVAL VIA PIVOTAL QUANTITY ===\n")
  cat(sprintf("Distribution   : %s\n", distribution))
  cat(sprintf("Parameter      : %s\n", parameter))
  cat(sprintf("Pivot Q(X,theta): %s\n", pivot))
  cat(sprintf("Confidence     : %.0f%%\n", 100*conf_level))
  cat("I will:\n")
  cat("  1. Verify Q(X, theta) is a pivot: its distribution is free of theta\n")
  cat("  2. Find constants (a, b) such that P(a <= Q <= b) = conf_level\n")
  cat("  3. Invert the inequality to isolate theta: L(X) <= theta <= U(X)\n")
  cat("  4. Return the confidence interval [L(X), U(X)]\n")
  cat("  5. Interpret: in repeated sampling, (1-alpha)x100% of such intervals\n")
  cat("     will contain the true parameter\n")
  invisible(NULL)
}

#' Find the Shortest Expected Length Confidence Interval
#'
#' @param distribution Character string naming the distribution
#' @param parameter    Character string: parameter of interest
#' @param conf_level   Numeric: confidence level
#' @return A message describing what will be computed
#' @export
ci_shortest <- function(distribution, parameter, conf_level = 0.95) {
  cat("=== SHORTEST EXPECTED LENGTH CONFIDENCE INTERVAL ===\n")
  cat(sprintf("Distribution : %s\n", distribution))
  cat(sprintf("Parameter    : %s\n", parameter))
  cat(sprintf("Confidence   : %.0f%%\n", 100*conf_level))
  cat("I will:\n")
  cat("  1. Parameterise the family of valid CIs by choice of (a, b)\n")
  cat("     with P(a <= Q <= b) = conf_level\n")
  cat("  2. Express expected length E[U(X) - L(X)] in terms of (a, b)\n")
  cat("  3. Minimise expected length subject to the coverage constraint\n")
  cat("  4. For symmetric unimodal pivots: equal-tail CI = shortest\n")
  cat("  5. Return the optimal (a*, b*) and the minimum expected length\n")
  invisible(NULL)
}

#' Construct a UMA (Uniformly Most Accurate) Confidence Interval
#'
#' @param distribution Character string naming the distribution
#' @param parameter    Character string: parameter of interest
#' @param conf_level   Numeric: confidence level
#' @return A message describing what will be computed
#' @export
ci_uma <- function(distribution, parameter, conf_level = 0.95) {
  cat("=== UMA (UNIFORMLY MOST ACCURATE) CONFIDENCE INTERVAL ===\n")
  cat(sprintf("Distribution : %s\n", distribution))
  cat(sprintf("Parameter    : %s\n", parameter))
  cat(sprintf("Confidence   : %.0f%%\n", 100*conf_level))
  cat("I will:\n")
  cat("  1. Recall: a CI C(X) is UMA if it minimises P(theta' in C(X)) for all theta' != theta\n")
  cat("     (smallest probability of covering FALSE parameter values)\n")
  cat("  2. Duality with hypothesis testing: invert the UMP test\n")
  cat("  3. For each theta_0: C(X) = { theta_0 : UMP test of H_0:theta=theta_0 does NOT reject }\n")
  cat("  4. Show this CI is UMA among all (1-alpha) confidence sets\n")
  invisible(NULL)
}

#' Construct a UMAU (Uniformly Most Accurate Unbiased) Confidence Interval
#'
#' @param distribution Character string naming the distribution
#' @param parameter    Character string: parameter of interest
#' @param conf_level   Numeric: confidence level
#' @return A message describing what will be computed
#' @export
ci_umau <- function(distribution, parameter, conf_level = 0.95) {
  cat("=== UMAU (UNIFORMLY MOST ACCURATE UNBIASED) CI ===\n")
  cat(sprintf("Distribution : %s\n", distribution))
  cat(sprintf("Parameter    : %s\n", parameter))
  cat(sprintf("Confidence   : %.0f%%\n", 100*conf_level))
  cat("I will:\n")
  cat("  1. Recall unbiased CI: P(theta' in C(X) | theta) <= (1-alpha) for theta' != theta\n")
  cat("  2. Duality: invert the UMPU test to construct the UMAU CI\n")
  cat("  3. For two-sided CI: solve for (c1(t), c2(t)) making the interval UMAU\n")
  cat("  4. Show this CI is optimal within the class of unbiased confidence sets\n")
  cat("  5. Compare UMAU vs UMCI: UMAU is preferred when UMA does not exist\n")
  invisible(NULL)
}
