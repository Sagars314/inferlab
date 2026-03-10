# ============================================================
#  MTH209 :: MODULE 6 — INFORMATION INEQUALITIES & LOWER BOUNDS
# ============================================================

#' Compute Fisher Information
#'
#' @param distribution Character string naming the distribution
#' @param parameter    Character string: the parameter theta
#' @return A message describing what will be computed
#' @export
info_fisher <- function(distribution, parameter) {
  cat("=== FISHER INFORMATION ===\n")
  cat(sprintf("Distribution : %s\n", distribution))
  cat(sprintf("Parameter    : %s\n", parameter))
  cat("I will:\n")
  cat("  1. Compute I(theta) = E_theta[ (d/dtheta log f(X;theta))^2 ]\n")
  cat("  2. Alternative formula (under regularity): I(theta) = -E[d^2/dtheta^2 log f(X;theta)]\n")
  cat("  3. For n iid observations: I_n(theta) = n * I(theta)\n")
  cat("  4. Return a formula for I(theta) as a function of theta\n")
  invisible(NULL)
}

#' Compute the Fisher Information Matrix (multi-parameter)
#'
#' @param distribution Character string naming the distribution
#' @param parameters   Character vector of parameter names
#' @return A message describing what will be computed
#' @export
info_fisher_matrix <- function(distribution, parameters) {
  cat("=== FISHER INFORMATION MATRIX ===\n")
  cat(sprintf("Distribution : %s\n", distribution))
  cat(sprintf("Parameters   : (%s)\n", paste(parameters, collapse=", ")))
  cat("I will:\n")
  cat("  1. Compute the score vector: s(theta) = d/dtheta log f(X;theta)\n")
  cat("  2. Build I(theta) = E[ s(theta) s(theta)' ]  (outer product)\n")
  cat("  3. Alternative: I_{ij}(theta) = -E[ d^2/dtheta_i dtheta_j log f(X;theta) ]\n")
  cat(sprintf("  4. Return the %dx%d information matrix\n",
              length(parameters), length(parameters)))
  invisible(NULL)
}

#' Apply the Cramer-Rao Lower Bound
#'
#' @param estimator    Character string: estimator T(X)
#' @param parameter    Character string: g(theta) being estimated
#' @param distribution Character string naming the distribution
#' @return A message describing what will be computed
#' @export
info_cramer_rao <- function(estimator, parameter, distribution) {
  cat("=== CRAMER-RAO LOWER BOUND ===\n")
  cat(sprintf("Estimator    : T(X) = %s\n", estimator))
  cat(sprintf("g(theta)     : %s\n", parameter))
  cat(sprintf("Distribution : %s\n", distribution))
  cat("I will:\n")
  cat("  1. Verify regularity conditions (interchange derivative and integral)\n")
  cat("  2. Compute I(theta) — Fisher information\n")
  cat("  3. CRLB for unbiased estimator of g(theta):\n")
  cat("     Var(T) >= [g'(theta)]^2 / I(theta)\n")
  cat("  4. Compute the efficiency: e(T) = CRLB / Var(T)\n")
  cat("  5. If e(T) = 1 -> T is a UMVUE (achieves the bound)\n")
  invisible(NULL)
}

#' Apply the Hammersley-Chapman-Robbins (HCR) Inequality
#'
#' @param estimator    Character string: estimator T(X)
#' @param parameter    Character string: g(theta) being estimated
#' @param distribution Character string naming the distribution
#' @return A message describing what will be computed
#' @export
info_hcr_bound <- function(estimator, parameter, distribution) {
  cat("=== HAMMERSLEY-CHAPMAN-ROBBINS INEQUALITY ===\n")
  cat(sprintf("Estimator    : T(X) = %s\n", estimator))
  cat(sprintf("g(theta)     : %s\n", parameter))
  cat(sprintf("Distribution : %s\n", distribution))
  cat("I will:\n")
  cat("  1. Note: HCR does NOT require differentiability — stronger than CRLB\n")
  cat("  2. For any delta != 0: Var_theta(T) >= [g(theta+delta) - g(theta)]^2\n")
  cat("                                          / E_theta[(f(X;theta+delta)/f(X;theta) - 1)^2]\n")
  cat("  3. Optimise over delta to get the tightest bound\n")
  cat("  4. As delta->0 this recovers the Cramer-Rao bound (when CRLB applies)\n")
  invisible(NULL)
}

#' Compute the Bhattacharya System of Lower Bounds
#'
#' @param distribution Character string naming the distribution
#' @param parameter    Character string: parameter theta
#' @param order        Integer: order of the Bhattacharya bound (1 = CRLB, 2, 3, ...)
#' @return A message describing what will be computed
#' @export
info_bhattacharya <- function(distribution, parameter, order = 2) {
  cat("=== BHATTACHARYA SYSTEM OF LOWER BOUNDS ===\n")
  cat(sprintf("Distribution : %s\n", distribution))
  cat(sprintf("Parameter    : %s\n", parameter))
  cat(sprintf("Order        : %d\n", order))
  cat("I will:\n")
  cat("  1. Compute higher-order score functions:\n")
  cat("     rho_k(x;theta) = (d^k/dtheta^k) log f(x;theta),  k=1,...,order\n")
  cat("  2. Build the Bhattacharya information matrix B_order(theta)\n")
  cat("  3. Compute the order-k lower bound: Var(T) >= [g^(k)]' * B_k^{-1} * g^(k)\n")
  cat("  4. Note: Higher order = tighter bound. Order 1 = standard CRLB\n")
  invisible(NULL)
}

#' Information inequality for s-parameter exponential family
#'
#' @param distribution Character string naming the s-parameter exponential family
#' @param s            Integer: number of parameters
#' @return A message describing what will be computed
#' @export
info_exponential_family <- function(distribution, s) {
  cat("=== INFORMATION INEQUALITY: s-PARAMETER EXPONENTIAL FAMILY ===\n")
  cat(sprintf("Distribution : %s\n", distribution))
  cat(sprintf("Parameters s : %d\n", s))
  cat("I will:\n")
  cat("  1. Write the family in canonical form with natural parameters (eta_1,...,eta_s)\n")
  cat("  2. Show I(eta) = Var_eta(T(X)) = Hessian of A(eta)  (the cumulant function)\n")
  cat("  3. State the CRLB for this family in terms of A''(eta)\n")
  cat("  4. Note: In full-rank exp family, CRLB is always achieved by the MLE\n")
  invisible(NULL)
}
