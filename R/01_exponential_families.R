# ============================================================
#  MTH209 :: MODULE 1 — EXPONENTIAL FAMILIES
# ============================================================

#' Check if a distribution belongs to an exponential family
#'
#' @param distribution Character string naming the distribution (e.g. "Normal", "Poisson")
#' @param params Named list of parameters
#' @return A message describing what will be computed
#' @export
exp_family_check <- function(distribution, params = list()) {
  cat("=== EXPONENTIAL FAMILY CHECK ===\n")
  cat(sprintf("Distribution supplied: %s\n", distribution))
  cat("I will:\n")
  cat("  1. Write the pdf/pmf in the form: h(x) * exp[ eta(theta) * T(x) - A(theta) ]\n")
  cat("  2. Identify the natural parameter  eta(theta)\n")
  cat("  3. Identify the sufficient statistic T(x)\n")
  cat("  4. Identify the log-partition function A(theta)\n")
  cat("  5. Identify the carrier measure h(x)\n")
  cat("  6. Return TRUE/FALSE: Is this an exponential family?\n")
  invisible(NULL)
}

#' Express a distribution in canonical (natural) form
#'
#' @param distribution Character string naming the distribution
#' @return A message describing what will be computed
#' @export
exp_family_canonical <- function(distribution) {
  cat("=== CANONICAL FORM ===\n")
  cat(sprintf("Distribution: %s\n", distribution))
  cat("I will:\n")
  cat("  1. Re-parameterise using the natural parameter eta = eta(theta)\n")
  cat("  2. Express the family as: p(x|eta) = h(x) * exp[ eta' T(x) - A(eta) ]\n")
  cat("  3. Identify the natural parameter space: Eta = { eta : A(eta) < infinity }\n")
  cat("  4. State whether the canonical form is minimal\n")
  invisible(NULL)
}

#' Determine if an exponential family is full rank
#'
#' @param distribution Character string naming the distribution
#' @return A message describing what will be computed
#' @export
exp_family_full_rank <- function(distribution) {
  cat("=== FULL RANK CHECK ===\n")
  cat(sprintf("Distribution: %s\n", distribution))
  cat("I will:\n")
  cat("  1. Identify the dimension k of the sufficient statistic T(x)\n")
  cat("  2. Check that no component of T(x) is an affine function of the others\n")
  cat("  3. Check that the natural parameter space Eta contains an open set in R^k\n")
  cat("  4. Return TRUE/FALSE: Is this family full rank (also called 'steep')?\n")
  invisible(NULL)
}
