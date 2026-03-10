# ============================================================
#  MTH209 :: MODULE 2 — SUFFICIENCY
# ============================================================

#' Apply the Neyman-Fisher Factorisation Criterion
#'
#' @param distribution Character string naming the distribution
#' @param statistic     Character string describing the candidate statistic T(x)
#' @return A message describing what will be computed
#' @export
sufficiency_nf_factorisation <- function(distribution, statistic) {
  cat("=== NEYMAN-FISHER FACTORISATION CRITERION ===\n")
  cat(sprintf("Distribution : %s\n", distribution))
  cat(sprintf("Statistic    : T(x) = %s\n", statistic))
  cat("I will:\n")
  cat("  1. Write the joint likelihood L(theta; x)\n")
  cat("  2. Attempt to factorise as: L(theta; x) = g(T(x), theta) * h(x)\n")
  cat("  3. If factorisation holds -> conclude T(x) is sufficient for theta\n")
  cat("  4. If factorisation fails  -> conclude T(x) is NOT sufficient\n")
  invisible(NULL)
}

#' Check minimal sufficiency of a statistic
#'
#' @param distribution Character string naming the distribution
#' @param statistic     Character string describing the candidate statistic T(x)
#' @return A message describing what will be computed
#' @export
sufficiency_minimal <- function(distribution, statistic) {
  cat("=== MINIMAL SUFFICIENCY CHECK ===\n")
  cat(sprintf("Distribution : %s\n", distribution))
  cat(sprintf("Statistic    : T(x) = %s\n", statistic))
  cat("I will:\n")
  cat("  1. Confirm T(x) is sufficient via NF factorisation\n")
  cat("  2. Use Lehmann-Scheffe ratio criterion:\n")
  cat("     Check that L(theta;x)/L(theta;y) is free of theta <=> T(x) = T(y)\n")
  cat("  3. Return TRUE/FALSE: Is T(x) minimally sufficient?\n")
  invisible(NULL)
}

#' Check if a statistic is ancillary
#'
#' @param distribution Character string naming the distribution
#' @param statistic     Character string describing the candidate statistic
#' @return A message describing what will be computed
#' @export
sufficiency_ancillary <- function(distribution, statistic) {
  cat("=== ANCILLARY STATISTIC CHECK ===\n")
  cat(sprintf("Distribution : %s\n", distribution))
  cat(sprintf("Statistic    : A(x) = %s\n", statistic))
  cat("I will:\n")
  cat("  1. Derive the distribution of A(x)\n")
  cat("  2. Check whether that distribution depends on the parameter theta\n")
  cat("  3. Return TRUE/FALSE: Is A(x) ancillary?\n")
  cat("  4. If TRUE, note: ancillary + complete sufficient => independent (Basu)\n")
  invisible(NULL)
}
