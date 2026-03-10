# ============================================================
#  MTH209 :: MODULE 10 — TESTING IN EXPONENTIAL FAMILIES
# ============================================================

#' UMP test for one-sided hypothesis in one-parameter exponential family
#'
#' @param distribution Character string: the exponential family distribution
#' @param null_hyp     Character string: H_0 (e.g. "theta <= theta_0")
#' @param alt_hyp      Character string: H_1 (e.g. "theta > theta_0")
#' @param alpha        Numeric: significance level
#' @return A message describing what will be computed
#' @export
expfam_test_onesided <- function(distribution, null_hyp, alt_hyp, alpha = 0.05) {
  cat("=== ONE-PARAMETER EXPONENTIAL FAMILY: ONE-SIDED TEST ===\n")
  cat(sprintf("Family  : %s\n", distribution))
  cat(sprintf("H_0     : %s\n", null_hyp))
  cat(sprintf("H_1     : %s\n", alt_hyp))
  cat(sprintf("alpha   : %.3f\n", alpha))
  cat("I will:\n")
  cat("  1. Write the family in canonical form: p(x|eta) = h(x) exp[eta T(x) - A(eta)]\n")
  cat("  2. Exploit MLR property in T(X) — automatic for exponential families\n")
  cat("  3. By Karlin-Rubin: Reject H_0 if T(X) > c  (or < c for left-sided)\n")
  cat("  4. Determine c from null distribution of T(X) at level alpha\n")
  cat("  5. Confirm this test is UMP (not just MP)\n")
  invisible(NULL)
}

#' UMP/UMPU tests for two-sided hypothesis in one-parameter exponential family
#'
#' @param distribution Character string: the exponential family distribution
#' @param null_hyp     Character string: H_0 (e.g. "theta = theta_0")
#' @param alt_hyp      Character string: H_1 (e.g. "theta != theta_0")
#' @param alpha        Numeric: significance level
#' @return A message describing what will be computed
#' @export
expfam_test_twosided <- function(distribution, null_hyp, alt_hyp, alpha = 0.05) {
  cat("=== ONE-PARAMETER EXPONENTIAL FAMILY: TWO-SIDED TEST ===\n")
  cat(sprintf("Family  : %s\n", distribution))
  cat(sprintf("H_0     : %s\n", null_hyp))
  cat(sprintf("H_1     : %s\n", alt_hyp))
  cat(sprintf("alpha   : %.3f\n", alpha))
  cat("I will:\n")
  cat("  1. Note: No UMP test exists for two-sided alternatives (by symmetry argument)\n")
  cat("  2. Restrict to unbiased tests -> seek UMPU test\n")
  cat("  3. Rejection region: C = { T < c1 } union { T > c2 }\n")
  cat("  4. Solve two equations for (c1, c2):\n")
  cat("     (a) E_{theta_0}[phi(X)] = alpha           (size constraint)\n")
  cat("     (b) E_{theta_0}[T(X) phi(X)] = alpha * E_{theta_0}[T(X)]  (unbiasedness)\n")
  cat("  5. This test is UMPU of size alpha for the two-sided hypothesis\n")
  invisible(NULL)
}

#' Tests with Neyman structure in multi-parameter exponential family
#'
#' @param distribution Character string: the multi-parameter exponential family
#' @param null_hyp     Character string: H_0 on the parameter of interest
#' @param nuisance     Character string: nuisance parameter(s)
#' @param alpha        Numeric: significance level
#' @return A message describing what will be computed
#' @export
expfam_test_neyman_structure <- function(distribution, null_hyp, nuisance, alpha = 0.05) {
  cat("=== MULTI-PARAMETER EXPONENTIAL FAMILY: NEYMAN STRUCTURE & SIMILAR TESTS ===\n")
  cat(sprintf("Family           : %s\n", distribution))
  cat(sprintf("H_0              : %s\n", null_hyp))
  cat(sprintf("Nuisance params  : %s\n", nuisance))
  cat(sprintf("alpha            : %.3f\n", alpha))
  cat("I will:\n")
  cat("  1. Partition parameters: (theta_1 [interest], theta_2 [nuisance])\n")
  cat("  2. Find the sufficient statistic for nuisance T_2(X)\n")
  cat("  3. Neyman structure: test has size alpha conditionally on T_2(X) = t_2\n")
  cat("  4. Construct the conditional test: reject if T_1 | T_2 = t_2 is extreme\n")
  cat("  5. The resulting test is UMPU similar of size alpha\n")
  invisible(NULL)
}

#' Likelihood Ratio Test (LRT)
#'
#' @param distribution Character string naming the distribution
#' @param null_hyp     Character string: H_0 (can be composite)
#' @param alt_hyp      Character string: H_1
#' @param alpha        Numeric: significance level
#' @return A message describing what will be computed
#' @export
test_likelihood_ratio <- function(distribution, null_hyp, alt_hyp, alpha = 0.05) {
  cat("=== LIKELIHOOD RATIO TEST (LRT) ===\n")
  cat(sprintf("Distribution : %s\n", distribution))
  cat(sprintf("H_0 : %s\n", null_hyp))
  cat(sprintf("H_1 : %s\n", alt_hyp))
  cat(sprintf("alpha        : %.3f\n", alpha))
  cat("I will:\n")
  cat("  1. Compute: Lambda(x) = sup_{theta in H_0} L(theta; x)\n")
  cat("                        / sup_{theta in Omega} L(theta; x)\n")
  cat("  2. Reject H_0 if Lambda(x) < c  (small Lambda = data supports H_1 more)\n")
  cat("  3. Find c: under H_0, -2 log Lambda ->  chi^2(r)  by Wilks' theorem\n")
  cat("     where r = dim(Omega) - dim(H_0)\n")
  cat("  4. Use chi^2 critical value: reject if -2 log Lambda > chi^2_{r, alpha}\n")
  cat("  5. Note: LRT is general but may not be UMPU — use when UMP/UMPU unavailable\n")
  invisible(NULL)
}
