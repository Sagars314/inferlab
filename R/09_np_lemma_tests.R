# ============================================================
#  MTH209 :: MODULE 9 — NEYMAN-PEARSON LEMMA & OPTIMAL TESTS
# ============================================================

#' Apply the Neyman-Pearson Lemma
#'
#' @param null_dist Character string: distribution under H_0
#' @param alt_dist  Character string: distribution under H_1
#' @param alpha     Numeric: significance level
#' @return A message describing what will be computed
#' @export
np_lemma <- function(null_dist, alt_dist, alpha = 0.05) {
  cat("=== NEYMAN-PEARSON LEMMA ===\n")
  cat(sprintf("H_0 distribution : %s\n", null_dist))
  cat(sprintf("H_1 distribution : %s\n", alt_dist))
  cat(sprintf("Significance     : alpha = %.3f\n", alpha))
  cat("I will:\n")
  cat("  1. Form the likelihood ratio: Lambda(x) = f(x; theta_1) / f(x; theta_0)\n")
  cat("  2. Find critical value k such that P(Lambda(X) > k | H_0) = alpha\n")
  cat("  3. Define the Most Powerful (MP) test: reject H_0 if Lambda(X) > k\n")
  cat("  4. State the NP Lemma: this test has maximum power among all size-alpha tests\n")
  cat("  5. Handle randomisation if P(Lambda = k | H_0) > 0\n")
  invisible(NULL)
}

#' Construct a Most Powerful (MP) test
#'
#' @param null_dist Character string: distribution under H_0 (simple)
#' @param alt_dist  Character string: distribution under H_1 (simple)
#' @param alpha     Numeric: significance level
#' @return A message describing what will be computed
#' @export
test_mp <- function(null_dist, alt_dist, alpha = 0.05) {
  cat("=== MOST POWERFUL TEST ===\n")
  cat(sprintf("H_0 : %s  (simple)\n", null_dist))
  cat(sprintf("H_1 : %s  (simple)\n", alt_dist))
  cat(sprintf("alpha = %.3f\n", alpha))
  cat("I will:\n")
  cat("  1. Compute likelihood ratio Lambda(x) via NP Lemma\n")
  cat("  2. Simplify Lambda > k to an equivalent condition on a natural statistic T\n")
  cat("  3. Determine the critical value from the null distribution of T\n")
  cat("  4. State the MP test in terms of T\n")
  cat("  5. Compute the power of this MP test\n")
  invisible(NULL)
}

#' Construct a Uniformly Most Powerful (UMP) test
#'
#' @param null_hyp     Character string: H_0
#' @param alt_hyp      Character string: H_1 (composite)
#' @param distribution Character string naming the distributional family
#' @param alpha        Numeric: significance level
#' @return A message describing what will be computed
#' @export
test_ump <- function(null_hyp, alt_hyp, distribution, alpha = 0.05) {
  cat("=== UNIFORMLY MOST POWERFUL TEST ===\n")
  cat(sprintf("H_0 : %s\n", null_hyp))
  cat(sprintf("H_1 : %s  (composite)\n", alt_hyp))
  cat(sprintf("Distribution : %s\n", distribution))
  cat(sprintf("alpha = %.3f\n", alpha))
  cat("I will:\n")
  cat("  1. Check if the NP rejection region is the SAME for all theta_1 in H_1\n")
  cat("  2. If YES: the MP test is also UMP (uniform over all alternatives)\n")
  cat("  3. Construct the UMP test and its rejection region\n")
  cat("  4. Plot the power function beta(theta) over all theta\n")
  cat("  5. Note: UMP tests often exist for one-sided hypotheses, rarely for two-sided\n")
  invisible(NULL)
}

#' Construct an Unbiased Test and UMPU Test
#'
#' @param null_hyp     Character string: H_0
#' @param alt_hyp      Character string: H_1
#' @param distribution Character string naming the distributional family
#' @param alpha        Numeric: significance level
#' @return A message describing what will be computed
#' @export
test_umpu <- function(null_hyp, alt_hyp, distribution, alpha = 0.05) {
  cat("=== UNIFORMLY MOST POWERFUL UNBIASED TEST (UMPU) ===\n")
  cat(sprintf("H_0 : %s\n", null_hyp))
  cat(sprintf("H_1 : %s\n", alt_hyp))
  cat(sprintf("Distribution : %s\n", distribution))
  cat(sprintf("alpha = %.3f\n", alpha))
  cat("I will:\n")
  cat("  1. Recall: An unbiased test requires beta(theta) >= alpha for all theta in H_1\n")
  cat("  2. Show that an unrestricted UMP test may not exist (two-sided case)\n")
  cat("  3. Restrict to unbiased tests and find the UMPU within this class\n")
  cat("  4. Construct the UMPU rejection region (often involves two critical values)\n")
  cat("  5. Plot the power function confirming unbiasedness\n")
  invisible(NULL)
}

#' Test using Monotone Likelihood Ratio (MLR) property
#'
#' @param distribution Character string naming the distribution
#' @param statistic    Character string: the statistic with MLR property
#' @param null_hyp     Character string: H_0
#' @param alt_hyp      Character string: H_1
#' @param alpha        Numeric: significance level
#' @return A message describing what will be computed
#' @export
test_mlr <- function(distribution, statistic, null_hyp, alt_hyp, alpha = 0.05) {
  cat("=== MONOTONE LIKELIHOOD RATIO (MLR) TEST ===\n")
  cat(sprintf("Distribution : %s\n", distribution))
  cat(sprintf("Statistic    : T(X) = %s\n", statistic))
  cat(sprintf("H_0 : %s\n", null_hyp))
  cat(sprintf("H_1 : %s\n", alt_hyp))
  cat("I will:\n")
  cat("  1. Verify the MLR property: Lambda(x; theta_2, theta_1) non-decreasing in T(X)\n")
  cat("     for all theta_2 > theta_1\n")
  cat("  2. By Karlin-Rubin theorem: the test 'Reject H_0 if T > c' is UMP\n")
  cat("  3. Determine c from the null distribution at level alpha\n")
  cat("  4. Show the power function is non-decreasing in theta\n")
  invisible(NULL)
}
