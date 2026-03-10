# ============================================================
#  MTH209 :: MODULE 8 — BASIC CONCEPTS IN HYPOTHESIS TESTING
# ============================================================

#' Define and classify a statistical hypothesis
#'
#' @param null_hyp  Character string: the null hypothesis H_0
#' @param alt_hyp   Character string: the alternative hypothesis H_1
#' @return A message describing what will be computed
#' @export
test_define_hypothesis <- function(null_hyp, alt_hyp) {
  cat("=== HYPOTHESIS DEFINITION & CLASSIFICATION ===\n")
  cat(sprintf("H_0 : %s\n", null_hyp))
  cat(sprintf("H_1 : %s\n", alt_hyp))
  cat("I will:\n")
  cat("  1. Classify H_0: Simple (single point) or Composite (set of values)?\n")
  cat("  2. Classify H_1: Simple or Composite?\n")
  cat("  3. Identify the type of test:\n")
  cat("     - One-sided (H_1: theta > theta_0  OR  theta < theta_0)\n")
  cat("     - Two-sided (H_1: theta != theta_0)\n")
  invisible(NULL)
}

#' Analyse the critical region of a test
#'
#' @param test_statistic  Character string: the test statistic T(X)
#' @param critical_region Character string: description of the rejection region C
#' @param distribution    Character string naming the distribution under H_0
#' @return A message describing what will be computed
#' @export
test_critical_region <- function(test_statistic, critical_region, distribution) {
  cat("=== CRITICAL REGION ANALYSIS ===\n")
  cat(sprintf("Test statistic   : T(X) = %s\n", test_statistic))
  cat(sprintf("Critical region  : C = { %s }\n", critical_region))
  cat(sprintf("Distribution H_0 : %s\n", distribution))
  cat("I will:\n")
  cat("  1. Identify the sample space Omega = C union C^c\n")
  cat("  2. Define: Reject H_0 if T(X) in C\n")
  cat("  3. Define: Fail to reject H_0 if T(X) not in C\n")
  cat("  4. Visualise the rejection region on the distribution\n")
  invisible(NULL)
}

#' Compute Type I error, Type II error, size and power
#'
#' @param null_hyp        Character string: H_0
#' @param alt_hyp         Character string: H_1
#' @param test_statistic  Character string: T(X)
#' @param critical_region Character string: rejection region
#' @param distribution    Character string naming the distribution
#' @return A message describing what will be computed
#' @export
test_errors_and_power <- function(null_hyp, alt_hyp, test_statistic,
                                   critical_region, distribution) {
  cat("=== TYPE I/II ERRORS, SIZE & POWER ===\n")
  cat(sprintf("H_0 : %s\n", null_hyp))
  cat(sprintf("H_1 : %s\n", alt_hyp))
  cat(sprintf("Test statistic  : T(X) = %s\n", test_statistic))
  cat(sprintf("Rejection region: C = { %s }\n", critical_region))
  cat("I will:\n")
  cat("  1. Type I error  (alpha): P(T in C | H_0 true) — 'false positive'\n")
  cat("  2. Type II error (beta) : P(T not in C | H_1 true) — 'false negative'\n")
  cat("  3. Size of the test: sup_{theta in H_0} P(T in C | theta)\n")
  cat("  4. Power function: beta(theta) = P(T in C | theta) for theta in H_1\n")
  cat("  5. Power at a specific alternative: 1 - beta\n")
  invisible(NULL)
}
