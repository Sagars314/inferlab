# ============================================================
#  MTH209 :: MODULE 7 — METHODS OF ESTIMATION
# ============================================================

#' Method of Moments Estimator (MoME)
#'
#' @param distribution Character string naming the distribution
#' @param parameters   Character vector of parameters to estimate
#' @param data         Optional numeric vector of observed data
#' @return A message describing what will be computed
#' @export
estimation_mome <- function(distribution, parameters, data = NULL) {
  cat("=== METHOD OF MOMENTS ESTIMATOR (MoME) ===\n")
  cat(sprintf("Distribution : %s\n", distribution))
  cat(sprintf("Parameters   : (%s)\n", paste(parameters, collapse=", ")))
  if (!is.null(data)) cat(sprintf("n            : %d observations supplied\n", length(data)))
  cat("I will:\n")
  cat("  1. Derive the population moments: mu_k = E[X^k], k = 1, 2, ...\n")
  cat("  2. Express each parameter in terms of population moments\n")
  cat("  3. Replace population moments with sample moments: m_k = (1/n) sum X_i^k\n")
  cat("  4. Solve the system of equations for the MoME estimators\n")
  cat("  5. State consistency: MoME is consistent under regularity (by LLN)\n")
  invisible(NULL)
}

#' Maximum Likelihood Estimator (MLE)
#'
#' @param distribution Character string naming the distribution
#' @param parameters   Character vector of parameters to estimate
#' @param data         Optional numeric vector of observed data
#' @return A message describing what will be computed
#' @export
estimation_mle <- function(distribution, parameters, data = NULL) {
  cat("=== MAXIMUM LIKELIHOOD ESTIMATOR (MLE) ===\n")
  cat(sprintf("Distribution : %s\n", distribution))
  cat(sprintf("Parameters   : (%s)\n", paste(parameters, collapse=", ")))
  if (!is.null(data)) cat(sprintf("n            : %d observations supplied\n", length(data)))
  cat("I will:\n")
  cat("  1. Write the log-likelihood: l(theta) = sum log f(x_i; theta)\n")
  cat("  2. Compute the score equation(s): dl/dtheta = 0\n")
  cat("  3. Solve for the MLE: theta_hat\n")
  cat("  4. Verify it is a maximum (second-order condition: d^2l/dtheta^2 < 0)\n")
  cat("  5. State asymptotic properties: theta_hat -> N(theta, I(theta)^{-1})\n")
  invisible(NULL)
}

#' Minimum Mean Squared Error (MinMSE) Estimator
#'
#' @param estimator_class Character string: class of estimators to search over
#' @param parameter       Character string: parameter of interest
#' @param distribution    Character string naming the distribution
#' @return A message describing what will be computed
#' @export
estimation_minmse <- function(estimator_class, parameter, distribution) {
  cat("=== MINIMUM MSE ESTIMATOR ===\n")
  cat(sprintf("Estimator class : %s\n", estimator_class))
  cat(sprintf("Parameter       : %s\n", parameter))
  cat(sprintf("Distribution    : %s\n", distribution))
  cat("I will:\n")
  cat("  1. Express MSE: MSE(T) = Var(T) + [Bias(T)]^2\n")
  cat("  2. Allow for biased estimators (no unbiasedness constraint)\n")
  cat("  3. Minimise MSE over the class of estimators\n")
  cat("  4. Show the bias-variance trade-off explicitly\n")
  cat("  5. Note: MinMSE estimator may be biased but has smaller MSE than UMVUE\n")
  invisible(NULL)
}
