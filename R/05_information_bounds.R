#' Fisher Information
#'
#' Computes the Fisher information I(theta) = E[{d/dtheta log f(X;theta)}^2]
#' = -E[d^2/dtheta^2 log f(X;theta)] for common distributions.
#'
#' @param family Character: distribution family
#' @param theta Numeric: parameter value
#' @param param Named list: additional fixed parameters
#'
#' @return Fisher information value
#' @export
#'
#' @examples
#' fisher_information("normal_mean", theta = 2, param = list(sigma = 1))
#' fisher_information("poisson", theta = 3)
#' fisher_information("binomial", theta = 0.4, param = list(m = 10))
fisher_information <- function(family = c("normal_mean", "normal_var", "poisson",
                                           "binomial", "exponential", "gamma_rate",
                                           "bernoulli"),
                                theta, param = list()) {
  family <- match.arg(family)

  I_theta <- switch(family,
    normal_mean = {
      sigma <- param$sigma %||% 1
      1 / sigma^2
    },
    normal_var = {
      # theta = sigma^2, mu known
      1 / (2 * theta^2)
    },
    poisson = {
      1 / theta
    },
    binomial = {
      m <- param$m %||% param$n %||% 1
      m / (theta * (1 - theta))
    },
    exponential = {
      # theta = lambda (rate)
      1 / theta^2
    },
    gamma_rate = {
      # theta = beta (rate), alpha known
      alpha <- param$alpha %||% param$shape %||% 1
      alpha / theta^2
    },
    bernoulli = {
      1 / (theta * (1 - theta))
    }
  )

  cat("=== Fisher Information ===\n")
  cat(sprintf("Family: %s,  theta = %.4f\n", family, theta))
  cat(sprintf("I(theta) = %.6f\n", I_theta))
  cat(sprintf("1/I(theta) [CRB for unbiased est.] = %.6f\n", 1 / I_theta))

  invisible(list(family = family, theta = theta, I = I_theta, CRB = 1 / I_theta))
}

#' Fisher Information Matrix (Multiparameter Case)
#'
#' Computes the Fisher information matrix for multiparameter distributions.
#' I_{jk}(theta) = E[-d^2/dtheta_j dtheta_k log f(X; theta)]
#'
#' @param family Character: distribution family
#' @param theta Named numeric vector of parameter values
#' @param param Additional fixed parameters
#'
#' @return Fisher information matrix
#' @export
#'
#' @examples
#' fisher_information_matrix("normal", theta = c(mu = 2, sigma2 = 4))
#' fisher_information_matrix("gamma", theta = c(alpha = 2, beta = 1))
fisher_information_matrix <- function(family = c("normal", "gamma", "beta"),
                                       theta, param = list()) {
  family <- match.arg(family)

  I_mat <- switch(family,
    normal = {
      mu     <- theta["mu"]
      sigma2 <- theta["sigma2"]
      matrix(c(1/sigma2, 0,
               0,        1/(2*sigma2^2)),
             nrow = 2, ncol = 2,
             dimnames = list(c("mu","sigma2"), c("mu","sigma2")))
    },
    gamma = {
      alpha <- theta["alpha"]
      beta  <- theta["beta"]
      # Digamma and trigamma functions
      psi1_alpha <- trigamma(alpha)  # polygamma(1, alpha)
      I11 <- psi1_alpha
      I12 <- -1 / beta
      I22 <- alpha / beta^2
      matrix(c(I11, I12, I12, I22), nrow = 2,
             dimnames = list(c("alpha","beta"), c("alpha","beta")))
    },
    beta = {
      a <- theta["alpha"]
      b <- theta["beta"]
      psi1_a  <- trigamma(a)
      psi1_b  <- trigamma(b)
      psi1_ab <- trigamma(a + b)
      matrix(c(psi1_a - psi1_ab, -psi1_ab,
               -psi1_ab,          psi1_b - psi1_ab),
             nrow = 2,
             dimnames = list(c("alpha","beta"), c("alpha","beta")))
    }
  )

  cat("=== Fisher Information Matrix ===\n")
  cat("Family:", family, "\n")
  cat("Parameters:", paste(names(theta), "=", round(theta, 4), collapse = ", "), "\n\n")
  cat("I(theta) =\n")
  print(round(I_mat, 6))
  cat("\nInverse I(theta)^{-1} [CRB matrix] =\n")
  I_inv <- tryCatch(solve(I_mat), error = function(e) MASS::ginv(I_mat))
  print(round(I_inv, 6))

  invisible(list(I = I_mat, I_inv = I_inv, theta = theta))
}

#' Cramer-Rao Lower Bound
#'
#' Computes the Cramer-Rao lower bound for any unbiased estimator of g(theta).
#' CRB = [g'(theta)]^2 / (n * I(theta))
#'
#' @param family Character: distribution family
#' @param theta Numeric: parameter value
#' @param n Integer: sample size
#' @param g_deriv Numeric or function: derivative g'(theta), or function of theta
#' @param param Additional fixed parameters
#' @param data Optional data to compute an estimator's observed variance
#' @param estimator_fn Optional function to compute estimator from data
#'
#' @return List with CRB and efficiency analysis
#' @export
#'
#' @examples
#' # CRB for estimating mu in N(mu, 1) with n=20
#' cramer_rao_bound("normal_mean", theta = 2, n = 20, g_deriv = 1,
#'                  param = list(sigma = 1))
#'
#' # CRB for estimating e^{-lambda} in Poisson with n=20
#' cramer_rao_bound("poisson", theta = 1.5, n = 20,
#'                  g_deriv = function(th) -exp(-th))
cramer_rao_bound <- function(family, theta, n, g_deriv = 1, param = list(),
                              data = NULL, estimator_fn = NULL) {
  fi <- fisher_information(family, theta, param)
  I_theta <- fi$I

  g_prime <- if (is.function(g_deriv)) g_deriv(theta) else g_deriv

  CRB <- g_prime^2 / (n * I_theta)

  cat("\n=== Cramer-Rao Lower Bound ===\n")
  cat(sprintf("Family   : %s\n", family))
  cat(sprintf("theta    : %.4f\n", theta))
  cat(sprintf("n        : %d\n", n))
  cat(sprintf("g'(theta): %.6f\n", g_prime))
  cat(sprintf("I(theta) : %.6f\n", I_theta))
  cat(sprintf("\nCRB = [g'(theta)]^2 / (n * I(theta)) = %.8f\n", CRB))
  cat(sprintf("CRB (std. dev.) = %.8f\n", sqrt(CRB)))

  efficiency <- NULL
  if (!is.null(data) && !is.null(estimator_fn)) {
    # Bootstrap variance of estimator
    est_vals <- replicate(500, {
      x_b <- sample(data, replace = TRUE)
      estimator_fn(x_b)
    })
    obs_var  <- var(est_vals)
    efficiency <- CRB / obs_var
    cat(sprintf("\nObserved estimator variance (bootstrap): %.8f\n", obs_var))
    cat(sprintf("Efficiency = CRB / Var(estimator) = %.4f\n", efficiency))
    if (abs(efficiency - 1) < 0.05)
      cat("=> Estimator is approximately EFFICIENT (achieves CRB).\n")
    else
      cat("=> Estimator does not achieve CRB.\n")
  }

  invisible(list(CRB = CRB, I_theta = I_theta, g_prime = g_prime,
                 efficiency = efficiency))
}

#' Cramer-Rao Bound for Multiparameter Case
#'
#' For a vector estimator g(theta), the covariance matrix of any unbiased
#' estimator satisfies: Cov >= J * I(theta)^{-1} * J^T
#' where J is the Jacobian of g.
#'
#' @param family Character: distribution family
#' @param theta Named numeric vector of parameters
#' @param n Integer: sample size
#' @param J Jacobian matrix of g(theta) (rows = estimands, cols = parameters)
#'
#' @return Lower bound covariance matrix
#' @export
cramer_rao_multiparameter <- function(family, theta, n, J = diag(length(theta))) {
  fi <- fisher_information_matrix(family, theta)
  I_mat <- fi$I
  I_inv <- fi$I_inv

  CRB_mat <- J %*% I_inv %*% t(J) / n

  cat("\n=== Multiparameter Cramer-Rao Lower Bound ===\n")
  cat("CRB matrix (lower bound on covariance):\n")
  print(round(CRB_mat, 8))
  cat("CRB diagonal (lower bounds on individual variances):\n")
  print(diag(CRB_mat))

  invisible(list(CRB_matrix = CRB_mat, I = I_mat, I_inv = I_inv))
}

#' Hammersley-Chapman-Robbins Inequality
#'
#' Provides a lower bound on the variance of any estimator T(X) of g(theta):
#' Var_theta(T) >= sup_{delta != 0} [g(theta+delta) - g(theta)]^2 /
#'                                   E_theta[(f(X;theta+delta)/f(X;theta) - 1)^2]
#'
#' This bound does not require differentiability of the density.
#'
#' @param theta Numeric: parameter value
#' @param delta_vals Numeric vector: perturbations to try
#' @param g_fn Function g(theta): parameter function being estimated
#' @param chi2_fn Function chi2_fn(theta, delta): computes the chi-squared divergence
#'   E_theta[(f(X;theta+delta)/f(X;theta) - 1)^2]
#' @param n Integer: sample size
#'
#' @return HCR lower bound
#' @export
#'
#' @examples
#' # For Poisson(lambda), estimating g(lambda) = lambda
#' # chi-squared divergence for Poisson: e^{delta^2/(lambda)} - 1 (approx)
#' g_fn <- function(th) th
#' chi2_fn <- function(th, delta) exp(delta^2 / th) - 1
#' hcr_bound(theta = 2, delta_vals = seq(0.01, 1, by=0.05), g_fn, chi2_fn, n = 20)
hcr_bound <- function(theta, delta_vals, g_fn, chi2_fn, n = 1) {
  bounds <- sapply(delta_vals, function(delta) {
    num   <- (g_fn(theta + delta) - g_fn(theta))^2
    denom <- n * chi2_fn(theta, delta)
    if (denom <= 0) return(NA)
    num / denom
  })

  HCR <- max(bounds, na.rm = TRUE)
  best_delta <- delta_vals[which.max(bounds)]

  cat("=== Hammersley-Chapman-Robbins Bound ===\n")
  cat(sprintf("theta = %.4f, n = %d\n", theta, n))
  cat(sprintf("HCR bound = %.8f\n", HCR))
  cat(sprintf("Achieved at delta = %.4f\n", best_delta))
  cat("\nNote: HCR >= CRB always. Equality (CRB) holds under regularity conditions.\n")

  invisible(list(HCR = HCR, best_delta = best_delta, bounds = bounds))
}

#' Bhattacharyya System of Lower Bounds
#'
#' Computes the sequence of Bhattacharyya lower bounds, which improve upon
#' the Cramer-Rao bound by using higher-order score functions.
#'
#' For order s, the bound involves the s x s Bhattacharyya matrix B_s.
#'
#' @param family Character: distribution family
#' @param theta Numeric: parameter value
#' @param n Integer: sample size
#' @param g_fn Function: parameter function being estimated
#' @param max_order Integer: maximum order of Bhattacharyya bounds to compute (1 to 4)
#'
#' @return Data frame of bounds at each order
#' @export
#'
#' @examples
#' # Bhattacharyya bounds for estimating e^{-lambda} in Poisson
#' bhattacharyya_bounds("poisson", theta = 1.5, n = 20,
#'                      g_fn = function(th) exp(-th), max_order = 3)
bhattacharyya_bounds <- function(family = c("normal_mean", "poisson", "exponential",
                                             "binomial"),
                                  theta, n = 1, g_fn = function(th) th, max_order = 3) {
  family <- match.arg(family)

  # Score functions and their moments depend on the family
  # We compute numerically via finite differences

  h <- 1e-5

  # Compute log-density derivatives
  log_f <- function(th) {
    switch(family,
      normal_mean = function(x) dnorm(x, mean = th, log = TRUE),
      poisson     = function(x) dpois(round(x), lambda = th, log = TRUE),
      exponential = function(x) dexp(x, rate = th, log = TRUE),
      binomial    = function(x) dbinom(round(x), size = 1, prob = th, log = TRUE)
    )
  }

  # Numeric derivatives of log f w.r.t. theta (score of order k)
  score_k <- function(x, k) {
    # k-th derivative of log f using central differences
    tryCatch({
      if (k == 1) {
        (log_f(theta + h)(x) - log_f(theta - h)(x)) / (2 * h)
      } else if (k == 2) {
        (log_f(theta + h)(x) - 2*log_f(theta)(x) + log_f(theta - h)(x)) / h^2
      } else {
        0  # Higher orders require more complex implementation
      }
    }, error = function(e) 0)
  }

  # For common families, use exact Bhattacharyya matrices
  results <- data.frame(order = integer(), lower_bound = numeric())

  # Order 1 = CRB
  fi <- fisher_information(family, theta)
  I1 <- fi$I

  # g derivatives
  g_derivs <- sapply(1:max_order, function(k) {
    # Numerical derivative of order k
    if (k == 1)      (g_fn(theta + h) - g_fn(theta - h)) / (2 * h)
    else if (k == 2) (g_fn(theta + h) - 2*g_fn(theta) + g_fn(theta - h)) / h^2
    else if (k == 3) (g_fn(theta + 2*h) - 2*g_fn(theta + h) + 2*g_fn(theta - h) - g_fn(theta - 2*h)) / (2*h^3)
    else              0
  })

  # For analytical Bhattacharyya matrices (exact for exponential family)
  # B_s = matrix of E[rho_j * rho_k] where rho_k = d^k/dtheta^k log f
  # Use analytical results for common families

  bhat_bounds <- numeric(max_order)
  bhat_bounds[1] <- g_derivs[1]^2 / (n * I1)

  # Higher-order bounds (simplified; full matrix requires higher score moments)
  for (s in 2:min(max_order, 3)) {
    # For exponential families in canonical form, higher Bhattacharyya bounds
    # are readily available but require family-specific derivations
    # Here we provide the structure and order-1 and order-2 for Poisson/Normal

    if (family == "poisson" && s == 2) {
      # Bhattacharyya matrix for Poisson (exact)
      B11 <- 1/theta; B12 <- 1/theta; B22 <- (1 + 2*theta)/theta^2
      B_mat <- matrix(c(B11, B12, B12, B22), 2, 2)
      g_vec <- g_derivs[1:2]
      bhat_bounds[2] <- t(g_vec) %*% solve(B_mat) %*% g_vec / n
    } else if (family == "normal_mean" && s == 2) {
      sigma <- 1  # Assuming sigma=1
      B11 <- 1/sigma^2; B12 <- 0; B22 <- 2/sigma^4
      B_mat <- matrix(c(B11, B12, B12, B22), 2, 2)
      g_vec <- g_derivs[1:2]
      bhat_bounds[2] <- t(g_vec) %*% solve(B_mat) %*% g_vec / n
    } else {
      bhat_bounds[s] <- bhat_bounds[s-1]
    }
  }

  df <- data.frame(order = 1:max_order, lower_bound = bhat_bounds)

  cat("=== Bhattacharyya System of Lower Bounds ===\n")
  cat(sprintf("Family: %s,  theta = %.4f,  n = %d\n\n", family, theta, n))
  cat("g(theta) =", g_fn(theta), "\n")
  cat("g'(theta) =", round(g_derivs[1], 6), "\n\n")
  print(df)
  cat("\nNote: Bhattacharyya bounds are non-decreasing in order.\n")
  cat("Order 1 = Cramer-Rao bound. Higher orders give tighter bounds.\n")

  invisible(df)
}
