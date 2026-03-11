#' Method of Moments Estimator (MOME)
#'
#' Estimates parameters by equating population moments to sample moments.
#' Solves the system: mu_k(theta) = m_k for k = 1, ..., p parameters.
#'
#' @param x Numeric vector of observations
#' @param family Character: distribution family
#' @param n_moments Integer: number of moments to use (default = number of params)
#'
#' @return List with MOM estimates and their standard errors
#' @export
#'
#' @examples
#' x <- rgamma(50, shape = 2, rate = 3)
#' mome(x, family = "gamma")
#'
#' x <- rnorm(40, mean = 3, sd = 2)
#' mome(x, family = "normal")
mome <- function(x, family = c("normal", "poisson", "exponential", "gamma",
                                "binomial", "beta", "uniform", "geometric"),
                  n_moments = NULL) {
  family <- match.arg(family)
  n <- length(x)

  # Sample moments
  m1 <- mean(x)
  m2 <- mean(x^2)
  m3 <- mean(x^3)
  mu2 <- mean((x - m1)^2)  # 2nd central moment

  result <- switch(family,
    normal = {
      mu_hat    <- m1
      sigma_hat <- sqrt(mu2)
      cat("=== MOM Estimator: Normal(mu, sigma^2) ===\n")
      cat(sprintf("n = %d\n", n))
      cat(sprintf("E[X] = mu           => mu_hat    = X_bar = %.6f\n", mu_hat))
      cat(sprintf("E[(X-mu)^2] = sigma^2 => sigma_hat = S_n   = %.6f\n", sigma_hat))
      cat(sprintf("  (True unbiased S = %.6f)\n", sd(x)))
      list(estimates = c(mu = mu_hat, sigma = sigma_hat),
           moments_used = 2)
    },
    poisson = {
      lambda_hat <- m1
      cat("=== MOM Estimator: Poisson(lambda) ===\n")
      cat(sprintf("E[X] = lambda => lambda_hat = X_bar = %.6f\n", lambda_hat))
      list(estimates = c(lambda = lambda_hat), moments_used = 1)
    },
    exponential = {
      lambda_hat <- 1 / m1
      cat("=== MOM Estimator: Exponential(lambda) ===\n")
      cat(sprintf("E[X] = 1/lambda => lambda_hat = 1/X_bar = %.6f\n", lambda_hat))
      list(estimates = c(lambda = lambda_hat), moments_used = 1)
    },
    gamma = {
      # E[X] = alpha/beta, Var[X] = alpha/beta^2
      # => alpha = m1^2 / mu2, beta = m1 / mu2
      alpha_hat <- m1^2 / mu2
      beta_hat  <- m1 / mu2
      cat("=== MOM Estimator: Gamma(alpha, beta) ===\n")
      cat(sprintf("n = %d,  X_bar = %.4f,  S^2 = %.4f\n", n, m1, mu2))
      cat(sprintf("E[X] = alpha/beta    => ]\n"))
      cat(sprintf("Var[X] = alpha/beta^2 =>\n"))
      cat(sprintf("alpha_hat = X_bar^2 / S^2 = %.6f\n", alpha_hat))
      cat(sprintf("beta_hat  = X_bar / S^2   = %.6f\n", beta_hat))
      list(estimates = c(alpha = alpha_hat, beta = beta_hat), moments_used = 2)
    },
    beta = {
      # E[X] = alpha/(alpha+beta), Var[X] = alpha*beta/[(alpha+beta)^2*(alpha+beta+1)]
      mu_hat  <- m1
      var_hat <- mu2
      common  <- mu_hat * (1 - mu_hat) / var_hat - 1
      alpha_hat <- mu_hat * common
      beta_hat  <- (1 - mu_hat) * common
      cat("=== MOM Estimator: Beta(alpha, beta) ===\n")
      cat(sprintf("alpha_hat = %.6f,  beta_hat = %.6f\n", alpha_hat, beta_hat))
      list(estimates = c(alpha = alpha_hat, beta = beta_hat), moments_used = 2)
    },
    uniform = {
      # Uniform(a, b): E[X] = (a+b)/2, Var[X] = (b-a)^2/12
      a_hat <- m1 - sqrt(3 * mu2)
      b_hat <- m1 + sqrt(3 * mu2)
      cat("=== MOM Estimator: Uniform(a, b) ===\n")
      cat(sprintf("a_hat = %.6f,  b_hat = %.6f\n", a_hat, b_hat))
      list(estimates = c(a = a_hat, b = b_hat), moments_used = 2)
    },
    binomial = {
      # Assuming n_trials given or estimating p from proportions
      p_hat <- m1 / round(max(x))
      cat("=== MOM Estimator: Binomial(m, p) ===\n")
      cat(sprintf("p_hat = X_bar / m_obs = %.6f\n", p_hat))
      list(estimates = c(p = p_hat), moments_used = 1)
    },
    geometric = {
      p_hat <- 1 / m1
      cat("=== MOM Estimator: Geometric(p) ===\n")
      cat(sprintf("E[X] = 1/p => p_hat = 1/X_bar = %.6f\n", p_hat))
      list(estimates = c(p = p_hat), moments_used = 1)
    }
  )

  cat(sprintf("\nSample size n = %d\n", n))
  cat("MOM Estimates:", paste(names(result$estimates), "=",
                               round(result$estimates, 6), collapse = ", "), "\n")
  invisible(result)
}

#' Minimum Mean Square Error Estimator
#'
#' Finds the estimator of the form c * T(X) that minimizes the MSE,
#' where T is a known unbiased estimator. The optimal constant is:
#' c* = E[theta * T] / E[T^2]
#'
#' For estimators of the form c * T where T is unbiased for theta:
#' MSE(cT) = c^2 * Var(T) + (c-1)^2 * theta^2
#' Minimized at c* = theta^2 / (theta^2 + Var(T))
#'
#' @param x Numeric vector of observations
#' @param family Character: distribution family
#' @param estimand Character: parameter being estimated
#' @param theta_prior Numeric: prior estimate or true value of theta (for c* computation)
#'
#' @return List with MinMSE estimator and comparison to UMVUE
#' @export
#'
#' @examples
#' x <- rnorm(20, mean = 3, sd = 1)
#' min_mse_estimator(x, family = "normal_mean", theta_prior = 3)
#'
#' x <- rexp(20, rate = 2)
#' min_mse_estimator(x, family = "exponential", theta_prior = 2)
min_mse_estimator <- function(x, family = c("normal_mean", "normal_var",
                                              "poisson", "exponential"),
                                theta_prior = NULL, estimand = NULL) {
  family <- match.arg(family)
  n <- length(x)

  info <- switch(family,
    normal_mean = {
      sigma  <- 1  # assume known or estimated
      T_val  <- mean(x)
      var_T  <- sigma^2 / n
      theta0 <- theta_prior %||% T_val
      c_star <- theta0^2 / (theta0^2 + var_T)
      est_minmse <- c_star * T_val
      list(T_name = "X_bar", T_val = T_val, var_T = var_T,
           c_star = c_star, estimate = est_minmse,
           umvue = T_val, estimand = "mu")
    },
    normal_var = {
      s2    <- var(x)
      # MLE = (n-1)/n * s2 = sum(Xi - Xbar)^2 / n
      # UMVUE = s2 (unbiased)
      # MinMSE: c* = n / (n+1) applied to s2
      c_star <- (n - 1) / (n + 1)
      T_val  <- s2
      est_minmse <- c_star * T_val
      sigma4_est <- s2^2
      var_T  <- 2 * sigma4_est / (n - 1)
      list(T_name = "S^2", T_val = T_val, var_T = var_T,
           c_star = c_star, estimate = est_minmse,
           umvue = s2, estimand = "sigma^2",
           note = "c* = (n-1)/(n+1) minimizes MSE for sigma^2")
    },
    poisson = {
      T_val  <- mean(x)
      theta0 <- theta_prior %||% T_val
      var_T  <- theta0 / n
      c_star <- theta0^2 / (theta0^2 + var_T)
      list(T_name = "X_bar", T_val = T_val, var_T = var_T,
           c_star = c_star, estimate = c_star * T_val,
           umvue = T_val, estimand = "lambda")
    },
    exponential = {
      # UMVUE of lambda is (n-1)/sum(X), mean is 1/lambda
      T_val  <- mean(x)  # UMVUE of 1/lambda
      theta0 <- theta_prior %||% (1/T_val)  # theta = lambda
      # For estimating lambda: UMVUE = (n-1)/sum(X)
      umvue_lambda <- (n - 1) / sum(x)
      # MinMSE of form c*(n-1)/sum(X):
      # MSE = Var(c*T) + Bias^2 minimized at c* = n/(n+2)
      c_star <- n / (n + 2)
      est_minmse <- c_star * umvue_lambda
      var_T <- theta0^2 / (n - 2)   # Var of (n-1)/sum(X) for lambda
      list(T_name = "(n-1)/sum(X)", T_val = umvue_lambda, var_T = var_T,
           c_star = c_star, estimate = est_minmse,
           umvue = umvue_lambda, estimand = "lambda",
           note = "MinMSE shrinks UMVUE by c* = n/(n+2)")
    }
  )

  # MSE comparison
  theta_est <- theta_prior %||% info$T_val
  bias_umvue    <- 0
  var_umvue     <- info$var_T
  mse_umvue     <- var_umvue

  bias_minmse   <- (info$c_star - 1) * theta_est
  var_minmse    <- info$c_star^2 * info$var_T
  mse_minmse    <- var_minmse + bias_minmse^2

  cat("=== Minimum MSE Estimator ===\n")
  cat(sprintf("Family: %s,  n = %d\n\n", family, n))
  cat(sprintf("Unbiased estimator T  = %s = %.6f\n", info$T_name, info$T_val))
  cat(sprintf("MinMSE estimator c*T  = %.4f * %.6f = %.6f\n",
              info$c_star, info$T_val, info$estimate))
  if (!is.null(info$note)) cat("Note:", info$note, "\n")
  cat("\nMSE Comparison (at theta =", round(theta_est, 4), "):\n")
  cat(sprintf("  UMVUE    : MSE = Var = %.8f\n", mse_umvue))
  cat(sprintf("  MinMSE   : MSE = %.8f  (Var=%.8f, Bias^2=%.8f)\n",
              mse_minmse, var_minmse, bias_minmse^2))
  if (mse_minmse < mse_umvue)
    cat(sprintf("  => MinMSE improves over UMVUE by %.2f%%\n",
                100 * (mse_umvue - mse_minmse) / mse_umvue))

  invisible(list(
    T_val       = info$T_val,
    c_star      = info$c_star,
    estimate    = info$estimate,
    umvue       = info$umvue,
    mse_umvue   = mse_umvue,
    mse_minmse  = mse_minmse
  ))
}
