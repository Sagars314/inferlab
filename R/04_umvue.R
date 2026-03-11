#' UMVUE for Normal Mean (sigma known)
#'
#' The UMVUE of mu for X_i ~ N(mu, sigma^2) is X_bar (sample mean).
#' By Lehmann-Scheffe: X_bar is a function of the complete sufficient statistic
#' and is unbiased, hence it is the UMVUE.
#'
#' @param x Numeric vector of observations
#' @param sigma Known standard deviation
#'
#' @return List with UMVUE estimate and properties
#' @export
#'
#' @examples
#' x <- rnorm(25, mean = 3, sd = 2)
#' umvue_normal_mean(x, sigma = 2)
umvue_normal_mean <- function(x, sigma = 1) {
  n     <- length(x)
  xbar  <- mean(x)
  se    <- sigma / sqrt(n)

  cat("=== UMVUE: Normal Mean (sigma known) ===\n")
  cat(sprintf("n = %d, sigma = %.4f\n", n, sigma))
  cat(sprintf("UMVUE(mu)  = X_bar = %.6f\n", xbar))
  cat(sprintf("Variance   = sigma^2/n = %.6f\n", sigma^2 / n))
  cat(sprintf("Std. Error = sigma/sqrt(n) = %.6f\n", se))
  cat("\nTheory: X_bar is the unique UMVUE by Lehmann-Scheffe theorem,\n")
  cat("since sum(X_i) is complete sufficient and X_bar = sum(X_i)/n is unbiased.\n")

  invisible(list(
    estimator = "X_bar",
    estimate  = xbar,
    variance  = sigma^2 / n,
    se        = se,
    n         = n,
    complete_suff_stat = sum(x)
  ))
}

#' UMVUE for Normal Variance (mu unknown)
#'
#' The UMVUE of sigma^2 for X_i ~ N(mu, sigma^2) is S^2 = sum((X_i-X_bar)^2)/(n-1).
#'
#' @param x Numeric vector of observations
#'
#' @return List with UMVUE estimate and properties
#' @export
#'
#' @examples
#' x <- rnorm(20, mean = 0, sd = 3)
#' umvue_normal_variance(x)
umvue_normal_variance <- function(x) {
  n   <- length(x)
  s2  <- var(x)  # (n-1) denominator
  xbar <- mean(x)

  # UMVUE of sigma^k and related quantities
  umvue_sigma2 <- s2
  umvue_sigma  <- sqrt(s2) * gamma((n-1)/2) / (gamma(n/2) / sqrt((n-1)/2))
  # More precisely:
  c4 <- sqrt(2/(n-1)) * gamma(n/2) / gamma((n-1)/2)
  umvue_sigma_corrected <- sqrt(s2) / c4

  cat("=== UMVUE: Normal Variance (mu unknown) ===\n")
  cat(sprintf("n = %d\n", n))
  cat(sprintf("UMVUE(sigma^2) = S^2 = %.6f\n", umvue_sigma2))
  cat(sprintf("UMVUE(sigma)   = S/c4 = %.6f  (c4 = %.6f)\n", umvue_sigma_corrected, c4))
  cat(sprintf("Variance of S^2 = 2*sigma^4/(n-1) [estimated: %.6f]\n",
              2 * umvue_sigma2^2 / (n - 1)))
  cat("\nTheory: (X_bar, S^2) is complete sufficient for (mu, sigma^2).\n")
  cat("S^2 = sum((X_i-X_bar)^2)/(n-1) is unbiased for sigma^2, hence UMVUE.\n")

  invisible(list(
    estimator            = "S^2",
    estimate_sigma2      = umvue_sigma2,
    estimate_sigma       = umvue_sigma_corrected,
    c4                   = c4,
    variance_of_s2       = 2 * umvue_sigma2^2 / (n - 1),
    n                    = n
  ))
}

#' UMVUE for Exponential Distribution
#'
#' For X_i ~ Exp(lambda) (rate parameterization), the UMVUE of:
#'   - lambda: does not have a simple form via T = sum(X)
#'   - 1/lambda (mean): T/n = X_bar
#'   - P(X > t): (1 - t/T)^{n-1} for T > t
#'
#' @param x Numeric vector of observations
#' @param estimand Character: "mean", "rate", "survival" or "cdf"
#' @param t0 Numeric: value for survival/CDF estimation (required for those)
#'
#' @return List with UMVUE estimate
#' @export
#'
#' @examples
#' x <- rexp(20, rate = 2)
#' umvue_exponential(x, estimand = "mean")
#' umvue_exponential(x, estimand = "survival", t0 = 0.3)
umvue_exponential <- function(x, estimand = c("mean", "rate", "survival", "cdf"), t0 = NULL) {
  estimand <- match.arg(estimand)
  n <- length(x)
  T_stat <- sum(x)
  xbar   <- mean(x)

  result <- switch(estimand,
    mean = {
      est <- xbar
      cat("=== UMVUE: Exponential Mean (1/lambda) ===\n")
      cat(sprintf("n = %d,  T = sum(X_i) = %.4f\n", n, T_stat))
      cat(sprintf("UMVUE(1/lambda) = X_bar = T/n = %.6f\n", est))
      cat("\nTheory: T = sum(X_i) ~ Gamma(n, lambda) is complete sufficient.\n")
      cat("E[T/n] = 1/lambda, so T/n is UMVUE for 1/lambda.\n")
      list(estimand = "1/lambda", estimate = est)
    },
    rate = {
      # UMVUE of lambda is (n-1)/T
      est <- (n - 1) / T_stat
      cat("=== UMVUE: Exponential Rate (lambda) ===\n")
      cat(sprintf("n = %d,  T = sum(X_i) = %.4f\n", n, T_stat))
      cat(sprintf("UMVUE(lambda) = (n-1)/T = %.6f\n", est))
      cat("\nTheory: E[(n-1)/T] = lambda since T ~ Gamma(n, lambda).\n")
      list(estimand = "lambda", estimate = est)
    },
    survival = {
      if (is.null(t0)) stop("Provide t0 for survival estimation.")
      if (T_stat <= t0) {
        cat("Warning: T <= t0, UMVUE not defined.\n")
        return(invisible(NULL))
      }
      est <- (1 - t0 / T_stat)^(n - 1)
      cat("=== UMVUE: P(X > t0) for Exponential ===\n")
      cat(sprintf("n = %d,  t0 = %.4f,  T = %.4f\n", n, t0, T_stat))
      cat(sprintf("UMVUE(P(X>t0)) = (1 - t0/T)^(n-1) = %.6f\n", est))
      list(estimand = sprintf("P(X > %.4f)", t0), estimate = est)
    },
    cdf = {
      if (is.null(t0)) stop("Provide t0 for CDF estimation.")
      if (T_stat <= t0) {
        cat("Warning: T <= t0, UMVUE not defined.\n")
        return(invisible(NULL))
      }
      est <- 1 - (1 - t0 / T_stat)^(n - 1)
      cat("=== UMVUE: P(X <= t0) for Exponential ===\n")
      cat(sprintf("n = %d,  t0 = %.4f,  T = %.4f\n", n, t0, T_stat))
      cat(sprintf("UMVUE(P(X<=t0)) = 1 - (1 - t0/T)^(n-1) = %.6f\n", est))
      list(estimand = sprintf("P(X <= %.4f)", t0), estimate = est)
    }
  )
  invisible(result)
}

#' UMVUE for Binomial Parameter Functions
#'
#' For X_i ~ Bin(m, p), the UMVUE of various functions of p.
#'
#' @param x Numeric vector of observations (successes in each trial)
#' @param m Integer: number of trials per observation
#' @param estimand Character: "p", "p_sq", "variance", or "p_k"
#' @param k Integer: power for p_k estimand
#'
#' @return List with UMVUE estimate
#' @export
#'
#' @examples
#' x <- rbinom(15, size = 10, prob = 0.4)
#' umvue_binomial(x, m = 10, estimand = "p")
umvue_binomial <- function(x, m = 1, estimand = c("p", "p_sq", "variance", "p_k"), k = 2) {
  estimand <- match.arg(estimand)
  n  <- length(x)
  T  <- sum(x)   # complete sufficient statistic
  N  <- n * m    # total trials

  result <- switch(estimand,
    p = {
      est <- T / N
      cat("=== UMVUE: Binomial p ===\n")
      cat(sprintf("n = %d observations, m = %d trials each, T = sum(X) = %d\n", n, m, T))
      cat(sprintf("UMVUE(p) = T/(nm) = %.6f\n", est))
      list(estimand = "p", estimate = est)
    },
    p_sq = {
      # UMVUE of p^2: T(T-1) / (N(N-1))
      est <- if (N < 2) NA else T * (T - 1) / (N * (N - 1))
      cat("=== UMVUE: Binomial p^2 ===\n")
      cat(sprintf("UMVUE(p^2) = T(T-1)/(N(N-1)) = %.6f\n", est))
      list(estimand = "p^2", estimate = est)
    },
    variance = {
      # Var(X) = mp(1-p); UMVUE
      p_hat <- T / N
      p2_hat <- T * (T - 1) / (N * (N - 1))
      est <- m * (p_hat - p2_hat)
      cat("=== UMVUE: Binomial Variance mp(1-p) ===\n")
      cat(sprintf("UMVUE(mp(1-p)) = %.6f\n", est))
      list(estimand = "mp(1-p)", estimate = est)
    },
    p_k = {
      # UMVUE of p^k: product form
      if (T < k || N < k) {
        cat("T or N < k; UMVUE = 0\n")
        return(invisible(list(estimand = paste0("p^",k), estimate = 0)))
      }
      est <- prod((T - 0:(k-1))) / prod((N - 0:(k-1)))
      cat(sprintf("=== UMVUE: Binomial p^%d ===\n", k))
      cat(sprintf("UMVUE(p^%d) = T(T-1)...(T-%d+1) / [N(N-1)...(N-%d+1)] = %.6f\n",
                  k, k, k, est))
      list(estimand = paste0("p^", k), estimate = est)
    }
  )
  invisible(result)
}

#' UMVUE for Poisson Distribution
#'
#' For X_i ~ Poisson(lambda), the UMVUE of various functions of lambda.
#'
#' @param x Numeric vector of observations
#' @param estimand Character: "lambda", "prob_zero", "prob_k", "lambda_sq"
#' @param k0 Integer: value k for P(X = k)
#'
#' @return List with UMVUE estimate
#' @export
#'
#' @examples
#' x <- rpois(20, lambda = 2.5)
#' umvue_poisson(x, estimand = "lambda")
#' umvue_poisson(x, estimand = "prob_zero")
umvue_poisson <- function(x, estimand = c("lambda", "prob_zero", "prob_k", "lambda_sq"), k0 = 0) {
  estimand <- match.arg(estimand)
  n <- length(x)
  T <- sum(x)

  result <- switch(estimand,
    lambda = {
      est <- T / n
      cat("=== UMVUE: Poisson lambda ===\n")
      cat(sprintf("n = %d, T = sum(X_i) = %d\n", n, T))
      cat(sprintf("UMVUE(lambda) = X_bar = T/n = %.6f\n", est))
      list(estimand = "lambda", estimate = est)
    },
    prob_zero = {
      # UMVUE of P(X=0) = e^{-lambda} is (1 - 1/n)^T
      est <- (1 - 1/n)^T
      cat("=== UMVUE: P(X=0) = e^{-lambda} ===\n")
      cat(sprintf("n = %d, T = %d\n", n, T))
      cat(sprintf("UMVUE(e^{-lambda}) = (1 - 1/n)^T = %.6f\n", est))
      list(estimand = "P(X=0)", estimate = est)
    },
    prob_k = {
      # UMVUE of P(X=k) = e^{-lambda} * lambda^k / k!
      # is choose(T, k) * (1/n)^k * (1-1/n)^(T-k) for T >= k
      if (T < k0) {
        cat("T < k0; UMVUE = 0\n")
        return(invisible(list(estimand = paste0("P(X=",k0,")"), estimate = 0)))
      }
      est <- choose(T, k0) * (1/n)^k0 * (1 - 1/n)^(T - k0)
      cat(sprintf("=== UMVUE: P(X=%d) ===\n", k0))
      cat(sprintf("UMVUE = C(T,%d) * (1/n)^%d * (1-1/n)^(T-%d) = %.6f\n",
                  k0, k0, k0, est))
      list(estimand = paste0("P(X=",k0,")"), estimate = est)
    },
    lambda_sq = {
      # UMVUE of lambda^2 is T(T-n) / n^2 = X_bar^2 - X_bar/n
      est <- T * (T - n) / n^2
      cat("=== UMVUE: lambda^2 ===\n")
      cat(sprintf("n = %d, T = %d\n", n, T))
      cat(sprintf("UMVUE(lambda^2) = T(T-n)/n^2 = %.6f\n", est))
      list(estimand = "lambda^2", estimate = est)
    }
  )
  invisible(result)
}

#' Lehmann-Scheffe Theorem Application
#'
#' Demonstrates the Lehmann-Scheffe theorem: if T is a complete sufficient
#' statistic and g(T) is unbiased for theta, then g(T) is the unique UMVUE.
#'
#' @param x Numeric vector of observations
#' @param family Character: distribution family
#' @param estimand Character: quantity to estimate
#'
#' @return Structured explanation and UMVUE value
#' @export
#'
#' @examples
#' x <- rnorm(25, mean = 2, sd = 1)
#' lehmann_scheffe(x, family = "normal_mean", estimand = "mu")
lehmann_scheffe <- function(x,
                             family = c("normal_mean", "poisson", "exponential",
                                        "binomial", "normal_both"),
                             estimand = "mu") {
  family <- match.arg(family)
  n <- length(x)

  cat("=== Lehmann-Scheffe Theorem ===\n\n")
  cat("Statement: Let T be a complete sufficient statistic.\n")
  cat("If g(T) is unbiased for g(theta), then g(T) is the unique UMVUE.\n\n")

  info <- switch(family,
    normal_mean = list(
      dist         = "X_i ~ N(mu, sigma^2), sigma^2 known",
      complete_ss  = "T = sum(X_i) [or equivalently X_bar]",
      estimand_str = "mu",
      g_T          = "g(T) = T/n = X_bar",
      estimate     = mean(x),
      unbiased_check = "E[X_bar] = mu  ✓"
    ),
    poisson = list(
      dist         = "X_i ~ Poisson(lambda)",
      complete_ss  = "T = sum(X_i)",
      estimand_str = "lambda",
      g_T          = "g(T) = T/n = X_bar",
      estimate     = mean(x),
      unbiased_check = "E[X_bar] = lambda  ✓"
    ),
    exponential = list(
      dist         = "X_i ~ Exp(lambda)",
      complete_ss  = "T = sum(X_i) ~ Gamma(n, lambda)",
      estimand_str = "1/lambda (mean)",
      g_T          = "g(T) = T/n = X_bar",
      estimate     = mean(x),
      unbiased_check = "E[T/n] = 1/lambda  ✓"
    ),
    binomial = list(
      dist         = "X_i ~ Bin(m, p)",
      complete_ss  = "T = sum(X_i)",
      estimand_str = "p",
      g_T          = "g(T) = T/(nm)",
      estimate     = mean(x),
      unbiased_check = "E[T/(nm)] = p  ✓"
    ),
    normal_both = list(
      dist         = "X_i ~ N(mu, sigma^2), both unknown",
      complete_ss  = "(sum X_i, sum X_i^2)  [two-dimensional]",
      estimand_str = "(mu, sigma^2)",
      g_T          = "g(T) = (X_bar, S^2)",
      estimate     = c(mu = mean(x), sigma2 = var(x)),
      unbiased_check = "E[X_bar] = mu, E[S^2] = sigma^2  ✓"
    )
  )

  cat("Distribution       :", info$dist, "\n")
  cat("Complete Suff. Stat:", info$complete_ss, "\n")
  cat("Estimand           :", info$estimand_str, "\n")
  cat("UMVUE = g(T)       :", info$g_T, "\n")
  cat("Unbiasedness       :", info$unbiased_check, "\n")
  cat("\nUMVUE estimate     :", info$estimate, "\n")
  cat("\nConclusion: By Lehmann-Scheffe, this estimator is the unique UMVUE.\n")

  invisible(info)
}

#' Hoeffding's U-Statistics
#'
#' Computes a U-statistic with a given kernel function h for estimating
#' a population parameter. U-statistics are UMVUE in non-parametric families.
#'
#' @param x Numeric vector of observations
#' @param kernel_fn Function of degree m: h(x_1, ..., x_m) symmetric kernel
#' @param m Integer: degree of the kernel
#' @param estimand Character description of what is being estimated
#'
#' @return U-statistic value with variance estimate
#' @export
#'
#' @examples
#' x <- rnorm(20, mean = 3, sd = 2)
#'
#' # U-statistic for population mean (m=1)
#' hoeffding_u_stat(x, kernel_fn = function(x) x, m = 1, estimand = "mu")
#'
#' # U-statistic for population variance (m=2 kernel)
#' kernel_var <- function(x1, x2) (x1 - x2)^2 / 2
#' hoeffding_u_stat(x, kernel_fn = kernel_var, m = 2, estimand = "sigma^2")
#'
#' # Wilcoxon kernel for P(X > 0) [m=1]
#' hoeffding_u_stat(x, kernel_fn = function(x) as.numeric(x > 0), m = 1,
#'                  estimand = "P(X > 0)")
hoeffding_u_stat <- function(x, kernel_fn, m = 1, estimand = "theta") {
  n <- length(x)
  if (n < m) stop("Need n >= m.")

  # Generate all subsets of size m
  idx_mat <- combn(n, m)
  n_terms <- ncol(idx_mat)

  h_vals <- apply(idx_mat, 2, function(idx) {
    args <- as.list(x[idx])
    do.call(kernel_fn, args)
  })

  U <- mean(h_vals)

  # Jackknife variance estimate
  jack_vals <- sapply(seq_len(n), function(i) {
    x_leave1 <- x[-i]
    if (length(x_leave1) < m) return(NA)
    idx2 <- combn(length(x_leave1), m)
    h2   <- apply(idx2, 2, function(idx) {
      args <- as.list(x_leave1[idx])
      do.call(kernel_fn, args)
    })
    mean(h2)
  })
  var_jack <- (n - 1) / n * sum((jack_vals - mean(jack_vals, na.rm = TRUE))^2, na.rm = TRUE)

  cat("=== Hoeffding U-Statistic ===\n")
  cat(sprintf("n = %d, kernel degree m = %d\n", n, m))
  cat(sprintf("Estimand     : %s\n", estimand))
  cat(sprintf("U-statistic  : %.6f\n", U))
  cat(sprintf("Number of (m-subsets): %d\n", n_terms))
  cat(sprintf("Jackknife SE : %.6f\n", sqrt(var_jack)))
  cat("\nTheory: U-statistics are UMVUE in non-parametric families\n")
  cat("(Hoeffding 1948). They are unbiased and have minimum variance\n")
  cat("among all unbiased estimators based on exchangeable observations.\n")

  invisible(list(
    U          = U,
    estimand   = estimand,
    m          = m,
    n_subsets  = n_terms,
    h_vals     = h_vals,
    var_jack   = var_jack,
    se_jack    = sqrt(var_jack)
  ))
}
