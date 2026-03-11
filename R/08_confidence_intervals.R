#' Pivot-Based Confidence Interval
#'
#' Constructs a confidence interval using a pivotal quantity Q(X, theta)
#' whose distribution is free of unknown parameters.
#'
#' @param x Numeric vector of observations
#' @param family Character: distribution family
#' @param conf_level Numeric: confidence level (e.g., 0.95)
#' @param param Named list: additional known parameters
#'
#' @return List with confidence interval and pivotal quantity description
#' @export
#'
#' @examples
#' x <- rnorm(20, mean = 3, sd = 2)
#' pivot_ci(x, family = "normal_mean", conf_level = 0.95, param = list(sigma = 2))
#'
#' x <- rexp(25, rate = 1.5)
#' pivot_ci(x, family = "exponential", conf_level = 0.95)
pivot_ci <- function(x, family = c("normal_mean", "normal_both", "exponential",
                                    "poisson", "binomial", "gamma_rate"),
                      conf_level = 0.95, param = list()) {
  family <- match.arg(family)
  n      <- length(x)
  alpha  <- 1 - conf_level

  result <- switch(family,
    normal_mean = {
      sigma <- param$sigma %||% sd(x)
      known_sigma <- !is.null(param$sigma)
      xbar  <- mean(x)
      if (known_sigma) {
        # Pivot: Z = sqrt(n)(X_bar - mu)/sigma ~ N(0,1)
        z <- qnorm(1 - alpha/2)
        ci <- c(xbar - z * sigma/sqrt(n), xbar + z * sigma/sqrt(n))
        pivot_str <- "Z = sqrt(n)(X_bar - mu)/sigma ~ N(0,1)"
      } else {
        # Pivot: T = sqrt(n)(X_bar - mu)/S ~ t(n-1)
        t_crit <- qt(1 - alpha/2, df = n - 1)
        ci <- c(xbar - t_crit * sd(x)/sqrt(n), xbar + t_crit * sd(x)/sqrt(n))
        pivot_str <- "T = sqrt(n)(X_bar - mu)/S ~ t(n-1)"
      }
      list(ci = ci, pivot = pivot_str, estimand = "mu",
           point_est = xbar)
    },
    normal_both = {
      xbar  <- mean(x)
      s     <- sd(x)
      t_c   <- qt(1 - alpha/2, df = n - 1)
      chi_l <- qchisq(alpha/2, df = n - 1)
      chi_u <- qchisq(1 - alpha/2, df = n - 1)
      ci_mu    <- c(xbar - t_c * s/sqrt(n), xbar + t_c * s/sqrt(n))
      ci_sigma2 <- c((n-1)*s^2 / chi_u, (n-1)*s^2 / chi_l)
      list(ci_mu = ci_mu, ci_sigma2 = ci_sigma2,
           pivot_mu = "T = sqrt(n)(X_bar-mu)/S ~ t(n-1)",
           pivot_sigma2 = "(n-1)S^2/sigma^2 ~ chi^2(n-1)",
           point_est = c(mu = xbar, sigma2 = s^2))
    },
    exponential = {
      # Pivot: 2*lambda*sum(X) ~ chi^2(2n)
      T_stat <- sum(x)
      chi_l  <- qchisq(alpha/2, df = 2*n)
      chi_u  <- qchisq(1 - alpha/2, df = 2*n)
      ci <- c(chi_l / (2 * T_stat), chi_u / (2 * T_stat))
      list(ci = ci, pivot = "2*lambda*sum(X) ~ chi^2(2n)",
           estimand = "lambda (rate)", point_est = (n-1)/T_stat)
    },
    poisson = {
      # Exact CI based on Poisson / gamma duality
      T_stat <- sum(x)
      ci_l   <- qgamma(alpha/2, shape = T_stat, rate = n)
      ci_u   <- qgamma(1 - alpha/2, shape = T_stat + 1, rate = n)
      list(ci = c(ci_l, ci_u),
           pivot = "Exact CI via gamma-Poisson duality",
           estimand = "lambda", point_est = mean(x))
    },
    binomial = {
      # Wilson confidence interval
      T_stat <- sum(x)
      m      <- param$m %||% 1
      N      <- n * m
      p_hat  <- T_stat / N
      z      <- qnorm(1 - alpha/2)
      # Wilson interval
      denom  <- 1 + z^2 / N
      center <- (p_hat + z^2 / (2*N)) / denom
      margin <- z * sqrt(p_hat*(1-p_hat)/N + z^2/(4*N^2)) / denom
      ci     <- c(center - margin, center + margin)
      list(ci = ci, pivot = "Wilson interval for p",
           estimand = "p", point_est = p_hat)
    },
    gamma_rate = {
      alpha_shape <- param$alpha %||% param$shape %||% 1
      T_stat <- sum(x)
      chi_l  <- qchisq(alpha/2, df = 2*n*alpha_shape)
      chi_u  <- qchisq(1 - alpha/2, df = 2*n*alpha_shape)
      ci <- c(chi_l / (2*T_stat), chi_u / (2*T_stat))
      list(ci = ci,
           pivot = "2*beta*sum(X) ~ chi^2(2*n*alpha)",
           estimand = "beta (rate)", point_est = n*alpha_shape/T_stat)
    }
  )

  cat("=== Pivot-Based Confidence Interval ===\n")
  cat(sprintf("Family: %s,  n = %d,  conf_level = %.2f%%\n\n", family, n, 100*conf_level))
  if (!is.null(result$pivot))
    cat("Pivotal quantity:", result$pivot, "\n")
  if (!is.null(result$ci)) {
    cat(sprintf("Point estimate  : %.6f\n", result$point_est))
    cat(sprintf("%.0f%% CI: [%.6f, %.6f]\n", 100*conf_level, result$ci[1], result$ci[2]))
    cat(sprintf("CI width        : %.6f\n", diff(result$ci)))
  } else {
    cat("mu CI   :", round(result$ci_mu, 6), "\n")
    cat("sigma^2 CI:", round(result$ci_sigma2, 6), "\n")
  }

  invisible(result)
}

#' Shortest Expected Length Confidence Interval
#'
#' Finds the confidence interval of minimum expected length among all intervals
#' with the given coverage probability, by optimizing the allocation of alpha.
#'
#' @param x Numeric vector of observations
#' @param family Character: distribution family
#' @param conf_level Numeric: confidence level
#' @param param Additional fixed parameters
#'
#' @return Shortest CI with comparison to equal-tailed CI
#' @export
#'
#' @examples
#' x <- rexp(20, rate = 2)
#' shortest_ci(x, family = "exponential", conf_level = 0.95)
#'
#' x <- rgamma(20, shape = 3, rate = 2)
#' shortest_ci(x, family = "gamma_rate", conf_level = 0.95, param = list(alpha = 3))
shortest_ci <- function(x, family = c("exponential", "gamma_rate", "normal_var"),
                         conf_level = 0.95, param = list()) {
  family <- match.arg(family)
  n      <- length(x)
  alpha  <- 1 - conf_level

  result <- switch(family,
    exponential = {
      T_stat <- sum(x)
      df_val <- 2 * n
      # Length = (chi_upper - chi_lower) / (2*T)
      # Minimize chi_u - chi_l subject to F(chi_u) - F(chi_l) = conf_level
      # Optimal (shortest): f(chi_l) = f(chi_u) (equal density at endpoints)
      # For chi^2, numerically optimize:
      obj <- function(a_l) {
        chi_l <- qchisq(a_l, df = df_val)
        chi_u <- qchisq(a_l + conf_level, df = df_val)
        (chi_u - chi_l) / (2 * T_stat)
      }
      opt <- optimize(obj, interval = c(0, alpha), tol = 1e-8)
      a_l <- opt$minimum
      chi_l <- qchisq(a_l, df = df_val)
      chi_u <- qchisq(a_l + conf_level, df = df_val)
      ci_short <- c(chi_l / (2 * T_stat), chi_u / (2 * T_stat))

      # Equal-tailed for comparison
      ci_equal <- c(qchisq(alpha/2, df_val) / (2*T_stat),
                    qchisq(1-alpha/2, df_val) / (2*T_stat))
      list(ci_short = ci_short, ci_equal = ci_equal,
           length_short = diff(ci_short), length_equal = diff(ci_equal))
    },
    gamma_rate = {
      alpha_shape <- param$alpha %||% param$shape %||% 1
      T_stat <- sum(x)
      df_val <- 2 * n * alpha_shape
      obj <- function(a_l) {
        chi_l <- qchisq(a_l, df = df_val)
        chi_u <- qchisq(a_l + conf_level, df = df_val)
        (chi_u - chi_l) / (2 * T_stat)
      }
      opt    <- optimize(obj, interval = c(0, alpha), tol = 1e-8)
      chi_l  <- qchisq(opt$minimum, df = df_val)
      chi_u  <- qchisq(opt$minimum + conf_level, df = df_val)
      ci_short <- c(chi_l / (2*T_stat), chi_u / (2*T_stat))
      ci_equal <- c(qchisq(alpha/2, df_val) / (2*T_stat),
                    qchisq(1-alpha/2, df_val) / (2*T_stat))
      list(ci_short = ci_short, ci_equal = ci_equal,
           length_short = diff(ci_short), length_equal = diff(ci_equal))
    },
    normal_var = {
      s2     <- var(x)
      df_val <- n - 1
      obj <- function(a_l) {
        chi_l <- qchisq(a_l, df = df_val)
        chi_u <- qchisq(a_l + conf_level, df = df_val)
        (1/chi_l - 1/chi_u) * (n-1) * s2
      }
      opt    <- optimize(obj, interval = c(0, alpha), tol = 1e-8)
      chi_l  <- qchisq(opt$minimum, df = df_val)
      chi_u  <- qchisq(opt$minimum + conf_level, df = df_val)
      ci_short <- c((n-1)*s2/chi_u, (n-1)*s2/chi_l)
      ci_equal <- c((n-1)*s2 / qchisq(1-alpha/2, df_val),
                    (n-1)*s2 / qchisq(alpha/2, df_val))
      list(ci_short = ci_short, ci_equal = ci_equal,
           length_short = diff(ci_short), length_equal = diff(ci_equal))
    }
  )

  cat("=== Shortest Expected Length CI ===\n")
  cat(sprintf("Family: %s,  n = %d,  conf_level = %.0f%%\n\n", family, n, 100*conf_level))
  cat(sprintf("Shortest CI  : [%.6f, %.6f]  (length = %.6f)\n",
              result$ci_short[1], result$ci_short[2], result$length_short))
  cat(sprintf("Equal-tailed : [%.6f, %.6f]  (length = %.6f)\n",
              result$ci_equal[1], result$ci_equal[2], result$length_equal))
  cat(sprintf("Length reduction: %.2f%%\n",
              100 * (result$length_equal - result$length_short) / result$length_equal))

  invisible(result)
}

#' Uniformly Most Accurate (UMA) Confidence Interval
#'
#' Constructs the UMA confidence interval obtained by inverting the UMP test.
#' For one-sided H0: theta <= theta0, the UMA CI gives the smallest upper
#' confidence bound (or largest lower bound) uniformly over all theta.
#'
#' @param x Numeric vector of observations
#' @param family Character: distribution family
#' @param conf_level Numeric: confidence level
#' @param side Character: "upper" (upper confidence bound) or "lower"
#' @param param Additional fixed parameters
#'
#' @return UMA confidence bound
#' @export
#'
#' @examples
#' x <- rpois(20, lambda = 3)
#' uma_ci(x, family = "poisson", conf_level = 0.95, side = "upper")
#'
#' x <- rexp(15, rate = 2)
#' uma_ci(x, family = "exponential", conf_level = 0.95, side = "lower")
uma_ci <- function(x, family = c("normal_mean", "poisson", "exponential", "binomial"),
                    conf_level = 0.95, side = c("upper", "lower"), param = list()) {
  family <- match.arg(family)
  side   <- match.arg(side)
  n      <- length(x)
  alpha  <- 1 - conf_level

  T_obs <- sum(x)

  # UMA CI obtained by inverting UMP test
  # For UMP test of H0: theta <= theta0 at level alpha (rejecting for large T):
  # Upper 1-alpha confidence bound: smallest theta0 such that T_obs is not in rejection region
  bound <- switch(family,
    normal_mean = {
      sigma <- param$sigma %||% 1
      z     <- qnorm(conf_level)
      if (side == "upper")
        T_obs/n + z * sigma/sqrt(n)
      else
        T_obs/n - z * sigma/sqrt(n)
    },
    poisson = {
      # Invert: reject theta0 when sum(X) >= c(theta0)
      # Upper bound: smallest theta0 s.t. P_{theta0}(T >= T_obs) >= alpha
      if (side == "upper")
        qgamma(conf_level, shape = T_obs + 1, rate = n)
      else
        qgamma(1 - conf_level, shape = T_obs, rate = n)
    },
    exponential = {
      if (side == "lower")
        qchisq(1 - conf_level, df = 2*(T_obs + 1)) / (2 * T_obs/n)
      else
        2 * T_obs / qchisq(1 - conf_level, df = 2*n)
    },
    binomial = {
      m <- param$m %||% 1
      N <- n * m
      T <- sum(x)
      if (side == "upper")
        qbeta(conf_level, T + 1, N - T)
      else
        qbeta(1 - conf_level, T, N - T + 1)
    }
  )

  cat("=== Uniformly Most Accurate (UMA) CI ===\n")
  cat(sprintf("Family: %s,  n = %d,  conf_level = %.0f%%\n", family, n, 100*conf_level))
  cat(sprintf("T(x) = sum(x) = %.4f\n\n", T_obs))
  cat(sprintf("UMA %s confidence bound: %.6f\n",
              ifelse(side=="upper","upper","lower"), bound))
  cat(sprintf("This gives the %.0f%% one-sided CI: %s\n", 100*conf_level,
              ifelse(side=="upper", paste0("(-Inf, ", round(bound,4), "]"),
                     paste0("[", round(bound,4), ", Inf)"))))
  cat("\nTheory: UMA CI inverts the UMP one-sided test.\n")
  cat("It minimizes P(theta in CI | theta = theta') for all theta' != true theta.\n")

  invisible(list(bound = bound, side = side, conf_level = conf_level))
}

#' Uniformly Most Accurate Unbiased (UMAU) Confidence Interval
#'
#' Constructs the UMAU confidence interval by inverting the UMPU test.
#' Two-sided CI that is unbiased (coverage >= conf_level) and UMA among unbiased CIs.
#'
#' @param x Numeric vector of observations
#' @param family Character: distribution family
#' @param conf_level Numeric: confidence level
#' @param param Additional fixed parameters
#'
#' @return UMAU two-sided confidence interval
#' @export
#'
#' @examples
#' x <- rpois(20, lambda = 3)
#' umau_ci(x, family = "poisson", conf_level = 0.95)
#'
#' x <- rnorm(25, mean = 2, sd = 1)
#' umau_ci(x, family = "normal_mean", conf_level = 0.95, param = list(sigma = 1))
umau_ci <- function(x, family = c("normal_mean", "poisson", "exponential", "binomial"),
                     conf_level = 0.95, param = list()) {
  family <- match.arg(family)
  n      <- length(x)
  alpha  <- 1 - conf_level

  T_obs  <- sum(x)
  half   <- conf_level + alpha/2  # = 1 - alpha/2

  # Invert UMPU two-sided test
  ci <- switch(family,
    normal_mean = {
      sigma <- param$sigma %||% 1
      z <- qnorm(1 - alpha/2)
      c(T_obs/n - z * sigma/sqrt(n), T_obs/n + z * sigma/sqrt(n))
    },
    poisson = {
      ci_l <- qgamma(alpha/2,   shape = T_obs,     rate = n)
      ci_u <- qgamma(1-alpha/2, shape = T_obs + 1, rate = n)
      c(ci_l, ci_u)
    },
    exponential = {
      c(qchisq(alpha/2,   df = 2*n) / (2*T_obs),
        qchisq(1-alpha/2, df = 2*n) / (2*T_obs))
    },
    binomial = {
      m <- param$m %||% 1
      N <- n * m
      T <- sum(x)
      c(qbeta(alpha/2,   T,   N - T + 1),
        qbeta(1-alpha/2, T+1, N - T))
    }
  )

  cat("=== Uniformly Most Accurate Unbiased (UMAU) CI ===\n")
  cat(sprintf("Family: %s,  n = %d,  conf_level = %.0f%%\n", family, n, 100*conf_level))
  cat(sprintf("T(x) = sum(x) = %.4f\n\n", T_obs))
  cat(sprintf("UMAU %.0f%% CI: [%.6f, %.6f]\n", 100*conf_level, ci[1], ci[2]))
  cat(sprintf("CI width     : %.6f\n", diff(ci)))
  cat("\nTheory: UMAU CI inverts the UMPU two-sided test.\n")
  cat("It is unbiased (coverage >= nominal) and UMA among all unbiased CIs.\n")

  invisible(list(ci = ci, conf_level = conf_level, T_obs = T_obs))
}
