#' Neyman-Pearson Lemma: Most Powerful Test
#'
#' Constructs the most powerful (MP) test for a simple H0 vs simple H1,
#' using the likelihood ratio criterion. Rejects H0 when L(theta1)/L(theta0) > k.
#'
#' @param x Numeric vector of observations
#' @param family Character: distribution family
#' @param theta0 Numeric: null hypothesis parameter value
#' @param theta1 Numeric: alternative hypothesis parameter value
#' @param alpha Numeric: significance level (0 < alpha < 1)
#' @param param Additional fixed parameters
#'
#' @return List with test statistic, critical value, rejection decision, and power
#' @export
#'
#' @examples
#' set.seed(1); x <- rnorm(20, mean = 1.5, sd = 1)
#' np_lemma(x, "normal_mean", theta0 = 0, theta1 = 2, alpha = 0.05,
#'          param = list(sigma = 1))
#'
#' x <- rpois(15, lambda = 4)
#' np_lemma(x, "poisson", theta0 = 2, theta1 = 5, alpha = 0.05)
np_lemma <- function(x, family = c("normal_mean", "poisson", "exponential",
                                    "binomial", "bernoulli"),
                      theta0, theta1, alpha = 0.05, param = list()) {
  family <- match.arg(family)
  n <- length(x)

  # Log-likelihood ratio statistic
  log_lr <- switch(family,
    normal_mean = {
      sigma <- param$sigma %||% 1
      sum(dnorm(x, mean = theta1, sd = sigma, log = TRUE)) -
        sum(dnorm(x, mean = theta0, sd = sigma, log = TRUE))
    },
    poisson = {
      sum(dpois(round(x), lambda = theta1, log = TRUE)) -
        sum(dpois(round(x), lambda = theta0, log = TRUE))
    },
    exponential = {
      sum(dexp(x, rate = theta1, log = TRUE)) -
        sum(dexp(x, rate = theta0, log = TRUE))
    },
    binomial = {
      m <- param$m %||% param$n %||% 1
      sum(dbinom(round(x), size = m, prob = theta1, log = TRUE)) -
        sum(dbinom(round(x), size = m, prob = theta0, log = TRUE))
    },
    bernoulli = {
      sum(dbinom(round(x), size = 1, prob = theta1, log = TRUE)) -
        sum(dbinom(round(x), size = 1, prob = theta0, log = TRUE))
    }
  )

  # Sufficient statistic and its critical value
  T_obs <- switch(family,
    normal_mean = mean(x),
    poisson     = sum(x),
    exponential = sum(x),
    binomial    = sum(x),
    bernoulli   = sum(x)
  )

  # Critical value and power via exact/normal distribution
  crit_info <- .mp_critical_value(family, theta0, theta1, n, alpha, param)
  reject     <- T_obs > crit_info$c_alpha | T_obs < crit_info$c_lower

  cat("=== Neyman-Pearson Most Powerful Test ===\n")
  cat(sprintf("H0: theta = %.4f  vs  H1: theta = %.4f\n", theta0, theta1))
  cat(sprintf("alpha = %.4f,  n = %d,  family = %s\n\n", alpha, n, family))
  cat(sprintf("Sufficient statistic T(x) = %.4f\n", T_obs))
  cat(sprintf("Log-likelihood ratio      = %.4f\n", log_lr))
  cat(sprintf("Critical value c(alpha)   = %.4f\n", crit_info$c_alpha))
  cat(sprintf("Power at theta1           = %.4f\n", crit_info$power))
  cat(sprintf("Decision: %s H0\n", ifelse(reject, "REJECT", "FAIL TO REJECT")))
  cat("\nN-P Lemma: This is the MOST POWERFUL test at level alpha.\n")
  cat("The MP test rejects when the likelihood ratio L(theta1)/L(theta0) > k_alpha.\n")

  invisible(list(
    T_obs    = T_obs,
    log_lr   = log_lr,
    c_alpha  = crit_info$c_alpha,
    power    = crit_info$power,
    reject   = reject,
    alpha    = alpha
  ))
}

# Internal helper: critical values for common MP tests
.mp_critical_value <- function(family, theta0, theta1, n, alpha, param) {
  sigma <- param$sigma %||% 1
  m     <- param$m %||% param$n %||% 1
  direction <- if (theta1 > theta0) "upper" else "lower"

  result <- switch(family,
    normal_mean = {
      # T = X_bar ~ N(theta0, sigma^2/n) under H0
      se <- sigma / sqrt(n)
      c_alpha <- if (direction == "upper")
        qnorm(1 - alpha, mean = theta0 * n, sd = sigma * sqrt(n))
      else
        qnorm(alpha, mean = theta0 * n, sd = sigma * sqrt(n))
      power <- if (direction == "upper")
        pnorm(c_alpha, mean = theta1 * n, sd = sigma * sqrt(n), lower.tail = FALSE)
      else
        pnorm(c_alpha, mean = theta1 * n, sd = sigma * sqrt(n))
      list(c_alpha = c_alpha / n, power = power, c_lower = -Inf)
    },
    poisson = {
      # T = sum(X_i) ~ Poisson(n*theta) under H0
      c_alpha <- if (direction == "upper")
        qpois(1 - alpha, lambda = n * theta0)
      else
        qpois(alpha, lambda = n * theta0)
      power <- if (direction == "upper")
        ppois(c_alpha, lambda = n * theta1, lower.tail = FALSE)
      else
        ppois(c_alpha, lambda = n * theta1)
      list(c_alpha = c_alpha, power = power, c_lower = -Inf)
    },
    exponential = {
      # T = sum(X_i) ~ Gamma(n, theta0) under H0
      c_alpha <- if (direction == "lower")
        qgamma(alpha, shape = n, rate = theta0)
      else
        qgamma(1 - alpha, shape = n, rate = theta0)
      power <- if (direction == "lower")
        pgamma(c_alpha, shape = n, rate = theta1)
      else
        pgamma(c_alpha, shape = n, rate = theta1, lower.tail = FALSE)
      list(c_alpha = c_alpha, power = power, c_lower = -Inf)
    },
    binomial = , bernoulli = {
      total_trials <- n * m
      c_alpha <- if (direction == "upper")
        qbinom(1 - alpha, size = total_trials, prob = theta0)
      else
        qbinom(alpha, size = total_trials, prob = theta0)
      power <- if (direction == "upper")
        pbinom(c_alpha, size = total_trials, prob = theta1, lower.tail = FALSE)
      else
        pbinom(c_alpha, size = total_trials, prob = theta1)
      list(c_alpha = c_alpha, power = power, c_lower = -Inf)
    }
  )
  result
}

#' Uniformly Most Powerful (UMP) Test
#'
#' Constructs a UMP test for one-sided hypotheses in exponential families.
#' For H0: theta <= theta0 vs H1: theta > theta0, the UMP test rejects
#' when T(X) > c_alpha (one-sided test based on sufficient statistic).
#'
#' @param x Numeric vector of observations
#' @param family Character: distribution family
#' @param theta0 Numeric: boundary value under H0
#' @param alternative Character: "greater" or "less"
#' @param alpha Numeric: significance level
#' @param theta_seq Numeric vector: theta values for power curve
#' @param param Additional fixed parameters
#'
#' @return List with test results and power function
#' @export
#'
#' @examples
#' x <- rpois(20, lambda = 3.5)
#' ump_test(x, "poisson", theta0 = 2, alternative = "greater", alpha = 0.05,
#'          theta_seq = seq(0.5, 6, by = 0.25))
ump_test <- function(x, family = c("normal_mean", "poisson", "exponential",
                                    "binomial", "bernoulli"),
                      theta0, alternative = c("greater", "less"),
                      alpha = 0.05, theta_seq = NULL, param = list()) {
  family      <- match.arg(family)
  alternative <- match.arg(alternative)
  n <- length(x)

  T_obs <- switch(family,
    normal_mean = sum(x),
    poisson     = sum(x),
    exponential = sum(x),
    binomial    = sum(x),
    bernoulli   = sum(x)
  )

  # Critical value and test
  crit <- .mp_critical_value(family, theta0, theta0, n, alpha, param)
  c_alpha <- switch(family,
    normal_mean = {
      sigma <- param$sigma %||% 1
      if (alternative == "greater")
        qnorm(1 - alpha, mean = theta0 * n, sd = sigma * sqrt(n))
      else
        qnorm(alpha, mean = theta0 * n, sd = sigma * sqrt(n))
    },
    poisson = {
      if (alternative == "greater")
        qpois(1 - alpha, lambda = n * theta0)
      else
        qpois(alpha, lambda = n * theta0)
    },
    exponential = {
      if (alternative == "greater")
        qgamma(alpha, shape = n, rate = theta0)  # large theta = small sum
      else
        qgamma(1 - alpha, shape = n, rate = theta0)
    },
    binomial = , bernoulli = {
      m <- param$m %||% param$n %||% 1
      total_n <- n * m
      if (alternative == "greater")
        qbinom(1 - alpha, size = total_n, prob = theta0)
      else
        qbinom(alpha, size = total_n, prob = theta0)
    }
  )

  reject <- if (alternative == "greater") T_obs > c_alpha else T_obs < c_alpha

  # Power function
  power_fn <- NULL
  if (!is.null(theta_seq)) {
    power_fn <- sapply(theta_seq, function(th) {
      .power_at_theta(family, th, n, c_alpha, alternative, param)
    })
  }

  cat("=== Uniformly Most Powerful (UMP) Test ===\n")
  cat(sprintf("H0: theta %s %.4f  vs  H1: theta %s %.4f\n",
              ifelse(alternative=="greater","<=",">="), theta0,
              ifelse(alternative=="greater",">","<"), theta0))
  cat(sprintf("alpha = %.4f,  n = %d,  family = %s\n\n", alpha, n, family))
  cat(sprintf("T(x) = sum(x) = %.4f\n", T_obs))
  cat(sprintf("Critical value c(alpha) = %.4f\n", c_alpha))
  cat(sprintf("Decision: %s H0\n", ifelse(reject, "REJECT", "FAIL TO REJECT")))
  cat("\nTheory: For one-sided hypotheses in exponential families with\n")
  cat("monotone likelihood ratio, the UMP test exists and is given by\n")
  cat("the one-sided test on the sufficient statistic.\n")

  out <- list(T_obs = T_obs, c_alpha = c_alpha, reject = reject,
              alpha = alpha, family = family, theta0 = theta0,
              alternative = alternative, theta_seq = theta_seq,
              power = power_fn)
  class(out) <- "power_curve"
  invisible(out)
}

#' Power Function Computation
#'
#' Computes and optionally plots the power function of a test at given theta values.
#'
#' @param x Numeric vector of observations
#' @param family Character: distribution family
#' @param theta0 Numeric: null hypothesis value
#' @param theta_seq Numeric vector: theta values for power curve
#' @param alternative Character: "greater", "less", or "two.sided"
#' @param alpha Numeric: significance level
#' @param param Additional fixed parameters
#' @param plot Logical: whether to plot the power curve
#'
#' @return Data frame of theta values and corresponding power
#' @export
#'
#' @examples
#' x <- rnorm(25, mean = 1, sd = 1)
#' power_function(x, "normal_mean", theta0 = 0,
#'                theta_seq = seq(-2, 3, by = 0.1), alpha = 0.05,
#'                param = list(sigma = 1), plot = TRUE)
power_function <- function(x, family, theta0,
                            theta_seq = seq(theta0 - 3, theta0 + 3, length.out = 100),
                            alternative = c("greater", "less", "two.sided"),
                            alpha = 0.05, param = list(), plot = TRUE) {
  alternative <- match.arg(alternative)
  n  <- length(x)
  sigma <- param$sigma %||% 1
  m     <- param$m %||% param$n %||% 1

  # Compute critical region
  c_upper <- switch(family,
    normal_mean = qnorm(1 - alpha/ifelse(alternative=="two.sided",2,1),
                        mean = theta0 * n, sd = sigma * sqrt(n)),
    poisson     = qpois(1 - alpha/ifelse(alternative=="two.sided",2,1),
                        lambda = n * theta0),
    exponential = qgamma(alpha/ifelse(alternative=="two.sided",2,1),
                         shape = n, rate = theta0),
    binomial    = qbinom(1 - alpha/ifelse(alternative=="two.sided",2,1),
                         size = n*m, prob = theta0),
    bernoulli   = qbinom(1 - alpha/ifelse(alternative=="two.sided",2,1),
                         size = n, prob = theta0)
  )
  c_lower <- switch(family,
    normal_mean = qnorm(alpha/ifelse(alternative=="two.sided",2,1),
                        mean = theta0 * n, sd = sigma * sqrt(n)),
    poisson     = qpois(alpha/ifelse(alternative=="two.sided",2,1),
                        lambda = n * theta0),
    exponential = qgamma(1 - alpha/ifelse(alternative=="two.sided",2,1),
                         shape = n, rate = theta0),
    binomial    = qbinom(alpha/ifelse(alternative=="two.sided",2,1),
                         size = n*m, prob = theta0),
    bernoulli   = qbinom(alpha/ifelse(alternative=="two.sided",2,1),
                         size = n, prob = theta0)
  )

  powers <- sapply(theta_seq, function(th) {
    p_upper <- .power_at_theta(family, th, n, c_upper, "greater", param)
    p_lower <- .power_at_theta(family, th, n, c_lower, "less",    param)
    if (alternative == "greater") p_upper
    else if (alternative == "less") p_lower
    else p_upper + p_lower
  })
  powers <- pmin(1, pmax(0, powers))

  df <- data.frame(theta = theta_seq, power = powers)

  if (plot) {
    plot(theta_seq, powers, type = "l", col = "steelblue", lwd = 2,
         xlab = expression(theta), ylab = "Power",
         main = sprintf("Power Function: %s, n=%d, alpha=%.2f", family, n, alpha),
         ylim = c(0, 1))
    abline(h = alpha, lty = 2, col = "red")
    abline(v = theta0, lty = 3, col = "gray40")
    legend("bottomright", c("Power", paste("alpha =", alpha), "theta0"),
           lty = c(1,2,3), col = c("steelblue","red","gray40"), bty = "n")
  }

  invisible(df)
}

# Internal: power at a given theta
.power_at_theta <- function(family, theta, n, c_alpha, direction, param) {
  sigma <- param$sigma %||% 1
  m     <- param$m %||% param$n %||% 1
  switch(family,
    normal_mean = {
      if (direction == "greater")
        pnorm(c_alpha, mean = theta * n, sd = sigma * sqrt(n), lower.tail = FALSE)
      else
        pnorm(c_alpha, mean = theta * n, sd = sigma * sqrt(n))
    },
    poisson = {
      if (direction == "greater")
        ppois(c_alpha, lambda = n * theta, lower.tail = FALSE)
      else
        ppois(c_alpha, lambda = n * theta)
    },
    exponential = {
      if (direction == "greater")
        pgamma(c_alpha, shape = n, rate = theta, lower.tail = FALSE)
      else
        pgamma(c_alpha, shape = n, rate = theta)
    },
    binomial = , bernoulli = {
      total_n <- n * m
      if (direction == "greater")
        pbinom(c_alpha, size = total_n, prob = theta, lower.tail = FALSE)
      else
        pbinom(c_alpha, size = total_n, prob = theta)
    }
  )
}

#' Uniformly Most Powerful Unbiased (UMPU) Test
#'
#' Constructs the UMPU test for two-sided hypotheses H0: theta = theta0 vs
#' H1: theta ≠ theta0, or for hypotheses in two-parameter exponential families
#' using the Neyman structure.
#'
#' @param x Numeric vector of observations
#' @param family Character: distribution family
#' @param theta0 Numeric: null hypothesis value
#' @param alpha Numeric: significance level
#' @param param Additional fixed parameters
#'
#' @return Test result with UMPU critical region and p-value
#' @export
#'
#' @examples
#' x <- rpois(15, lambda = 2.8)
#' umpu_test(x, "poisson", theta0 = 2, alpha = 0.05)
#'
#' x <- rnorm(20, mean = 1.2, sd = 1)
#' umpu_test(x, "normal_mean", theta0 = 0, alpha = 0.05, param = list(sigma = 1))
umpu_test <- function(x, family = c("normal_mean", "poisson", "exponential",
                                     "binomial"),
                       theta0, alpha = 0.05, param = list()) {
  family <- match.arg(family)
  n <- length(x)

  T_obs <- sum(x)

  # Two-sided critical values
  c1 <- switch(family,
    normal_mean = {
      sigma <- param$sigma %||% 1
      qnorm(alpha/2, mean = n * theta0, sd = sigma * sqrt(n))
    },
    poisson  = qpois(alpha/2, lambda = n * theta0),
    exponential = qgamma(alpha/2, shape = n, rate = theta0),
    binomial = qbinom(alpha/2, size = n * (param$m %||% 1), prob = theta0)
  )
  c2 <- switch(family,
    normal_mean = {
      sigma <- param$sigma %||% 1
      qnorm(1 - alpha/2, mean = n * theta0, sd = sigma * sqrt(n))
    },
    poisson  = qpois(1 - alpha/2, lambda = n * theta0),
    exponential = qgamma(1 - alpha/2, shape = n, rate = theta0),
    binomial = qbinom(1 - alpha/2, size = n * (param$m %||% 1), prob = theta0)
  )

  reject <- T_obs < c1 | T_obs > c2

  # Approximate p-value
  p_val <- switch(family,
    normal_mean = {
      sigma <- param$sigma %||% 1
      z <- (T_obs - n * theta0) / (sigma * sqrt(n))
      2 * pnorm(-abs(z))
    },
    poisson     = 2 * min(ppois(T_obs, n * theta0),
                          ppois(T_obs - 1, n * theta0, lower.tail = FALSE)),
    exponential = 2 * min(pgamma(T_obs, n, theta0),
                          pgamma(T_obs, n, theta0, lower.tail = FALSE)),
    binomial    = {
      total_n <- n * (param$m %||% 1)
      2 * min(pbinom(T_obs, total_n, theta0),
              pbinom(T_obs - 1, total_n, theta0, lower.tail = FALSE))
    }
  )

  cat("=== UMPU Test (Two-Sided) ===\n")
  cat(sprintf("H0: theta = %.4f  vs  H1: theta != %.4f\n", theta0, theta0))
  cat(sprintf("alpha = %.4f,  n = %d,  family = %s\n\n", alpha, n, family))
  cat(sprintf("T(x) = sum(x) = %.4f\n", T_obs))
  cat(sprintf("Critical region: T < %.4f  or  T > %.4f\n", c1, c2))
  cat(sprintf("p-value approx.: %.6f\n", p_val))
  cat(sprintf("Decision: %s H0\n", ifelse(reject, "REJECT", "FAIL TO REJECT")))
  cat("\nTheory: UMPU test has unbiased power function (power >= alpha for all\n")
  cat("theta != theta0) and is optimal within unbiased tests.\n")

  invisible(list(T_obs = T_obs, c1 = c1, c2 = c2, p_value = p_val,
                 reject = reject, alpha = alpha))
}

#' Monotone Likelihood Ratio Test
#'
#' Tests hypotheses using the MLR property. If the family has MLR in T(X),
#' then UMP tests exist for one-sided hypotheses based on T.
#'
#' @param x Numeric vector of observations
#' @param family Character: distribution family
#' @param theta0 Numeric: null hypothesis value
#' @param theta1 Numeric: alternative hypothesis value (for direction)
#' @param alpha Numeric: significance level
#' @param param Additional fixed parameters
#'
#' @return Test result with MLR verification
#' @export
mlr_test <- function(x, family = c("normal_mean", "poisson", "exponential",
                                    "binomial"),
                      theta0, theta1, alpha = 0.05, param = list()) {
  family <- match.arg(family)
  direction <- if (theta1 > theta0) "greater" else "less"

  cat("=== Monotone Likelihood Ratio Test ===\n")
  cat(sprintf("Family: %s\n", family))
  cat("MLR property: The likelihood ratio f(x;theta2)/f(x;theta1) is\n")
  cat("non-decreasing in T(x) for theta2 > theta1.\n\n")

  # MLR statistic descriptions
  mlr_info <- switch(family,
    normal_mean = "T(X) = X_bar (or sum Xi); MLR in X_bar for normal location",
    poisson     = "T(X) = sum(Xi); LR = exp((theta2-theta1)*T) * C(theta1,theta2)\n  => monotone in sum(Xi)",
    exponential = "T(X) = sum(Xi); LR = (theta2/theta1)^n * exp(-(theta2-theta1)*T)\n  => monotone (decreasing) in sum(Xi)",
    binomial    = "T(X) = sum(Xi); MLR in T for binomial with varying p"
  )
  cat("MLR statistic:", mlr_info, "\n\n")

  result <- ump_test(x, family, theta0, alternative = direction,
                      alpha = alpha, param = param)
  invisible(result)
}

#' Likelihood Ratio Test
#'
#' Computes the generalized likelihood ratio test statistic:
#' Lambda = sup_{theta in H0} L(theta) / sup_{theta in Omega} L(theta)
#' Rejects when -2 log(Lambda) > chi-squared critical value.
#'
#' @param x Numeric vector of observations
#' @param family Character: distribution family
#' @param theta0 Null hypothesis value (or vector for composite H0)
#' @param alpha Numeric: significance level
#' @param param Additional fixed parameters
#'
#' @return LRT statistic, p-value, and decision
#' @export
#'
#' @examples
#' x <- rnorm(25, mean = 1.5, sd = 1)
#' lrt(x, "normal_mean", theta0 = 0, param = list(sigma = 1))
#'
#' x <- rpois(20, lambda = 3)
#' lrt(x, "poisson", theta0 = 2)
lrt <- function(x, family = c("normal_mean", "poisson", "exponential", "binomial"),
                 theta0, alpha = 0.05, param = list()) {
  family <- match.arg(family)
  n <- length(x)

  # MLE (unrestricted)
  theta_mle <- switch(family,
    normal_mean = mean(x),
    poisson     = mean(x),
    exponential = 1 / mean(x),
    binomial    = {m <- param$m %||% 1; mean(x) / m}
  )

  # Log-likelihoods
  logL_H0 <- switch(family,
    normal_mean = {
      sigma <- param$sigma %||% 1
      sum(dnorm(x, mean = theta0, sd = sigma, log = TRUE))
    },
    poisson     = sum(dpois(round(x), lambda = theta0, log = TRUE)),
    exponential = sum(dexp(x, rate = theta0, log = TRUE)),
    binomial    = {
      m <- param$m %||% 1
      sum(dbinom(round(x), size = m, prob = theta0, log = TRUE))
    }
  )

  logL_MLE <- switch(family,
    normal_mean = {
      sigma <- param$sigma %||% 1
      sum(dnorm(x, mean = theta_mle, sd = sigma, log = TRUE))
    },
    poisson     = sum(dpois(round(x), lambda = theta_mle, log = TRUE)),
    exponential = sum(dexp(x, rate = theta_mle, log = TRUE)),
    binomial    = {
      m <- param$m %||% 1
      sum(dbinom(round(x), size = m, prob = theta_mle, log = TRUE))
    }
  )

  lambda_stat  <- -2 * (logL_H0 - logL_MLE)
  df           <- 1  # testing one parameter
  p_val        <- pchisq(lambda_stat, df = df, lower.tail = FALSE)
  reject       <- lambda_stat > qchisq(1 - alpha, df = df)

  cat("=== Likelihood Ratio Test ===\n")
  cat(sprintf("H0: theta = %.4f  vs  H1: theta != %.4f\n", theta0, theta0))
  cat(sprintf("Family: %s,  n = %d,  alpha = %.4f\n\n", family, n, alpha))
  cat(sprintf("MLE theta_hat       = %.6f\n", theta_mle))
  cat(sprintf("log L(theta0)       = %.6f\n", logL_H0))
  cat(sprintf("log L(theta_hat)    = %.6f\n", logL_MLE))
  cat(sprintf("\n-2 log Lambda       = %.6f\n", lambda_stat))
  cat(sprintf("Chi-sq critical val = %.6f  (df = %d)\n",
              qchisq(1 - alpha, df), df))
  cat(sprintf("p-value             = %.6f\n", p_val))
  cat(sprintf("Decision            : %s H0\n", ifelse(reject, "REJECT", "FAIL TO REJECT")))
  cat("\nAsymptotic theory: -2 log Lambda -> chi^2(df) under H0 (Wilks' theorem).\n")

  invisible(list(
    lambda_stat = lambda_stat,
    theta_mle   = theta_mle,
    logL_H0     = logL_H0,
    logL_mle    = logL_MLE,
    p_value     = p_val,
    reject      = reject,
    df          = df
  ))
}

#' @export
plot.power_curve <- function(x, ...) {
  if (!is.null(x$theta_seq) && !is.null(x$power)) {
    plot(x$theta_seq, x$power, type = "l", col = "steelblue", lwd = 2,
         xlab = expression(theta), ylab = "Power",
         main = sprintf("Power Curve: %s, n = %d", x$family, length(x$T_obs)),
         ylim = c(0, 1), ...)
    abline(h = x$alpha, lty = 2, col = "red")
    abline(v = x$theta0, lty = 3, col = "gray50")
    legend("bottomright", c("Power", paste("alpha =", x$alpha)),
           lty = c(1, 2), col = c("steelblue", "red"), bty = "n")
  }
  invisible(x)
}

#' Neyman-Pearson Critical Region
#'
#' Visualizes and computes the rejection region for a NP test.
#'
#' @param family Character: distribution family
#' @param theta0 Numeric: null parameter value
#' @param theta1 Numeric: alternative parameter value
#' @param n Integer: sample size
#' @param alpha Numeric: significance level
#' @param param Additional fixed parameters
#'
#' @return List with critical region and visualization
#' @export
neyman_pearson_region <- function(family, theta0, theta1, n = 20, alpha = 0.05, param = list()) {
  direction <- if (theta1 > theta0) "greater" else "less"
  crit <- .mp_critical_value(family, theta0, theta1, n, alpha, param)

  cat("=== Neyman-Pearson Critical Region ===\n")
  cat(sprintf("H0: theta = %.4f  vs  H1: theta = %.4f\n", theta0, theta1))
  cat(sprintf("Family: %s,  n = %d\n\n", family, n))
  cat(sprintf("Reject H0 when T(x) %s %.4f\n",
              ifelse(direction=="greater", ">", "<"), crit$c_alpha))
  cat(sprintf("Size (actual alpha) ≈ %.4f\n", alpha))
  cat(sprintf("Power at theta1 = %.4f : %.4f\n", theta1, crit$power))

  invisible(crit)
}
