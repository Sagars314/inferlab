#' Exponential Family Representation
#'
#' Creates an object representing a one-parameter or multi-parameter exponential family
#' in canonical form: f(x; theta) = h(x) * exp(eta(theta) * T(x) - A(theta))
#'
#' @param family Character string: "normal", "binomial", "poisson", "exponential",
#'   "gamma", "bernoulli", or "custom"
#' @param param Named list of fixed parameters (e.g., list(sigma=1) for normal)
#' @param custom_fns For custom families: list with h, eta, T_stat, A functions
#'
#' @return An object of class "expfam" containing canonical components
#' @export
#'
#' @examples
#' ef_norm <- exponential_family("normal", param = list(sigma = 1))
#' ef_pois <- exponential_family("poisson")
#' print(ef_norm)
exponential_family <- function(family = c("normal", "binomial", "poisson",
                                           "exponential", "gamma", "bernoulli",
                                           "custom"),
                                param = list(),
                                custom_fns = NULL) {
  family <- match.arg(family)

  ef <- switch(family,
    normal = {
      sigma <- param$sigma %||% 1
      list(
        family      = "normal",
        description = sprintf("N(mu, sigma^2=%.2f), canonical param eta=mu/sigma^2", sigma^2),
        h           = function(x) exp(-x^2 / (2 * sigma^2)) / (sqrt(2 * pi) * sigma),
        eta         = function(theta) theta / sigma^2,
        T_stat      = function(x) x,
        A           = function(theta) theta^2 / (2 * sigma^2),
        nat_param   = "eta = mu / sigma^2",
        suff_stat   = "T(X) = X (sample mean for iid sample)",
        support     = "(-Inf, Inf)",
        param_fixed = param
      )
    },
    binomial = {
      n <- param$n %||% 1
      list(
        family      = "binomial",
        description = sprintf("Bin(%d, p), canonical param eta=log(p/(1-p))", n),
        h           = function(x) choose(n, x),
        eta         = function(theta) log(theta / (1 - theta)),
        T_stat      = function(x) x,
        A           = function(theta) n * log(1 + exp(log(theta / (1 - theta)))),
        nat_param   = "eta = log(p/(1-p))  [log-odds]",
        suff_stat   = "T(X) = sum(X_i)  [total successes]",
        support     = paste0("{0,1,...,", n, "}"),
        param_fixed = param
      )
    },
    poisson = {
      list(
        family      = "poisson",
        description = "Poisson(lambda), canonical param eta=log(lambda)",
        h           = function(x) 1 / factorial(x),
        eta         = function(theta) log(theta),
        T_stat      = function(x) x,
        A           = function(theta) theta,
        nat_param   = "eta = log(lambda)",
        suff_stat   = "T(X) = sum(X_i)  [total count]",
        support     = "{0, 1, 2, ...}",
        param_fixed = param
      )
    },
    exponential = {
      list(
        family      = "exponential",
        description = "Exp(lambda), canonical param eta=-lambda",
        h           = function(x) ifelse(x >= 0, 1, 0),
        eta         = function(theta) -theta,
        T_stat      = function(x) x,
        A           = function(theta) -log(theta),
        nat_param   = "eta = -lambda",
        suff_stat   = "T(X) = sum(X_i)  [total time]",
        support     = "[0, Inf)",
        param_fixed = param
      )
    },
    gamma = {
      alpha <- param$alpha %||% param$shape %||% 1
      list(
        family      = "gamma",
        description = sprintf("Gamma(alpha=%.2f, beta), canonical param eta=-beta", alpha),
        h           = function(x) ifelse(x > 0, x^(alpha - 1) / gamma(alpha), 0),
        eta         = function(theta) -theta,
        T_stat      = function(x) x,
        A           = function(theta) -alpha * log(theta),
        nat_param   = "eta = -beta",
        suff_stat   = "T(X) = sum(X_i)",
        support     = "(0, Inf)",
        param_fixed = param
      )
    },
    bernoulli = {
      list(
        family      = "bernoulli",
        description = "Bernoulli(p), canonical param eta=log(p/(1-p))",
        h           = function(x) 1,
        eta         = function(theta) log(theta / (1 - theta)),
        T_stat      = function(x) x,
        A           = function(theta) log(1 + exp(log(theta / (1 - theta)))),
        nat_param   = "eta = log(p/(1-p))  [log-odds]",
        suff_stat   = "T(X) = sum(X_i)  [number of successes]",
        support     = "{0, 1}",
        param_fixed = param
      )
    },
    custom = {
      if (is.null(custom_fns)) stop("Provide 'custom_fns' list for custom family.")
      required <- c("h", "eta", "T_stat", "A")
      missing_fns <- setdiff(required, names(custom_fns))
      if (length(missing_fns) > 0)
        stop("custom_fns must contain: ", paste(missing_fns, collapse = ", "))
      c(list(family = "custom",
             description = custom_fns$description %||% "Custom exponential family"),
        custom_fns)
    }
  )

  class(ef) <- "expfam"
  ef
}

#' @export
print.expfam <- function(x, ...) {
  cat("=== Exponential Family ===\n")
  cat("Family     :", x$family, "\n")
  cat("Description:", x$description, "\n")
  if (!is.null(x$nat_param))  cat("Nat. Param :", x$nat_param, "\n")
  if (!is.null(x$suff_stat))  cat("Suff. Stat :", x$suff_stat, "\n")
  if (!is.null(x$support))    cat("Support    :", x$support, "\n")
  invisible(x)
}

#' @export
summary.expfam <- function(object, theta_vals = NULL, ...) {
  print(object)
  if (!is.null(theta_vals) && !is.null(object$eta) && !is.null(object$A)) {
    cat("\nCanonical parameter and log-partition values:\n")
    df <- data.frame(
      theta = theta_vals,
      eta   = sapply(theta_vals, object$eta),
      A     = sapply(theta_vals, object$A)
    )
    print(df)
  }
  invisible(object)
}

#' Extract Sufficient Statistic from Exponential Family
#'
#' For an iid sample from an exponential family, returns the (minimal) sufficient statistic.
#'
#' @param ef An "expfam" object
#' @param x Numeric vector of observations
#'
#' @return Value of the sufficient statistic T(x)
#' @export
#'
#' @examples
#' ef <- exponential_family("poisson")
#' x <- rpois(20, lambda = 3)
#' exp_family_sufficient_stat(ef, x)
exp_family_sufficient_stat <- function(ef, x) {
  if (!inherits(ef, "expfam")) stop("ef must be of class 'expfam'")
  T_vals <- sapply(x, ef$T_stat)
  list(
    T_sum  = sum(T_vals),
    T_mean = mean(T_vals),
    n      = length(x),
    family = ef$family
  )
}

#' Check Full Rank Condition for Exponential Family
#'
#' Verifies whether an exponential family is of full rank by checking that
#' eta(theta) is not constant over the natural parameter space.
#'
#' @param ef An "expfam" object
#' @param theta_range Numeric vector of theta values to evaluate over
#'
#' @return Logical: TRUE if the family appears to be full rank
#' @export
#'
#' @examples
#' ef <- exponential_family("normal", param = list(sigma = 1))
#' check_full_rank(ef, theta_range = seq(-3, 3, by = 0.5))
check_full_rank <- function(ef, theta_range = seq(0.1, 5, by = 0.1)) {
  if (!inherits(ef, "expfam")) stop("ef must be of class 'expfam'")
  eta_vals <- tryCatch(sapply(theta_range, ef$eta), error = function(e) NULL)
  if (is.null(eta_vals)) {
    message("Could not evaluate eta over provided range.")
    return(NA)
  }
  is_full_rank <- length(unique(round(eta_vals, 10))) > 1
  cat("Full rank check for", ef$family, "family:\n")
  cat("  eta range: [", round(min(eta_vals), 4), ",", round(max(eta_vals), 4), "]\n")
  cat("  Full rank:", is_full_rank, "\n")
  invisible(is_full_rank)
}

# Helper: null-coalescing operator
`%||%` <- function(a, b) if (!is.null(a)) a else b
