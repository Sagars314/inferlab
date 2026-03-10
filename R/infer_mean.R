#' Infer Mean using t-test
#'
#' Performs a one-sample t-test and returns a tidy result.
#'
#' @param x Numeric vector
#' @param mu Hypothesized mean
#'
#' @return A list containing test statistic, p-value and confidence interval
#' @export
infer_mean <- function(x, mu = 0) {

  test <- t.test(x, mu = mu)

  list(
    statistic = test$statistic,
    p_value = test$p.value,
    conf_int = test$conf.int
  )
}
