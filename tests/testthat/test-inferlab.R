library(testthat)
library(inferlab)

test_that("exponential_family creates correct objects", {
  ef <- exponential_family("normal", param = list(sigma = 1))
  expect_s3_class(ef, "expfam")
  expect_equal(ef$family, "normal")

  ef2 <- exponential_family("poisson")
  expect_equal(ef2$family, "poisson")

  ef3 <- exponential_family("binomial", param = list(n = 10))
  expect_equal(ef3$family, "binomial")
})

test_that("sufficient_statistic returns correct values", {
  set.seed(1)
  x <- rnorm(20, mean = 3, sd = 1)
  res <- sufficient_statistic("normal_mean", data = x)
  expect_equal(res$value, mean(x))

  x_pois <- rpois(15, lambda = 2)
  res2 <- sufficient_statistic("poisson", data = x_pois)
  expect_equal(res2$value, sum(x_pois))
})

test_that("umvue_normal_mean is correct", {
  set.seed(42)
  x <- rnorm(30, mean = 5, sd = 2)
  res <- umvue_normal_mean(x, sigma = 2)
  expect_equal(res$estimate, mean(x))
  expect_equal(res$variance, 4/30)
})

test_that("umvue_poisson gives correct estimates", {
  set.seed(1)
  x <- rpois(20, lambda = 2)
  res1 <- umvue_poisson(x, estimand = "lambda")
  expect_equal(res1$estimate, mean(x))

  res2 <- umvue_poisson(x, estimand = "prob_zero")
  expect_equal(res2$estimate, (1 - 1/20)^sum(x))
})

test_that("fisher_information returns positive values", {
  fi_norm <- fisher_information("normal_mean", theta = 2, param = list(sigma = 1))
  expect_gt(fi_norm$I, 0)

  fi_pois <- fisher_information("poisson", theta = 3)
  expect_equal(fi_pois$I, 1/3)

  fi_binom <- fisher_information("binomial", theta = 0.4, param = list(m = 10))
  expect_equal(fi_binom$I, 10 / (0.4 * 0.6))
})

test_that("cramer_rao_bound is positive", {
  crb <- cramer_rao_bound("normal_mean", theta = 2, n = 20, g_deriv = 1,
                           param = list(sigma = 1))
  expect_equal(crb$CRB, 1/20)
})

test_that("mome gives correct MOM estimates", {
  set.seed(1)
  x <- rpois(100, lambda = 3)
  res <- mome(x, family = "poisson")
  expect_equal(res$estimates["lambda"], c(lambda = mean(x)))

  x_exp <- rexp(100, rate = 2)
  res2 <- mome(x_exp, family = "exponential")
  expect_equal(as.numeric(res2$estimates["lambda"]), 1/mean(x_exp))
})

test_that("np_lemma gives valid output", {
  set.seed(42)
  x <- rnorm(20, mean = 1.5, sd = 1)
  res <- np_lemma(x, "normal_mean", theta0 = 0, theta1 = 2, alpha = 0.05,
                   param = list(sigma = 1))
  expect_true(is.logical(res$reject))
  expect_gt(res$power, 0)
  expect_lt(res$power, 1 + 1e-10)
})

test_that("lrt gives valid statistics", {
  set.seed(1)
  x <- rpois(20, lambda = 3)
  res <- lrt(x, "poisson", theta0 = 2)
  expect_gt(res$lambda_stat, 0)
  expect_true(res$p_value >= 0 && res$p_value <= 1)
})

test_that("pivot_ci produces valid intervals", {
  set.seed(1)
  x <- rnorm(25, mean = 2, sd = 1)
  res <- pivot_ci(x, family = "normal_mean", conf_level = 0.95,
                   param = list(sigma = 1))
  expect_true(res$ci[1] < res$ci[2])
  expect_true(res$ci[1] < mean(x) && mean(x) < res$ci[2])
})

test_that("umau_ci is valid", {
  set.seed(42)
  x <- rpois(20, lambda = 3)
  res <- umau_ci(x, family = "poisson", conf_level = 0.95)
  expect_true(res$ci[1] < res$ci[2])
  expect_gt(res$ci[1], 0)
})

test_that("hoeffding_u_stat for mean matches sample mean", {
  set.seed(1)
  x <- rnorm(20, mean = 3, sd = 1)
  res <- hoeffding_u_stat(x, kernel_fn = function(xi) xi, m = 1, estimand = "mu")
  expect_equal(res$U, mean(x))
})
