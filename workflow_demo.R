###############################################################
#  MTH209 PACKAGE — FULL WORKFLOW DEMO
#  Run this script to see the ENTIRE package in action
###############################################################

# Install the package first (run once):
# install.packages("devtools")
# devtools::install("path/to/MTH209")
# Then load it:
# library(MTH209)

cat("\n")
cat("##############################################################\n")
cat("##   MTH209 — ESTIMATION THEORY & MATHEMATICAL INFERENCE   ##\n")
cat("##            Full Package Workflow Demo                    ##\n")
cat("##############################################################\n\n")

# ==== CHAPTER 1: Exponential Families ====
cat("\n========== CHAPTER 1: EXPONENTIAL FAMILIES ==========\n\n")
exp_family_check("Normal", params = list(mu="unknown", sigma2="known"))
cat("\n")
exp_family_canonical("Poisson")
cat("\n")
exp_family_full_rank("Gamma")

# ==== CHAPTER 2: Sufficiency ====
cat("\n========== CHAPTER 2: SUFFICIENCY ==========\n\n")
sufficiency_nf_factorisation("Poisson(lambda)", "sum(X_i)")
cat("\n")
sufficiency_minimal("Normal(mu, sigma^2)", "( sum(X_i), sum(X_i^2) )")
cat("\n")
sufficiency_ancillary("Normal(mu, 1)", "X_1 - X_2")

# ==== CHAPTER 3: Completeness ====
cat("\n========== CHAPTER 3: COMPLETENESS ==========\n\n")
completeness_family("Poisson family, lambda > 0")
cat("\n")
completeness_statistic("Binomial(n, p)", "sum(X_i)")
cat("\n")
completeness_basu("Normal(mu, 1)", "X_bar", "S^2 / X_bar")
cat("\n")
completeness_rao_blackwell("X_1", "sum(X_i)", "Bernoulli(p)")

# ==== CHAPTER 4: Unbiasedness & UMVUE ====
cat("\n========== CHAPTER 4: UNBIASEDNESS & UMVUE ==========\n\n")
unbiased_check("X_bar", "mu", "Normal(mu, sigma^2)")
cat("\n")
umvue_find("p^2", "Binomial(n, p)")
cat("\n")
umvue_lehmann_scheffe("(sum(X_i)*(sum(X_i)-1)) / (n*(n-1))", "sum(X_i)", "Binomial(n,p)")

# ==== CHAPTER 5: UMVUE Methods ====
cat("\n========== CHAPTER 5: METHODS FOR FINDING UMVUE ==========\n\n")
umvue_method_direct("e^{-lambda}", "sum(X_i)", "Poisson(lambda)")
cat("\n")
umvue_method_rao_blackwell("I(X_1 = 0)", "sum(X_i)", "Poisson(lambda)")
cat("\n")
umvue_ustatistic("I(x_1 <= x_2)", degree=2, parameter="P(X1 <= X2)", n=20)

# ==== CHAPTER 6: Information Bounds ====
cat("\n========== CHAPTER 6: INFORMATION INEQUALITIES ==========\n\n")
info_fisher("Poisson(lambda)", "lambda")
cat("\n")
info_fisher_matrix("Bivariate Normal", c("mu_1", "mu_2", "sigma_11", "sigma_12", "sigma_22"))
cat("\n")
info_cramer_rao("X_bar", "mu", "Normal(mu, sigma^2=1)")
cat("\n")
info_hcr_bound("X_bar", "mu", "Uniform(0, theta)")
cat("\n")
info_bhattacharya("Normal(0, sigma^2)", "sigma^2", order=3)

# ==== CHAPTER 7: Methods of Estimation ====
cat("\n========== CHAPTER 7: METHODS OF ESTIMATION ==========\n\n")
estimation_mome("Gamma(alpha, beta)", c("alpha", "beta"))
cat("\n")
estimation_mle("Normal(mu, sigma^2)", c("mu", "sigma^2"))
cat("\n")
estimation_minmse("Linear estimators aX_bar + b", "mu", "Normal(mu, sigma^2)")

# ==== CHAPTER 8: Hypothesis Testing Basics ====
cat("\n========== CHAPTER 8: HYPOTHESIS TESTING BASICS ==========\n\n")
test_define_hypothesis("mu = 0", "mu > 0")
cat("\n")
test_errors_and_power("mu = 0", "mu = 1", "X_bar", "X_bar > 1.645/sqrt(n)", "Normal(mu,1)")

# ==== CHAPTER 9: NP Lemma & Optimal Tests ====
cat("\n========== CHAPTER 9: NEYMAN-PEARSON & OPTIMAL TESTS ==========\n\n")
np_lemma("Normal(mu=0, sigma=1)", "Normal(mu=1, sigma=1)", alpha=0.05)
cat("\n")
test_ump("lambda <= 1", "lambda > 1", "Poisson(lambda)", alpha=0.05)
cat("\n")
test_umpu("mu = 0", "mu != 0", "Normal(mu, sigma^2=1)", alpha=0.05)
cat("\n")
test_mlr("Exponential(lambda)", "sum(X_i)", "lambda <= 1", "lambda > 1", alpha=0.05)

# ==== CHAPTER 10: Testing in Exponential Families ====
cat("\n========== CHAPTER 10: TESTING IN EXPONENTIAL FAMILIES ==========\n\n")
expfam_test_onesided("Poisson(lambda)", "lambda <= 2", "lambda > 2", alpha=0.05)
cat("\n")
expfam_test_twosided("Normal(mu, 1)", "mu = 0", "mu != 0", alpha=0.05)
cat("\n")
expfam_test_neyman_structure("Normal(mu, sigma^2)", "mu = 0", "sigma^2", alpha=0.05)
cat("\n")
test_likelihood_ratio("Normal(mu, sigma^2)", "mu = 0, sigma^2 = 1",
                      "mu != 0 or sigma^2 != 1", alpha=0.05)

# ==== CHAPTER 11: Confidence Intervals ====
cat("\n========== CHAPTER 11: CONFIDENCE INTERVALS ==========\n\n")
ci_pivot("Normal(mu, 1)", "mu", "(X_bar - mu)*sqrt(n)", conf_level=0.95)
cat("\n")
ci_shortest("Normal(mu, sigma^2)", "mu", conf_level=0.95)
cat("\n")
ci_uma("Poisson(lambda)", "lambda", conf_level=0.95)
cat("\n")
ci_umau("Normal(mu, sigma^2)", "mu", conf_level=0.95)

cat("\n\n")
cat("##############################################################\n")
cat("##  END OF WORKFLOW DEMO — All 11 modules covered!         ##\n")
cat("##  Next step: replace each cat() stub with real maths :)  ##\n")
cat("##############################################################\n")
