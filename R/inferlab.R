#' inferlab: Statistical Inference Toolkit
#'
#' A comprehensive R package for classical statistical inference,
#' covering the full MTH418 curriculum.
#'
#' @section Exponential Families:
#' \itemize{
#'   \item \code{\link{exponential_family}} - Create exponential family objects in canonical form
#'   \item \code{\link{exp_family_sufficient_stat}} - Extract sufficient statistics
#'   \item \code{\link{check_full_rank}} - Verify full rank condition
#' }
#'
#' @section Sufficiency:
#' \itemize{
#'   \item \code{\link{sufficient_statistic}} - Sufficient statistics for common families
#'   \item \code{\link{neyman_fisher_factorization}} - Verify Neyman-Fisher factorization
#'   \item \code{\link{minimal_sufficient}} - Minimal sufficiency via ratio criterion
#'   \item \code{\link{ancillary_statistic_test}} - Test ancillarity via Monte Carlo
#' }
#'
#' @section Completeness:
#' \itemize{
#'   \item \code{\link{completeness_check}} - Completeness assessment
#'   \item \code{\link{basu_theorem_demo}} - Demonstrate Basu's theorem
#'   \item \code{\link{rao_blackwell}} - Rao-Blackwell improvement
#' }
#'
#' @section UMVUE Estimation:
#' \itemize{
#'   \item \code{\link{umvue_normal_mean}} - UMVUE for normal mean
#'   \item \code{\link{umvue_normal_variance}} - UMVUE for normal variance
#'   \item \code{\link{umvue_exponential}} - UMVUE for exponential parameters
#'   \item \code{\link{umvue_binomial}} - UMVUE for binomial parameter functions
#'   \item \code{\link{umvue_poisson}} - UMVUE for Poisson parameter functions
#'   \item \code{\link{lehmann_scheffe}} - Lehmann-Scheffe theorem application
#'   \item \code{\link{hoeffding_u_stat}} - Hoeffding's U-statistics
#' }
#'
#' @section Information Inequalities:
#' \itemize{
#'   \item \code{\link{fisher_information}} - Fisher information for common families
#'   \item \code{\link{fisher_information_matrix}} - Fisher information matrix (multiparameter)
#'   \item \code{\link{cramer_rao_bound}} - Cramer-Rao lower bound
#'   \item \code{\link{cramer_rao_multiparameter}} - CRB for multiparameter case
#'   \item \code{\link{hcr_bound}} - Hammersley-Chapman-Robbins bound
#'   \item \code{\link{bhattacharyya_bounds}} - Bhattacharyya system of bounds
#' }
#'
#' @section Methods of Estimation:
#' \itemize{
#'   \item \code{\link{mome}} - Method of Moments Estimation
#'   \item \code{\link{min_mse_estimator}} - Minimum MSE estimator
#' }
#'
#' @section Hypothesis Testing:
#' \itemize{
#'   \item \code{\link{np_lemma}} - Neyman-Pearson Most Powerful test
#'   \item \code{\link{ump_test}} - Uniformly Most Powerful test
#'   \item \code{\link{umpu_test}} - Uniformly Most Powerful Unbiased test
#'   \item \code{\link{mlr_test}} - Monotone Likelihood Ratio test
#'   \item \code{\link{lrt}} - Likelihood Ratio Test
#'   \item \code{\link{power_function}} - Power function computation and plotting
#'   \item \code{\link{neyman_pearson_region}} - NP critical region
#' }
#'
#' @section Confidence Intervals:
#' \itemize{
#'   \item \code{\link{pivot_ci}} - Pivot-based confidence intervals
#'   \item \code{\link{shortest_ci}} - Shortest expected length CI
#'   \item \code{\link{uma_ci}} - Uniformly Most Accurate CI
#'   \item \code{\link{umau_ci}} - Uniformly Most Accurate Unbiased CI
#' }
#'
#' @docType package
#' @name inferlab
"_PACKAGE"
