
#' Calcluates point estimates for target estimands
#'
#' This function calculates the estimated values in a direct plug-in method.
#'
#' @inheritParams estimateTV_IPTW
#' @param num_alphas The number of allocation parameters
#' @param num_clusters The number of unique clusters (i.e., i.i.d. sample units)
#' @param unique_clusters The ID values for the unique clusters
#' @param grouping_vector The vector of cluster ID for all individuals
#' @param treatment_vector The vector of treatment identifiers for all
#'   individuals
#' @param outcome_vector The vector of observed outcome for all individuals
#' @param ps_model_matrix The matrix of pre-treatment variables in the
#'   propensity score model for all individuals
#' @param fixefs The estimated values of the fixed effects parameters from the
#'   propensity score model
#' @param sigma The estimated value of the (single) random effect variance
#'   component from the propensity score model
#'
estimateTargets <- function(
  alphas, num_alphas,
  num_clusters, unique_clusters, grouping_vector,

  treatment_vector,
  ps_model_matrix,
  outcome_vector,

  fixefs,
  sigma,

  integrate_alphas,
  randomization_probability,
  weight_type
){

  pi_ipw_alphas <- list()
  mu_alphas_ests <- list()

  for(alp_num in 1:num_alphas) {
    alpha <- alphas[alp_num]


    pi_ipw_vec <- rep(NA, num_clusters)
    summand_mat <- matrix(NA, ncol = 3, nrow = num_clusters)
    denom_mat <- matrix(NA, ncol = 3, nrow = num_clusters)

    for (cluster_number in 1:num_clusters){

      this_cluster <- unique_clusters[cluster_number]
      cluster_indices <- (grouping_vector == this_cluster)

      treatment <- treatment_vector[cluster_indices]
      model_matrix <- ps_model_matrix[cluster_indices, , drop=FALSE]
      outcome <- outcome_vector[cluster_indices]
      ind_z1 <- treatment==1
      ind_z0 <- treatment==0

      pi_ipw  <- calcPiIPW(
        treatment = treatment,
        model_matrix = model_matrix,
        fixefs = fixefs,
        sigma = sigma, ## perhaps NULL (for GLM's)
        integrate_alphas = integrate_alphas,
        randomization_probability = randomization_probability,
        alpha = alpha
      )
      pi_ipw_vec[cluster_number] <- pi_ipw


      pi_ipw_Y <- outcome*pi_ipw

      ## adjusting for pi(A_i_notj, alpha) for conditional estimands
      sum_y_ipw1 <- sum( pi_ipw_Y[ind_z1] )/alpha
      sum_y_ipw0 <- sum( pi_ipw_Y[ind_z0] )/(1-alpha)

      summand_mat[cluster_number, ] <- c(
        sum_y_ipw1, sum_y_ipw0, (alpha*sum_y_ipw1 + (1-alpha)*sum_y_ipw0)
      )

      ## The conditional treatments are the sum times alpha

      if (weight_type == "HT") {
        denom_mat[cluster_number, ] <- length(treatment)
      } else
        if (weight_type == "Hajek1") {
          ## calculate individual prop scores
          indiv_prop_scores_ipw <- rep(NA, length(treatment))

          for (jj in 1:length(indiv_prop_scores_ipw)){
            indiv_prop_scores_ipw[jj] <-
              calcPiIPW(
                treatment = treatment[jj],
                model_matrix = model_matrix[jj, , drop = FALSE],
                fixefs = fixefs,
                sigma = sigma, ## perhaps NULL
                integrate_alphas = FALSE, ##ALPHA NOT USED HERE
                randomization_probability = randomization_probability,
                alpha = NULL ##alpha is not used in individual prop scores
              )

          }

          sum_ipw_1 <- sum(indiv_prop_scores_ipw[ind_z1])
          sum_ipw_0 <- sum(indiv_prop_scores_ipw[ind_z0])
          denom_mat[cluster_number, ] <- c(sum_ipw_1, sum_ipw_0, sum_ipw_1 + sum_ipw_0)

        } else
          if (weight_type == "Hajek2") {

            sum_ipw_1 <- sum(ind_z1)*pi_ipw/alpha
            sum_ipw_0 <- sum(ind_z0)*pi_ipw/(1-alpha)
            denom_mat[cluster_number, ] <-
              c(sum_ipw_1, sum_ipw_0, alpha*sum_ipw_1 + (1-alpha)*sum_ipw_0)
          } else {
            stop("weight_type not recognized")
          }
    } ## looping over cluster_number
    stopifnot(!any(is.na(pi_ipw_vec)))
    stopifnot(!any(is.na(denom_mat)))
    stopifnot(!any(is.na(summand_mat)))

    denom_vec <- colSums(denom_mat)
    summand_vec <- colSums(summand_mat)
    mu_ests_vec <- summand_vec / denom_vec

    mu_alphas_ests[[alp_num]] <- mu_ests_vec
    pi_ipw_alphas[[alp_num]] <- pi_ipw_vec

    names(mu_alphas_ests)[alp_num] <- alpha
    names(pi_ipw_alphas)[alp_num] <- alpha
  } ## looping over alphas

  cluster_propensity_scores <- pi_ipw_alphas

  out_list <- list(
    mu_alphas_ests = mu_alphas_ests,
    pi_ipw_alphas = pi_ipw_alphas,
    cluster_propensity_scores = cluster_propensity_scores
  )
}
