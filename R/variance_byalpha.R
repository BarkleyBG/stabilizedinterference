
estimateVarianceByAlpha <- function(
  alphas,
  num_alphas,

  mu_alphas_ests,

  data,

  num_fixefs,
  fixefs,
  sigma,
  trt_model_obj,


  var_names,
  target_grids,
  # pop_mean_alphas_list,
  randomization_probability,
  weight_type,

  verbose,

  keep_components,
  compute_roots,
  integrate_alphas,
  deriv_control,
  contrast_type
){

  # average_treatment <-  mean(data[[var_names$treatment]])

  if (verbose){message('starting variance calcs')}


  # new_data <- data[names(data)%in%var_names]
  # new_data$model_mat <-
  #   split(ps_model_matrix , f = 1:NROW(new_data), drop = FALSE)
  ## variance estimates

  pop_mean_alphas_list <- list()
  for (alp_num in 1:num_alphas){

    mu_alpha_vec <- mu_alphas_ests[[alp_num]]
    theta_hat <- c(mu_alpha_vec, fixefs, sigma)

    alpha <- alphas[alp_num]
    # }
    # if (compute_roots){theta_hat <- NULL}
    geex_args_alpha <- list(
      estFUN = eeFunTV_IPTW,
      ##geex
      data = data,
      # data = new_data,
      units = var_names$grouping,
      compute_roots = compute_roots,
      roots = theta_hat,
      # root_control = setup_root_control(start = c(coef(trt_model_obj), 0)),
      # inner_args = list(
      #   weight_type = weight_type ### use hajek type here?
      # ),
      outer_args = list(
        # trt_model_obj = trt_model_obj,
        # outcome_var_name = outcome_var_name,
        # calcFunIPTW = calcFunIPTW
        trt_model_obj = trt_model_obj,
        num_fixefs = num_fixefs,
        var_names = var_names,
        alpha = alpha,
        integrate_alphas = integrate_alphas,
        x_levels = geex::grab(from = trt_model_obj, "design_levels"),
        randomization_probability  = randomization_probability,
        weight_type = weight_type
        # average_treatment = average_treatment


      )
    )
    if (!is.null(deriv_control)){
      geex_args_alpha$deriv_control <- deriv_control
    }

    # saveRDS(geex_args_alpha, file = quickLookup("geex_args_alpha.Rds"))
    # saveRDS(geex_args_alpha, file = quickLookup("geex_args_alpha_glm_hajek2.Rds"))

    variance_from_geex <- do.call(geex::m_estimate, geex_args_alpha)

    if (verbose){
      message(paste("done with alpha variance ", alp_num, "of", num_alphas))}

    pop_mean_alphas_list[[alp_num]] <- list(
      estimates = variance_from_geex@estimates,
      vcov = variance_from_geex@vcov
    )
    if (keep_components){
      pop_mean_alphas_list[[alp_num]]$sandwich_components <-
        variance_from_geex@sandwich_components
    }
    # pop_mean_alphas_list[[alp_num]] <- variance_from_geex
  }

  ## #' @param pop_means_grid "difference" or others (perhaps, later).
  ## #' @param effects_grid "difference" or others (perhaps, later).
  ## #' @param num_alphas number of allocations

  #
  tidy_args <- list(
    target_grids = target_grids,
    pop_mean_alphas_list = pop_mean_alphas_list,
    # num_alphas,
    contrast_type = contrast_type
  )
  # saveRDS(tidy_args, file = quickLookup("test_tidyTargets_glm.Rds"))

  target_ests <- do.call(tidyTargets,tidy_args)

  output <- list(
    target_ests = target_ests,
    pop_mean_alphas_list = pop_mean_alphas_list
  )
  output
}


# options(error='recover')

#' Estimating Function for IPTW
#'
#' This function is to be passed into geex::m_estimate
#'
#' @inheritParams estimateTV_IPTW
#' @param alpha One of the allocations at a time
#' @param num_fixefs Number of fixed effect parameters from treatment model.
#'   Perhaps unncessaary coding.
#' @param var_names A list of names for outcome, treatment, clustering, and
#'   perhaps participation.
#' @param trt_model_obj The fitted model object (usually a glm).
#' @param x_levels default NULL unless there are factos in design matrix. From
#'   \code{\link[geex]{grab_design_levels}}.
#' @param randomization_probability usually 1. e.g. 2/3  in Perez-Heydrich et
#'   al. (2014) Biometrics
#' @param integrate_alphas true
#'
#' @export
eeFunTV_IPTW <- function(
  data,
  trt_model_obj,
  num_fixefs,
  var_names,
  alpha,
  x_levels = NULL,
  integrate_alphas ,#= integrate_allocations,
  randomization_probability ,#= randomization_probability,
  weight_type ##HT, Hajek1, Hajek2
  # average_treatment
  # outcome_var_name,
  # calcFunTVIPTW
){



  # if ( "glm" %in% class(trt_model_obj) ){
  model_matrix <- geex::grab_design_matrix(
    geex::grab_fixed_formula(trt_model_obj),
    # rhs_formula = trt_model_obj,
    xlev = x_levels,
    data = data
  )
  # model_matrix <- do.call(rbind, data$model_mat)
  # model_matrix <- stats::model.matrix(trt_model_obj$formula, data = data)
  # } else
  #   if ("glmerMod" %in% class(trt_model_obj) ) {
  # model_matrix <- lme4::getME(trt_model_obj, "X", data = data)
  # model_matrix <- geex::grab_design_matrix(
  #   geex::grab_fixed_formula(trt_model_obj),
  #   data = data
  #   )
  # geex::gra
  # }
  ## Put this at the end
  closureModel <- geex::grab_psiFUN(
    data=data,
    object = trt_model_obj,
    xlev = x_levels
  )
  ## Or perhaps re-do data to be a split list with these pre-specified.

  # # treatment <- stats::model.response(stats::model.frame(trt_model_obj$formula, data = data))
  # # model_DV <- data[[var_names$model_DV]]
  outcome <- data[[var_names$outcome]]
  treatment <- data[[var_names$treatment]]
  # participation <- data[[var_names$participation]] ##perhaps NULL

  ind_z1 <- treatment==1
  ind_z0 <- treatment==0


  # sum_trt <- sum(treatment)
  # clust_size <- length(treatment)



  ## Estimating function for IPTW (unstabilized)
  closureTV_IPTW <- function(theta){

    mu1 <- theta[1]
    mu0 <- theta[2]
    mu_marg <- theta[3]
    fixefs <- theta[3+(1:num_fixefs)]
    if ("glmerMod" %in% class(trt_model_obj) ) {
      sigma <- theta[4+num_fixefs]

    } else
      if ( "glm" %in% class(trt_model_obj) ) {


        sigma <- NULL

      } else {
        stop("treatment model type not recognized")
      }


    ipw_args <- list(
      # logit_integrand_fun,
      fixefs = fixefs,
      treatment = treatment,
      sigma = sigma, ## perhaps NULL
      # participation = participation,
      model_matrix = model_matrix,
      integrate_alphas = integrate_alphas,
      randomization_probability = randomization_probability,
      alpha = alpha
    )

    # saveRDS(ipw_args, file = quickLookup("calcIPW_avgtrt_args.Rds"))
    # ipw_args$alpha <- average_treatment
    # pi_ipw_basic <- do.call(calcPiIPW, ipw_args)
    # alpha_here <- alpha
    # alpha_ratio <- alpha_here/average_treatment
    # alpha_inv_ratio <- (1-alpha_here)/(1-average_treatment)
    # pi_ratio <- (alpha_ratio)^(sum_trt)*(alpha_inv_ratio)^(clust_size - sum_trt)
    # pi_ipw <- pi_ratio*pi_ipw_basic

    pi_ipw <- do.call(calcPiIPW, ipw_args)

    if (weight_type == "HT") {
      pi_ipw_Y <- outcome*pi_ipw


      # sum1 <- sum( pi_ipw_Y[ind_z1] )/alpha
      # sum0 <- sum( pi_ipw_Y[ind_z0] )/(1-alpha)
      # out_vec <- c(
      #   sum1 - mu1, ## for mu(1,\alpha)
      #   sum0 - mu0,## for mu(0,\alpha)
      #   (alpha*sum1 + (1-alpha)*sum0) - mu_marg ## for mu( \alpha) marginal
      # )

      out_vec <- c(
        sum( pi_ipw_Y*ind_z1/(alpha) - mu1 ) , ## for mu(1,\alpha)
        sum( pi_ipw_Y*ind_z0/(1-alpha) - mu0 ) , ## for mu(1,\alpha)
        sum( pi_ipw_Y - mu_marg )   ## for mu( \alpha) marginal
      )

    } else
      if (weight_type == "Hajek1") {
        ## calculate individual prop scores
        indiv_prop_scores <- rep(NA, length(treatment))

        for (jj in 1:length(indiv_prop_scores)){
          indiv_prop_scores[jj] <-
            calcPiIPW(
              # logit_integrand_fun,
              fixefs = fixefs,
              treatment = treatment[jj],
              sigma = sigma, ## perhaps NULL
              # participation = participation[jj],
              model_matrix = model_matrix[jj, , drop = FALSE],
              integrate_alphas = integrate_alphas,
              randomization_probability = randomization_probability,
              alpha = 1######alpha is not used in individual prop scores
            )

        }
        pi_ipw_Y <- outcome*pi_ipw

        out_vec <- c(
          (pi_ipw_Y[ind_z1]/alpha) - ( mu1 / (indiv_prop_scores[ind_z1]) ),
          (pi_ipw_Y[ind_z0]/(1-alpha)) - ( mu0 / (indiv_prop_scores[ind_z0]) ),
          pi_ipw_Y  - ( mu_marg / (indiv_prop_scores))
        )


      } else
        if (weight_type == "Hajek2") {

          out_vec <-  c(
            sum(outcome[ind_z1] - mu1) * pi_ipw/alpha,
            sum(outcome[ind_z0] - mu0) * pi_ipw/(1-alpha),
            sum(outcome - mu_marg) * pi_ipw
          )
        }
    out_vec
  }

  closureStacked <- function(theta){

    out_full <- c(
      closureTV_IPTW(theta), ## Returns vector of length 3 for: mu1, mu0, mu_marg
      closureModel(theta[-(1:3)]) ## Returns vector of len=length(theta)-1
    ) ## Returns a vector of len=length(theta)


    out_full

  }
}
#
