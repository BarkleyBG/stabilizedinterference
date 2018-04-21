
#' Estimates of the asymptotic variance of the estimators
#'
#' These functions carry out the M-estimation procedures for the "sandwich"
#' variance estimators
#'
#' @inheritParams estimateEffects
#' @inheritParams estimateTargets
#' @inheritParams getModel
#' @inheritParams tidyTargetsCombined
#' @inheritParams eeFunCombined
#' @param mu_alphas_ests Point estimates of population mean target estimands
#' @param num_fixefs The number of fixed effect parameters from treatment model
#' @param deriv_control Optional for \code{\link[geex]{m_estimate}}
#' @param verbose Optional argument from \code{\link{estimateEffects}}
#' @param keep_components Optional argument from \code{\link{estimateEffects}}
#' @param compute_roots Optional argument from \code{\link{estimateEffects}}
#'
estimateVarianceCombined <- function(
  alphas,
  num_alphas,

  mu_alphas_ests,

  data,

  num_fixefs,
  fixefs,
  sigma,
  x_levels,
  trt_model_obj,


  var_names,
  target_grids,
  # pop_mean_alphas_list,
  randomization_probability,
  weight_type,
  # average_treatment,

  verbose,

  keep_components,
  compute_roots,
  integrate_alphas,
  deriv_control,
  contrast_type
){


  if (verbose){message('starting variance calcs')}


  # new_data <- data[names(data)%in%var_names]
  # new_data$model_mat <-
  #   split(ps_model_matrix , f = 1:NROW(new_data), drop = FALSE)
  ## variance estimates

  # mu_alpha_array <- array(
  #   NA,
  #   dim = c(
  #     length(unique(data[[var_names$grouping]])),
  #     3,
  #     num_alphas
  #   )
  # )
  if (any(c("lm", "glm", "glmerMod") %in% class(trt_model_obj))){
    theta_hat <- c(unlist(mu_alphas_ests), fixefs, sigma)

  } else
    if (
      length(class(trt_model_obj))==1 &&
      class(trt_model_obj) == "list" ) {
      # if (model_method == "oracle"){
      theta_hat <- c(unlist(mu_alphas_ests))#, fixefs, sigma)

    }

  if (compute_roots){
    roots <-  NULL
    # root_control <-  geex::setup_root_control(start = theta_hat)
  } else{
    roots <- theta_hat
    # root_control <- NULL
  }
  average_treatment <-  mean(data[[var_names$treatment]])

  geex_args_alphas <- list(
    estFUN = eeFunCombined,
    ##geex
    data = data,
    units = var_names$grouping,
    compute_roots = compute_roots,
    roots = roots,
    # root_control = root_control,
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
      alphas = alphas,
      num_alphas = num_alphas,
      integrate_alphas = integrate_alphas,
      x_levels = x_levels,
      randomization_probability  = randomization_probability,
      weight_type = weight_type,
      average_treatment = average_treatment


    )
  )

  if (!any(c("lm", "glm", "glmerMod") %in% class(trt_model_obj)) &&
      length(class(trt_model_obj))==1 &&
      class(trt_model_obj) == "list" ) {
  # if (model_method == "oracle"){
    geex_args_alphas$outer_args$oracle_fixefs <- fixefs
    geex_args_alphas$outer_args$oracle_sigma <- sigma
  }
  if (!is.null(deriv_control)){
    geex_args_alphas$deriv_control <- deriv_control
  }
  if (compute_roots){
    # roots <-  NULL
    geex_args_alphas$root_control <-
      geex::setup_root_control(start = theta_hat)
  }
  # message("about to geex")

  # saveRDS(geex_args_alpha, file = quickLookup("geex_args_alpha.Rds"))
  # saveRDS(geex_args_alpha, file = quickLookup("geex_args_alpha_glm_hajek2.Rds"))

  geex_output <- do.call(geex::m_estimate, geex_args_alphas)

  alpha_grid <- data.frame(
    alpha1_num  = rep(1:num_alphas, each = 3),
    alpha1 = rep(alphas, each = 3),
    trt1 = rep(c(1,0,NA), num_alphas),
    alpha2_num = NA,
    alpha2 = NA,
    trt2 = NA
  )
  # ests <- geex_output@roots
  alpha_grid$estimate <- geex_output@estimates[1:(3*num_alphas)]
  alpha_grid$variance <- diag(geex_output@vcov)[1:(3*num_alphas)]

  tidy_grid_args <- list(
  # tidy_args <- list(
    target_grids = target_grids,
    alpha_grid = alpha_grid,
    geex_output = geex_output,
    contrast_type = contrast_type
  )
  # saveRDS(tidy_grid_args, file = quickLookup("test_tidyTargetsCombined_glm.Rds"))

  target_ests <- do.call(tidyTargetsCombined,tidy_grid_args)
  #
  output <- list(
    target_ests = target_ests,
    geex_object = geex_output
  )
  output
}


# options(error='recover')

#' Estimating Function for IPTW
#'
#' This function is to be passed into geex::m_estimate
#'
#' @inheritParams estimateEffects
#' @inheritParams estimateTargets
#' @inheritParams getModel
#' @param num_fixefs Number of fixed effect parameters from treatment model.
#'   Perhaps unncessaary coding.
#' @param var_names A list of names for outcome, treatment, clustering, and
#'   perhaps participation.
#' @param x_levels Default NULL unless there are factos in design matrix. From
#'   \code{\link[geex]{grab_design_levels}}.
#' @param average_treatment e.g. 0.37
#' @param oracle_fixefs When \code{model_method=="Oracle"} then the oracle
#'   values of fixed effects are passed through this argument. Default
#'   \code{NULL}.
#' @param oracle_sigma When \code{model_method=="Oracle"} then the oracle values
#'   of the random effect component is passed through this argument. Default
#'   \code{NULL}.
#'
#' @export
eeFunCombined <- function(
  data,
  trt_model_obj,
  num_fixefs,
  var_names,
  alphas,
  num_alphas,
  x_levels,
  integrate_alphas ,#= integrate_allocations,
  randomization_probability ,#= randomization_probability,
  weight_type, ##HT, Hajek1, Hajek2
  average_treatment,
  oracle_fixefs = NULL,
  oracle_sigma = NULL
){


  # dots <- list(...)
  model_class <- class(trt_model_obj)
 if (any(c("lm", "glm", "glmerMod") %in% model_class)) {
  # if ( "glm" %in% class(trt_model_obj) ){
  model_matrix <- geex::grab_design_matrix(
    geex::grab_fixed_formula(trt_model_obj),
    # rhs_formula = trt_model_obj,
    xlev = x_levels,
    data = data
  )

    ## Put this at the end
    closureModel <- geex::grab_psiFUN(
      data = data,
      object = trt_model_obj,
      xlev = x_levels
    )
 } else
   if ( length(model_class)==1 && model_class == "list" ) {

   # if (length(class(trt_model_obj))==1 && class(trt_model_obj)=="list"){
   closureModel <- function(x){NULL}
   model_matrix <- geex::grab_design_matrix(
     trt_model_obj$modeling_formula,
     # rhs_formula = trt_model_obj,
     xlev = x_levels,
     data = data
   )
   fixefs <- oracle_fixefs
   sigma <- oracle_sigma

 }  else {stop("model_class not recognized")}

  ## Or perhaps re-do data to be a split list with these pre-specified.

  # # treatment <- stats::model.response(stats::model.frame(trt_model_obj$formula, data = data))
  # # model_DV <- data[[var_names$model_DV]]
  outcome <- data[[var_names$outcome]]
  treatment <- data[[var_names$treatment]]
  # participation <- data[[var_names$participation]] ##perhaps NULL

  ind_z1 <- treatment==1
  ind_z0 <- treatment==0

  ## Estimating function for IPTW (unstabilized)
  closureTV_IPTW <- function(theta){

    if (any(c("lm", "glm", "glmerMod") %in% model_class)) {

      fixefs <- theta[(3*num_alphas)+(1:num_fixefs)]
      if ("glmerMod" %in% class(trt_model_obj) ) {
        sigma <- theta[ (3*num_alphas) +num_fixefs + 1]
        stopifnot(sigma>0)
      } else
        if ( "lm" %in% class(trt_model_obj) ) {
          sigma <- NULL
        } else {
          stop("treatment model type not recognized")
        }
    } else
      if ( length(model_class)==1 && model_class == "list" ) {
      ##nothing
    }

    pi_ipw_basic <- calcPiIPW(
      # logit_integrand_fun,
      fixefs = fixefs,
      treatment = treatment,
      sigma = sigma, ## perhaps NULL
      # participation = participation,
      model_matrix = model_matrix,
      integrate_alphas = integrate_alphas,
      randomization_probability = randomization_probability,
      alpha = average_treatment
    )

    # pi_ipws <- lapply(alphas, function(alpha){
    #   alpha_ratio <- alpha/average_treatment
    #   alpha_inv_ratio <- (1-alpha)/(1-average_treatment)
    #   pi_ratio <- (alpha_ratio)^treatment*(alpha_inv_ratio)^(1-treatment)
    #   pi_ipw <- pi_ratio*pi_ipw_basic
    # })
    sum_trt <- sum(treatment)
    clust_size <- length(treatment)

    # mu_list <- list()
    out_vec_list <- list()
    for (alpha_num in 1:num_alphas){

      alpha_here <- alphas[alpha_num]
      alpha_ratio <- alpha_here/average_treatment
      alpha_inv_ratio <- (1-alpha_here)/(1-average_treatment)
      pi_ratio <- (alpha_ratio)^(sum_trt)*(alpha_inv_ratio)^(clust_size - sum_trt)
      pi_ipw <- pi_ratio*pi_ipw_basic


      mu_vec_here <- theta[ (alpha_num-1)*3 + (1:3)]
      # mu_list[[alpha_num]] <- mu_vec_here

      if (weight_type =="HT"){
        pi_ipw_Y <- outcome*pi_ipw
        out_vec <- c(
          sum( pi_ipw_Y*ind_z1/(alpha_here) - mu_vec_here[1] ) , ## for mu(1,\alpha)
          sum( pi_ipw_Y*ind_z0/(1-alpha_here) - mu_vec_here[2] ) , ## for mu(1,\alpha)
          sum( pi_ipw_Y - mu_vec_here[3] )   ## for mu( \alpha) marginal
        )
        out_vec_list[[alpha_num]] <- out_vec
      } else
        if (weight_type == "Hajek1"){
          stop("Hajek1 not implemented yet")
        } else
          if (weight_type == "Hajek2"){
            out_vec <-  c(
              sum(outcome[ind_z1] - mu_vec_here[1]) * pi_ipw/alpha_here,
              sum(outcome[ind_z0] - mu_vec_here[2]) * pi_ipw/(1-alpha_here),
              sum(outcome - mu_vec_here[3]) * pi_ipw
            )
            out_vec_list[[alpha_num]] <- out_vec
          } else
            if (weight_type == "HT_TV"){
              pi_ipw_Y <- outcome*pi_ipw
              out_vec <- c(
                mean( pi_ipw_Y*ind_z1/(alpha_here) - mu_vec_here[1] ) , ## for mu(1,\alpha)
                mean( pi_ipw_Y*ind_z0/(1-alpha_here) - mu_vec_here[2] ) , ## for mu(1,\alpha)
                mean( pi_ipw_Y - mu_vec_here[3] )   ## for mu( \alpha) marginal
              )
              out_vec_list[[alpha_num]] <- out_vec
            } else {stop("Please supply weight type = HT, Hajek1 or Hajek2")}

    }


    long_output_vector <- unlist(out_vec_list)
    long_output_vector
  }

  closureStacked <- function(theta){

    out_full <- c(
      closureTV_IPTW(theta), ## Returns vector of length alpha*3 for: mu1, mu0, mu_marg
      closureModel( theta[-(1:(3*num_alphas))] ) ## Returns vector of len= num_fixefs + 1*(has_sigma)
    ) ## Returns a vector of len=length(theta)
    out_full
  }
}
#
