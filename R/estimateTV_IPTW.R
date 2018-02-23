
#
#' Estimate ATE with IPTW
#   #'
#' Fits a parametric model and estimates ATE via IPTW with Wald-type confidence
#' intervals from the empirical sandwich standard error estimates.
#'
#' @param formula Three-part formula: Outcome | Treatment ~ model_predictors | cluster_ID. Will be coerced to object of type Formula.
#' @param data the dataframe. Will be coerced from "tbl_df" to data.frame.
#' @param model_method currently only supported "logistic" for logit-link binomial GLM.
#' @param weight_type Currently only supports "unstabilized"
#' @param ... additional args
#' @param alphas the range of allocations or policies from 0 to 1.
#'
#' @export
estimateTV_IPTW <- function(
  data,
  formula,
  alphas,
  weight_type  = c("HT", "Hajek1", "Hajek2")[3],
  model_method  = c("glm", "glmer")[2],
  model_options = list(nAGQ = 5, family = "binomial"),
  ...
){
  # stop("take geex object out of the loop; save space!")

  dots <- list(...)
  dots_names <- names(dots)
  randomization_probability <-
    ifelse("randomization_probability" %in% dots_names,
           dots$randomization_probability, 1)
  compute_roots <-
    ifelse("compute_roots" %in% dots_names,
           dots$compute_roots, FALSE)
  integrate_alphas <-
    ifelse("integrate_alphas" %in% dots_names,
           dots$integrate_alphas, TRUE)
  contrast_type <-
    ifelse("contrast_type" %in% dots_names,
           dots$contrast_type, "difference")
  verbose <-
    ifelse("verbose" %in% dots_names,
           dots$verbose, FALSE)
  if (contrast_type != "difference"){stop("contrast_type must be 'difference'")}
  if ("target_grids" %in% dots_names){
    target_grids <- dots$target_grids
  } else {
    target_grids <- makeTargetGrids(alphas = alphas)
  }
  if ("deriv_control" %in% dots_names){
    deriv_control <- dots$deriv_control
  } else {
    deriv_control <- NULL
  }

  # target_grids <- makeTargetGrids(alphas = alphas)
  if(max(target_grids$pop_means$alpha1)>=randomization_probability){
    warning(paste0(
      "Largest alpha allocation is greater than randomization_probability.",
      "Inferences may be suspect.")
    )
  }
  stopifnot(
    length(weight_type)==1 && weight_type %in% c("HT", "Hajek1", "Hajek2"))


  ## tibbles not allowed
  if ( "tbl_df" %in% class(data) ) {
    data <- as.data.frame(data, stringsAsFactors = FALSE)
  }
  # if (weight_type!= "unstabilized") {stop("only unstabilized weights implemented")}

  formula <- Formula::as.Formula(formula)
  len_lhs_formula <- length(formula)[1]
  # len_rhs_formula <- length(formula)[2]
  formula_terms <- stats::terms(formula)
  formula_factors <- attr(formula_terms, "factors")
  outcome_var_name <- row.names(formula_factors)[1]
  treatment_var_name <- row.names(formula_factors)[2]
  grouping_var_name <- row.names(formula_factors)[nrow(formula_factors)]

  var_names <- list(
    outcome = outcome_var_name,
    treatment = treatment_var_name,
    grouping = grouping_var_name
    # participation = model_DV_var_name
  )
  if (len_lhs_formula==3){
    # stop("participation not implemented yet")
    # model_DV_var_name <- row.names(formula_factors)[3]
    # var_names$participation <- model_DV_var_name
    var_names$participation <- row.names(formula_factors)[3]
  } #else { it stays NULL
    # model_DV_var_name <- treatment_var_name
    # var_names$participation <- model_DV_var_name
  # }
  modeling_formula <- formula(
    stats::terms(formula, lhs = len_lhs_formula, rhs = -2)
  )

  treatment_vector <- data[[var_names$treatment]]
  # model_DV_vector <- data[[model_DV_var_name]]
  outcome_vector <- data[[var_names$outcome]]
  grouping_vector <- data[[var_names$grouping]]

  if (!is.null(dots$todos)){
    message("reorder groups in dataset")
    message("pass deriv options to geex")
  }
  unique_groups <- unique(grouping_vector)
  num_groups <- length(unique_groups)

  alphas <- sort(alphas)
  stopifnot(all(alphas>=0) & all(alphas<=1))

  ## Fit GLM model
  if (model_method == "glm") {
    trt_model_obj  <- stats::glm(
      formula = modeling_formula,
      data = data,
      family = stats::binomial()
    )
    # treatment_vector <- stats::model.response(stats::model.frame(
    #   trt_model_obj$formula, data = data))
    # stopifnot(all(treatment_vector == data[[treatment_var_name]]))


    fixefs <- trt_model_obj$coefficients
    sigma <- NULL

    ##perhaps instructive?
    ps_model_matrix <- stats::model.matrix(trt_model_obj$formula,data = data)
    treatment_param_ests <- trt_model_obj$coefficients
    prob_treated <- stats::plogis(ps_model_matrix %*% treatment_param_ests)

  } else
    if ( model_method == "glmer") {

      model_options$data <- data
      model_options$formula <- modeling_formula
      # model_options$family <- "binomial"

      trt_model_obj <- do.call(lme4::glmer, model_options)

      fixefs <- lme4::getME(trt_model_obj, "beta")
      sigma <- lme4::getME(trt_model_obj, "theta")
      ps_model_matrix <- lme4::getME(trt_model_obj, "X")
      # model_matrix_Z <- lme4::getME(trt_model_obj, "Z") ##for later

    } else {
      stop("only model_method='glm' or 'glmer' implemented")
    }

  pi_ipw_alphas <- list()
  mu_alphas_ests <- list()


  num_fixefs <- length(fixefs)
  num_alphas <- length(alphas)

  for(alp_num in 1:num_alphas) {
    alpha <- alphas[alp_num]
    # }



    pi_ipw_vec <- rep(NA, num_groups)
    summand_mat <- matrix(NA, ncol = 3, nrow = num_groups)
    denom_mat <- matrix(NA, ncol = 3, nrow = num_groups)
    # group_num <- 0
    for (gg_num in 1:num_groups){
      this_group <- unique_groups[gg_num]
      group_indices <- (grouping_vector == this_group)
      # group_num <- group_num
      treatment <- treatment_vector[group_indices]
      model_matrix <- ps_model_matrix[group_indices, , drop=FALSE]
      outcome <- outcome_vector[group_indices]
      ind_z1 <- treatment==1
      ind_z0 <- treatment==0

      pi_ipw  <- calcPiIPW(
        treatment = treatment,
        model_matrix = model_matrix,
        # logit_integrand_fun,
        fixefs = fixefs,
        sigma = sigma, ## perhaps NULL
        # participation = participation,
        integrate_alphas = integrate_alphas,
        randomization_probability = randomization_probability,
        alpha = alpha
      )
      pi_ipw_vec[gg_num] <- pi_ipw


      pi_ipw_Y <- outcome*pi_ipw

      ## adjusting for pi(A_i_notj, alpha)
      sum_y_ipw1 <- sum( pi_ipw_Y[ind_z1] )/alpha
      sum_y_ipw0 <- sum( pi_ipw_Y[ind_z0] )/(1-alpha)

      summand_mat[gg_num, ] <- c(
        sum_y_ipw1, sum_y_ipw0, (alpha*sum_y_ipw1 + (1-alpha)*sum_y_ipw0)
        )

      # stop("lurking alpha in the pi term.  How to get rid of it?")
      ## The conditional treatments are the sum times alpha

    if (weight_type == "HT") {
      denom_mat[gg_num, ] <- length(treatment)
    } else
      if (weight_type == "Hajek1") {
        ## calculate individual prop scores
        indiv_prop_scores_ipw <- rep(NA, length(treatment))

        for (jj in 1:length(indiv_prop_scores_ipw)){
          indiv_prop_scores_ipw[jj] <-
            calcPiIPW(
              treatment = treatment[jj],
              model_matrix = model_matrix[jj, , drop = FALSE],
              # logit_integrand_fun,
              fixefs = fixefs,
              sigma = sigma, ## perhaps NULL
              # participation = participation[jj],
              integrate_alphas = FALSE, ##ALPHA NOT USED HERE
              randomization_probability = randomization_probability,
              alpha = NULL
              # alpha = treatment[jj]######alpha is not used in individual prop scores
            )

        }

        sum_ipw_1 <- sum(indiv_prop_scores_ipw[ind_z1])
        sum_ipw_0 <- sum(indiv_prop_scores_ipw[ind_z0])
        denom_mat[gg_num, ] <- c(sum_ipw_1, sum_ipw_0, sum_ipw_1 + sum_ipw_0)

      } else
        if (weight_type == "Hajek2") {

          sum_ipw_1 <- sum(ind_z1)*pi_ipw/alpha
          sum_ipw_0 <- sum(ind_z0)*pi_ipw/(1-alpha)
          # sum_ipw_1 <- sum(pi_ipw[ind_z1])
          # sum_ipw_0 <- sum(pi_ipw[ind_z0])
          denom_mat[gg_num, ] <-
            c(sum_ipw_1, sum_ipw_0, alpha*sum_ipw_1 + (1-alpha)*sum_ipw_0)
        } else {
          stop("weight_type not recognized")
        }
    } ## looping over gg_num
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
  ## make a way to get out the raw CPS's

  # stop("make IPW ests")
  # cluster_propensity_scores <- calcPiIPW(
  #   ##stuff
  # )
  # IPTWs <- calcFunTVIPTW(
  #   outcome = outcome_vector,
  #   treatment = treatment_vector,
  #   prob_treated = prob_treated
  # )
  # estimated_delta <- mean(IPTWs)
  # theta_hat <- c(stats::coef(trt_model_obj), estimated_delta)


  if (verbose){message('starting variance calcs')}

  ## variance estimates

  pop_mean_alphas_list <- list()
  for (alp_num in 1:num_alphas){

    mu_alpha_vec <- mu_alphas_ests[[alp_num]]
    theta_hat <- c(mu_alpha_vec, fixefs, sigma)

    alpha <- alphas[alp_num]
    # }
    geex_args_alpha <- list(
      estFUN = eeFunTV_IPTW,
      ##geex
      data = data,
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
        randomization_probability  = randomization_probability,
        weight_type = weight_type


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

    pop_mean_alphas_list[[alp_num]] <- variance_from_geex
  }

  ## #' @param pop_means_grid "difference" or others (perhaps, later).
  ## #' @param effects_grid "difference" or others (perhaps, later).
  ## #' @param num_alphas number of allocations

#
#   tidy_args <- list(
#     target_grids = target_grids,
#     pop_mean_alphas_list = pop_mean_alphas_list,
#     # num_alphas,
#     contrast_type = contrast_type
#   )
#   saveRDS(tidy_args, file = quickLookup("test_tidyTargets_glm.Rds"))

  target_ests <- tidyTargets(
    target_grids,
    pop_mean_alphas_list,
    # num_alphas,
    contrast_type
  )

  output <- list(
    estimates = target_ests,
    prop_scores = cluster_propensity_scores,
    models = list(
      propensity_model = trt_model_obj
    ),
    geex_obj_list = pop_mean_alphas_list,
    misc = list(
      randomization_probability = randomization_probability,
      integrate_alphas = integrate_alphas,
      contrast_type = contrast_type
    ),
    formulas = list(
      full = formula,
      model = modeling_formula
    )
  )
  class(output) <- "ipw_interference"
  output
}


# options(error='recover')

#' Estimating Function for IPTW
#'
#' This function is to be passed into geex::m_estimate
#'
#' @param trt_model_obj The fitted model object (usually a glm).
#' @param outcome_var_name The name of the column in the dataframe indicating outcome of interest
#' @inheritParams estimateTV_IPTW
#' @param calcFunTVIPTW this is a function object specified in the weight_type argument
#'
#' @export
eeFunTV_IPTW <- function(
  data,
  trt_model_obj,
  num_fixefs,
  var_names,
  alpha,
  integrate_alphas ,#= integrate_allocations,
  randomization_probability ,#= randomization_probability,
  weight_type ##HT, Hajek1, Hajek2
  # outcome_var_name,
  # calcFunTVIPTW
){

  ## Put this at the end
  closureModel <- geex::grab_psiFUN(data=data,object = trt_model_obj)

  # if ( "glm" %in% class(trt_model_obj) ){
    model_matrix <- geex::grab_design_matrix(
      geex::grab_fixed_formula(trt_model_obj),
      data = data)

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

    mu1 <- theta[1]
    mu0 <- theta[2]
    mu_marg <- theta[3]
    fixefs <- theta[3+(1:num_fixefs)]
    if ("glmerMod" %in% class(trt_model_obj) ) {
      sigma <- theta[4+num_fixefs]

      # pi_ipw <- stats::integrate(
      #   fun = logit_integrand_fun,
      #   lower = -5*sigma,
      #   upper = 5*sigma,
      #   fixefs = fixefs,
      #   treatment = treatment,
      #   participation = participation, ##perhaps NULL,
      #   model_matrix = model_matrix,
      #   alpha = alpha
      # )
      #
      # pi_ipw <- calcPiIPW(
      #   # logit_integrand_fun,
      #   fixefs = fixefs,
      #   treatment = treatment,
      #   sigma = sigma, ## perhaps NULL
      #   participation = participation,
      #   model_matrix = model_matrix,
      #   integrate_alphas = integrate_alphas,
      #   randomization_probability = randomization_probability,
      #   alpha = alpha
      # )
    } else
      if ( "glm" %in% class(trt_model_obj) ) {



        #
        sigma <- NULL
        # pi_ipw <-   ## make calcPiIPW funtion that handles this switch?
        #   # stats::integrate(
        #   # fun =
        #   logit_integrand_fun(
        #     # lower = -5*sigma,
        #     # upper = 5*sigma,
        #     fixefs = fixefs,
        #     treatment = treatment,
        #     participation = participation,## NULL,
        #     model_matrix = model_matrix,
        #     alpha = alpha,
        #     # integrate_alphas = TRUE,
        #     randomization_probability = randomization_probability
        #   )
        # # linear_predictor  <- model_matrix %*% fixefs
        # prob_treated <- stats::plogis(linear_predictor)

      } else {
        stop("treatment model type not recognized")
      }

    pi_ipw <- calcPiIPW(
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
    # warning("did you get the logit_integrand fun reciprocal IPW right?")
    # # num_params  <- length(theta)
    # # num_model_params <- length(coef(trt_model_obj)) ##num_params-1
    # linear_predictor  <- model_matrix %*% theta[-num_params]
    # prob_treated <- stats::plogis(linear_predictor)

    # if (weight_type == "unstabilized") {
    ## Y*A/PS - Y*(1-A)/(1-PS)
    # IPTWs <-
    #   calcFunTVIPTW(
    #     # unstabilizedIPTW(
    #     outcome = outcome,
    #     treatment = treatment,
    #     prob_treated = prob_treated
    #   )
    # # } else {
    # #   stop("weight_type == 'unstabilized' only implemented. Not stabilized yet")
    # # }
    #
    # IPTWs - theta[num_params]
    if (weight_type == "HT") {
      pi_ipw_Y <- outcome*pi_ipw

      # c(
      #   # sum( pi_ipw_Y * () )
      # )
      sum1 <- sum( pi_ipw_Y[ind_z1] )/alpha
      sum0 <- sum( pi_ipw_Y[ind_z0] )/(1-alpha)
      out_vec <- c(
        sum1 - mu1, ## for mu(1,\alpha)
        sum0 - mu0,## for mu(0,\alpha)
        (alpha*sum1 + (1-alpha)*sum0) - mu_marg ## for mu( \alpha) marginal
      )

      # score_eqns <- apply(X, 2, function(x) sum((A - rho) * x))
      #
      # ce0 <- mean(Y * (A == 0))   / (1 - alpha)
      # ce1 <- mean(Y * (A == 1)) * IPW / (alpha)
      # #
      # c(score_eqns,
      #   ce0 - theta[p - 1],
      #   ce1 - theta[p])
      # out_vec <- c(
      #   ( mean(pi_ipw_Y*ind_z1)/ (alpha) ) - mu1, ## for mu(1,\alpha)
      #   ( mean(pi_ipw_Y*ind_z0)/ (1-alpha) ) - mu0, ## for mu(1,\alpha)
      #   ( mean(pi_ipw_Y) ) - mu_marg    ## for mu( \alpha) marginal
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

          # sum_ipw_1 <- sum(ind_z1)*pi_ipw/alpha
          # sum_ipw_0 <- sum(ind_z0)*pi_ipw/(1-alpha)
          # # sum_ipw_1 <- sum(pi_ipw[ind_z1])
          # # sum_ipw_0 <- sum(pi_ipw[ind_z0])
          # denom_mat[gg_num, ] <-
          #   c(sum_ipw_1, sum_ipw_0, alpha*sum_ipw_1 + (1-alpha)*sum_ipw_0)
          #

          # out_vec <-  c(
          #   (outcome[ind_z1] - mu1) * pi_ipw[ind_z1]/alpha,
          #   (outcome[ind_z0] - mu0) * pi_ipw[ind_z0]/(1-alpha),
          #   (outcome - mu_marg) * pi_ipw
          # )
          out_vec <-  c(
            sum(outcome[ind_z1] - mu1) * pi_ipw/alpha,
            sum(outcome[ind_z0] - mu0) * pi_ipw/(1-alpha),
            sum(outcome - mu_marg) * pi_ipw
          )
        }
    # sum_Y_Pi_IPW1 <- sum(outcome * (treatment==1)* pi_ipw )
    # sum_Y_Pi_IPW0 <- sum(outcome * (treatment==1)* pi_ipw )
    # sum_Y_Pi_IPW1 <- sum(outcome * (treatment==1)* pi_ipw )

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

