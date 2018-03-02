
#' Estimate ATE with IPTW
#'
#' Fits a parametric model and estimates ATE via IPTW with Wald-type confidence
#' intervals from the empirical sandwich standard error estimates.
#'
#' @param formula Multi-part formula: \code{Outcome | Treatment ~
#'   model_predictors | cluster_ID}. This will be coerced to object of type
#'   \code{\link[Formula]{Formula}}. When using \code{model_method = "glmer"},
#'   then the random intercept term is supplied in the style preferred by
#'   \code{lme4}'s \code{merMod}: e.g., \code{Outcome | Treatment ~
#'   model_predictors + ( 1 | cluster_ID ) | cluster_ID}.
#' @param data the dataframe. Will be coerced from "tbl_df" to data.frame.
#' @param model_method \code{"glm"} for logistic, \code{"glmer"} for logistic
#'   with single random intercept. \code{"oracle"} supported; see
#'   \code{\link[inferference]{interference}} for syntax.
#' @param weight_type Estimators as presented in Liu, Hudgens, and Becker-Dreps
#'   (2016) Biometrika. Select \code{"HT"} for unstabilized weights. Select
#'   \code{"Hajek1"} or \code{"Hajek2"} for stabilized weights. Select
#'   \code{"HT_TV"} for the estimators presented in Tchetgen Tchetgen and
#'   VanderWeele (2012) SMMR and Perez-Heydrich et al. (2014) Biometrics, which
#'   in general target estimands different from those in Liu, Hudgens, and
#'   Becker-Dreps (2016) Biometrika.
#' @param ... additional args. See details.
#' @param model_options passed to \code{\link[lme4]{glmer}} or perhaps glm(in
#'   future).
#' @param alphas the range of allocations or policies from 0 to 1.
#'
#'   When \code{model_method == "oracle"} then \code{model_options} must be a
#'   list with named numeric vectors \code{fixefs} and \code{var_comp}. See
#'   \code{\link{prepareOracle}}. For mixed effect model, note that that random
#'   intercept's term in the modeling formula (e.g., \code{ ( 1 | cluster_ID )
#'   }) must be omitted from \code{formula}.
#'
#'   Arguments that can be passed through \code{...} include: \itemize{ \item
#'   \code{integrate_alphas}. Not yet supported. \item \code{verbose}. Set to
#'   \code{TRUE} for more verbose messaging. Default \code{FALSE}. \item
#'   \code{contrast_type}. Not yet supported. \item \code{keep_components}. Set
#'   to \code{TRUE} for more verbose output. Default \code{FALSE}. \item
#'   \code{target_grids}. User can supply target estimands with
#'   \code{makeTargetGrids}. \item \code{deriv_control}. User can supply the
#'   \code{deriv_control} argument to \code{\link[geex]{m_estimate}} with
#'   \code{\link[geex]{setup_deriv_control}}. }
#'
#'
#' @export
estimateTV_IPTW <- function(
  data,
  formula,
  alphas,
  weight_type  = c("HT", "Hajek1", "Hajek2", "HT_TV")[3],
  model_method  = c("glm", "glmer", "oracle")[2],
  model_options = list(nAGQ = 5, family = "binomial"),
  ...
){
  # stop("take geex object out of the loop; save space!")

  dots <- list(...)
  dots_names <- names(dots)
  randomization_probability <-
    ifelse("randomization_probability" %in% dots_names,
           dots$randomization_probability, 1)
  # compute_roots <- FALSE
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
  keep_components <-
    ifelse("keep_components" %in% dots_names,
           dots$keep_components, FALSE)
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
    length(weight_type)==1 &&
      weight_type %in% c("HT", "Hajek1", "Hajek2"))#, "HT_TV"))


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
  avg_treatment <- mean(treatment_vector)
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

  ## Fitting treatment model (or getting oracle estimates)
  trt_model_obj <- fitModel(
    data = data,
    modeling_formula = modeling_formula,
    model_options = model_options,
    model_method = model_method
  )

  model_info <- getModel(trt_model_obj = trt_model_obj, data = data)

  fixefs <- model_info$fixefs
  sigma <- model_info$sigma
  ps_model_matrix <- model_info$ps_model_matrix
  x_levels <- model_info$x_levels
  # trt_model_obj <- model_fit$trt_model_obj

  num_fixefs <- length(fixefs)
  num_alphas <- length(alphas)

  pi_ipw_alphas <- list()
  mu_alphas_ests <- list()

  ## Getting point estimates - - needs cleanup
  for(alp_num in 1:num_alphas) {
    alpha <- alphas[alp_num]


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


  # if (model_method != "oracle") {
    ## Prep for variance
    var_args <- list(
      alphas =   alphas,
      num_alphas =   num_alphas,
      mu_alphas_ests =   mu_alphas_ests,
      data =   data,
      num_fixefs =   num_fixefs,
      fixefs =   fixefs,
      sigma = sigma,
      x_levels = x_levels,
      trt_model_obj = trt_model_obj,
      var_names = var_names,
      target_grids =   target_grids,
      # pop_mean_alphas_list =   # pop_mean_alphas_list,
      randomization_probability = randomization_probability,
      weight_type = weight_type,
      verbose = verbose,

      keep_components = keep_components,
      compute_roots = compute_roots,
      # roots = c(unlist(mu_alphas_ests))
      integrate_alphas =   integrate_alphas,
      deriv_control = deriv_control,
      contrast_type =   contrast_type
    )


    var_list <- do.call(estimateVarianceCombined, var_args)
    target_ests <- var_list$target_ests
  # } else {
  #
  #     ## Prep for variance
  #     var_args <- list(
  #       alphas =   alphas,
  #       num_alphas =   num_alphas,
  #       mu_alphas_ests =   mu_alphas_ests,
  #       data =   data,
  #       num_fixefs =   num_fixefs,
  #       fixefs =   fixefs,
  #       sigma = sigma,
  #       trt_model_obj = trt_model_obj,
  #       var_names = var_names,
  #       target_grids =   target_grids,
  #       # pop_mean_alphas_list =   # pop_mean_alphas_list,
  #       randomization_probability = randomization_probability,
  #       weight_type = weight_type,
  #       verbose = verbose,
  #
  #       keep_components = keep_components,
  #       compute_roots = compute_roots,
  #       integrate_alphas =   integrate_alphas,
  #       deriv_control = deriv_control,
  #       contrast_type =   contrast_type
  #     )
  #
  #
  #     var_list <- do.call(estimateVarianceCombined, var_args)
  #     target_ests <- var_list$target_ests
  #
  # }

  # saveRDS(var_args, file = quickLookup("var_args.Rds"))

  # var_list <- do.call(estimateVarianceByAlpha, var_args)
  # # var_list <- do.call(estimateVarianceCombined, var_args)
  # target_ests <- var_list$target_ests
  # pop_mean_alphas_list <- var_list$pop_mean_alphas_list
  #
  #
  # #   target_grids,
  # #   pop_mean_alphas_list,
  # #   # num_alphas,
  # #   contrast_type
  # # )
  #
  # output <- list(
  #   estimates = target_ests,
  #   prop_scores = cluster_propensity_scores,
  #   models = list(
  #     propensity_model = trt_model_obj
  #   ),
  #   geex_obj_list = pop_mean_alphas_list,
  #   misc = list(
  #     randomization_probability = randomization_probability,
  #     integrate_alphas = integrate_alphas,
  #     contrast_type = contrast_type
  #   ),
  #   formulas = list(
  #     full = formula,
  #     model = modeling_formula
  #   )
  # )

  output <- list(
    estimates = target_ests,
    prop_scores = cluster_propensity_scores,
    models = list(
      propensity_model = trt_model_obj
    ),
    geex_obj = var_list$geex_object,
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





