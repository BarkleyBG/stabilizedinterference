
#' Stabilized IPTW Estimators for Causal Effects Assuming Partial Interference
#'
#' This fits some of the IPW estimators introduced in Liu, Hudgens, and
#' Becker-Dreps (2016) Biometrika. These estimators estimate causal effects in
#' the presence of partial interference, with estimates of the asymptotic
#' variance from standard M-estimation theory.
#'
#' Note that these estimators estimate different causal estimands than those in
#' Tchetgen Tchetgen and VanderWeele (2012) SMMR that were applied in
#' Perez-Heydrich et al. (2014) Biometrics and implemented with
#' \code{\link[inferference]{interference}} by Saul and Hudgens (2017) JSS.
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
    message("reorder clusters in dataset")
    message("pass deriv options to geex")
  }
  unique_clusters <- unique(grouping_vector)
  num_clusters <- length(unique_clusters)

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


  ## Getting point estimates - - needs cleanup

  target_estimates <- estimateTargets(
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
  )



  mu_alphas_ests <- target_estimates$mu_alphas_ests
  pi_ipw_alphas <- target_estimates$pi_ipw_alphas
  cluster_propensity_scores <- target_estimates$cluster_propensity_scores


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
