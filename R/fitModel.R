
#' Estimate the treatment model parameters.
#'
#' Output of this function will be handled by \code{\link{getModel}}
#'
#' @inheritParams estimateEffects
#' @param modeling_formula This is just the 'inner' part of the multi-part
#'   formula for estimating treatment on predictors.
#'
#'  When \code{model_method = "oracle"} see \code{\link{prepareOracle}}.
#'
#'  @return \code{trt_model_obj}.
#'
#'  @export
fitModel <- function(
  data,
  modeling_formula,
  model_options,
  model_method
){
  ## Fit GLM model
  if (model_method == "glm") {
    trt_model_obj  <- stats::glm(
      formula = modeling_formula,
      data = data,
      family = stats::binomial()
    )
  } else
    if (model_method == "glmer") {

      model_options$data <- data
      model_options$formula <- modeling_formula

      trt_model_obj <- do.call(lme4::glmer, model_options)

    } else
      if (model_method == "oracle") {


        prep_oracle <- prepareOracle(
          data = data,
          model_options = model_options,
          modeling_formula = modeling_formula
        )

        trt_model_obj <- prep_oracle

      } else {
        stop("only model_method='glm', 'glmer', 'oracle' supported.")
      }

  trt_model_obj
}

#' Grabs necessary information from a treatment model.
#'
#' @inheritParams estimateEffects
#' @param trt_model_obj Fitted model object from \code{\link{fitModel}}
#'
#'  When \code{model_method = "oracle"} see \code{\link{prepareOracle}}.
#'
#'  @return A list including \itemize{
#'  \item \code{fixefs}
#'  \item \code{sigma}
#'  \item \code{ps_model_matrix}
#'  \item \code{x_levels}
#'  }
getModel <- function(
  trt_model_obj, data
){
  if ("lm" %in% class(trt_model_obj)) {
    fixefs <- trt_model_obj$coefficients
    sigma <- NULL
    ps_model_matrix <- stats::model.matrix(trt_model_obj$formula,data = data)

    x_levels <- geex::grab(from = trt_model_obj, "design_levels")

  } else
     if ("glmerMod" %in% class(trt_model_obj)) {


      fixefs <- lme4::getME(trt_model_obj, "beta")
      sigma <- lme4::getME(trt_model_obj, "theta")
      ps_model_matrix <- lme4::getME(trt_model_obj, "X")

      x_levels <- geex::grab(from = trt_model_obj, "design_levels")

    } else
      # if (model_method == "oracle") {
      if (class(trt_model_obj)=="list") {
        stopifnot(
          all(
            c("fixefs", "sigma", "model_matrix") %in% names(trt_model_obj)
          )
        )
        fixefs <- trt_model_obj$fixefs
        sigma <- trt_model_obj$sigma
        ps_model_matrix <- trt_model_obj$model_matrix

        x_levels <- trt_model_obj$x_levels


      } else {
        stop("only model_method='glm', 'glmer', 'oracle' supported.")
      }
  out_list <- list(
    fixefs = fixefs,
    sigma = sigma,
    ps_model_matrix = ps_model_matrix,
    x_levels = x_levels
  )
  out_list
}
