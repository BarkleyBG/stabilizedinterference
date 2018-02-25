
calcPiIPW <- function(
  # integrand_fun,# = logitIntegrandFun,
  treatment,# = treatment[jj],
  # participation,# = participation[jj], ##perhaps NULL
  model_matrix ,#= model_matrix[jj, , drop = FALSE],
  fixefs ,#= fixefs,
  sigma ,#= sigma, ## perhaps NULL
  alpha ,#= 1######alpha is not used in individual prop scores

  # parameters,
  # integrand,
  # allocation,
  integrate_alphas = TRUE,
  randomization_probability, #as in Caro
  ...
){
# stop('fix model matrix')
  stopifnot((is.null(alpha) && integrate_alphas)|| length(alpha)==1)
  ## Necessary pieces ##
  # integrand         <- match.fun(integrand_fun)
  # integrand.formals <- names(formals(integrand))
  dots              <- list(...)
  # dot.names         <- names(dots)
  # A                 <- dots[[match.arg('A', dot.names)]]
  stopifnot(integrate_alphas)

  ipw_args <- list(
    treatment = treatment,
    # participation = participation, ##will this be a problem?
    model_matrix = model_matrix,
    fixefs = fixefs,
    sigma = sigma,
    # integrate_alphas = integrate_alphas,
    randomization_probability = randomization_probability#,
    # alpha = ifelse(integrate_alphas,alpha,treatment)
  )

  if (integrate_alphas) { ipw_args$alpha <- alpha } else {ipw_args$alpha <- treatment}

  if ( is.null(sigma) ) {
    ## GLM
    ipw_args$raneff <- NA
    down_weight <- do.call(logitIntegrandFun, ipw_args)

  } else {
    ## glmer

    ipw_args$lower <- -5*sigma ##bounds for raneff or "b" random intercept
    ipw_args$upper <- 5*sigma

    ipw_args$f <- logitIntegrandFun


    logit_integral <- try(do.call(stats::integrate, args = ipw_args), silent = TRUE)

    if ( 'try-error' %in% class(logit_integral) ) {
      down_weight <- NA
    } else {
      down_weight <- logit_integral$value
    }


  }

  ## Compute the weight ##
  # weight <- 1/PrA

  if(integrate_alphas){
    up_weight <- 1/down_weight
  } else {
    alpha_prod    <- prod( alpha^treatment * (1-alpha)^(1-treatment) )
    up_weight <- alpha_prod/down_weight
  }

  up_weight ## up_weight = pi_ipw  = pi_term / cluster_prop_score
}


# ## Warnings ##
# # if(!'A' %in% dot.names){
# #   stop("The argument 'A' (treatment assignment) must be specified")
# # }
# #
# # ## Integrate() arguments ##
# # if(!'lower' %in% dot.names){
# #   dots$lower <- -Inf
# # }
# #
# # if(!'upper' %in% dot.names){
# #   dots$upper <- Inf
# # }
#
# int.args <- append(get_args(stats::integrate, dots),
#                    list(f = integrand, parameters = parameters))
#
# args <- append(get_args(integrand, dots), int.args)
#
# # alphas is optional in user-defined integrands. Include this argument
# # when necessary. Note that alphas will either be used in this function
# # or passed to the integrand function.
# if("alphas" %in% integrand.formals){
#   args$alphas <- alphas
# }
#
# ## Compute the integral ##
# # if any of the products within the integrand return Inf, then return NA
# # else return the result of integration

