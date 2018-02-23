
logitIntegrandFun <- function(
  raneff,
  treatment,
  # participation, ## perhaps NULL
  model_matrix,
  fixefs,
  sigma,
  # integrate_alphas,
  alpha,
  randomization_probability, ## i.e. 2/3 in Perez-Heydrich
  ...
){

  if(!is.matrix(model_matrix)){
    stop("why not a matrix?")
    model_matrix <- as.matrix(model_matrix)
  }

  # model_DV <- ifelse(is.null(participation[1]),treatment, participation)

  # if ( length(model_DV) != nrow(model_matrix) ) {
  #   stop('Length of treatment vector is not equal to number of observations in
  #        model_matrix')
  # }

  if ( length(treatment) != nrow(model_matrix) ) {
    stop('Length of treatment vector is not equal to number of observations in
         model_matrix')
  }


  ### calculations for glm - no random effect // no sigma
  if ( is.null(sigma) ) { ##Ignore random effect; it's a glm() model
    stopifnot(is.na(raneff))

    indiv_propensities <- randomization_probability *
      (stats::plogis(model_matrix %*% fixefs))

    bernoulli_pre_product <- (indiv_propensities/alpha)^treatment *
      ((1-indiv_propensities)/(1 - alpha))^(1-treatment)

    # if ( is.null(sigma) ) {  ##ignore random effect because it's a glm()
    # in this way dnorm integrates to one when integrating from -Inf to Inf
    # integrand <- exp(sum(log(bernoulli_pre_product))) * stats::dnorm(raneff, mean=0, sd = 1)
    cluster_propensity_score <- exp(sum(log(bernoulli_pre_product))) ## do not integrate over raneff

    return(cluster_propensity_score)
  }

  # # Check whether to ignore random effect
  # ignore_re <- (length(theta) == p || theta[p + 1] <= 0)

  ## Calculations ##
  ## Glmer calcs
  # } else {
  if ( nrow(model_matrix)==1 ) {
    linpred_vec <- drop(outer(model_matrix %*% fixefs, raneff, '+'))
    linpred_mat <- matrix(linpred_vec, byrow=TRUE,
                          nrow=1, ncol = length(linpred_vec))
    indiv_propensities <- randomization_probability *
      (stats::plogis(linpred_mat))
  } else {
    indiv_propensities <- randomization_probability *
      (stats::plogis(drop(outer(model_matrix %*% fixefs, raneff, '+'))))
  }
  # }

  bernoulli_pre_product <- (indiv_propensities/alpha)^treatment *
    ((1-indiv_propensities)/(1 - alpha))^(1-treatment)

  ## take Bernoulli product
  bernoulli_product <- apply(bernoulli_pre_product, 2,
                           function(x) exp(sum(log(x))))
  integrand <- bernoulli_product *
    stats::dnorm(raneff, mean=0, sd = sigma)#theta[p + 1])
  # }

  integrand
}


#
# logit_integrand <- function(b, X, A,
#                             parameters,
#                             allocation = A,
#                             randomization = 1)
# {
#   ## In the case of an intercept-only model, X needs to be converted to matrix
#   # for the warning to work
#   if(!is.matrix(X)){
#     X <- as.matrix(X)
#   }
#
#   theta <- parameters
#   p <- ncol(X)
#
#   ## Warnings ##
#   # if(p != ncol(X)){
#   #   stop('The number of fixed effect parameters is not equal to the number \n
#   #        of columns in the covariate matrix')
#   # }
#
#   if(length(A) != nrow(X)){
#     stop('Length of treatment vector is not equal to number of observations in
#          X matrix')
#   }
#
#   # Check whether to ignore random effect
#   ignore_re <- (length(theta) == p || theta[p + 1] <= 0)
#
#   ## Calculations ##
#   if(ignore_re){
#     pr.b <- randomization * (stats::plogis(X %*% theta[1:p]))
#   } else {
#     if (nrow(X)==1) {
#       linpred_vec <- drop(outer(X %*% theta[1:p], b, '+'))
#       ##drop() will return a vector if X has one row - BGB 2017-02-12
#       ##solution: create a matrix of one row from that vector - BGB 2017-02-12
#       linpred_mat <- matrix(linpred_vec, byrow=TRUE,
#                             nrow=1, ncol = length(linpred_vec))
#       pr.b <- randomization * (stats::plogis(linpred_mat))
#       ##pr.b should not throw errors in the apply() fun below. - BGB 2017-02-12
#     } else {
#       pr.b <- randomization * (stats::plogis(drop(outer(X %*% theta[1:p], b, '+'))))
#     }
#   }
#
#   hh <- (pr.b/allocation)^A * ((1-pr.b)/(1 - allocation))^(1-A)
#
#   if(ignore_re){
#     # in this way dnorm integrates to one when integrating from -Inf to Inf
#     out <- exp(sum(log(hh))) * stats::dnorm(b, mean=0, sd = 1)
#   } else {
#     hha <- apply(hh, 2, function(x) exp(sum(log(x))))
#     out <- hha * stats::dnorm(b, mean=0, sd = theta[p + 1])
#   }
#
#   return(out)
# }
