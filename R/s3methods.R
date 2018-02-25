


#' Prints small amount of output from a ipw_interference object
#'
#' @param x object of class "ipw_interference"
#' @param ... dots
#'
#' @method print ipw_interference
#'
#' @export
print.ipw_interference <- function(x,...){
  stopifnot(class(x) == "ipw_interference")

  ests <- x$estimates
  good_names <- c("alpha1", "trt1", "alpha2", "trt2", "estimate", "std_error")

  cat("--- estimates ---- \n")
  # cat(head(x[,good_names]))
  # cat(utils::tail(x[,good_names]))
  print(ests[1:min(7, NROW(ests)), good_names], row.names = FALSE, digits = 3)
  if (NROW(ests)>7){
    cat("     ...        ...      \n")
    ests2 <- ests[ests$effect_type != 'outcome' &
              (ests$alpha1 > ests$alpha2), good_names]
    print(ests2[utils::tail(seq.int(NROW(ests2))), ],row.names = FALSE, digits = 3)
  }
  cat("--- formula ---- \n")
  print((x$formulas$full))
  cat("--- done ---- \n")
}


#' Returns the model coefficients from ipw_interference object
#'
#' @param object object of class "ipw_interference"
#' @param ... dots
#'
#' @method coef ipw_interference
#'
#' @export
coef.ipw_interference <- function(object, ...){
  model <- object$models$propensity_model
  if ("lm" %in% class(model)){
    coeffs <- stats::coef(model)
  }
  if ("glmerMod" %in% class(model)){
    coeffs <- lme4::getME(model, c("beta", "theta"))
  }
  coeffs

}



# print.policyFX <- function(x,...){
#
#   ests <- x$estimates
#   keep_rows <- (ests$alpha1 != ests$alpha2)
#   keep_rows[is.na(keep_rows)] <- 1
#   ests <- ests[which(keep_rows==1), ]
#   tot_rows <- NROW(ests)
#
#   dots <- list(...)
#   if (!("nrows" %in% names(dots))){
#     nrows <- NROW(ests[is.na(ests$alpha2),]) + 2
#   } else{
#     nrows <- dots$nrows
#   }
#   nrows <- min(nrows, tot_rows)
#
#   nrows_more <- tot_rows - nrows
#
#   out_dfm <- ests[
#     1:nrows,
#     c("estimand", "estimate", "se", "LCI", "UCI")
#     ]
#
#
#   cat("------------- causal estimates --------------\n")
#   print(out_dfm, row.names = FALSE, digits = 3)
#   # cat('\n')
#   if (nrows_more>0){
#     cat('          ... and', nrows_more, 'more rows ...', "\n")
#   }
#   cat("---------------------------------------------\n")
# }


# #' Prints a summary of a policyFX object
# #'
# #' @param object object of class "policyFX"
# #' @param ... User may specify integer \code{nrows}.
# #'
# #' @method summary policyFX
# #'
# #' @author Brian G. Barkley
# #'
# #' @export
# summary.policyFX <- function(object, ...){
#
#
#   ests <- object$estimates
#   keep_rows <- (ests$alpha1 != ests$alpha2)
#   keep_rows[is.na(keep_rows)] <- 1
#   ests <- ests[which(keep_rows==1), ]
#   tot_rows <- NROW(ests)
#
#   dots <- list(...)
#   if (!("nrows" %in% names(dots))){
#     nrows <- NROW(ests[is.na(ests$alpha2),]) + 2
#   } else{
#     nrows <- dots$nrows
#   }
#   nrows <- min(nrows, tot_rows)
#
#   nrows_more <- tot_rows - nrows
#
#   out_dfm <- ests[
#     1:nrows,
#     c("estimand", "estimate", "se", "LCI", "UCI")
#     ]
#
#   model <- object$model
#   # ps <- as.vector(object$propensity_scores$CPS)
#   ps <- sprintf("%.3g",object$propensity_scores$CPS)
#   names(ps) <- object$propensity_scores$cluster_name
#
#   cat("------------- causal estimates --------------\n")
#   print(out_dfm, row.names = FALSE, digits = 3)
#   cat('\n')
#   if (nrows_more>0){
#     cat('          ... and', nrows_more, 'more rows ...', "\n")
#   }
#   cat('\n')
#   cat("-------------- treatment model -------------\n")
#
#   print(model)
#   cat('\n')
#   cat("------------- propensity scores -------------\n")
#
#   print(ps, quote = FALSE)
#   cat("---------------------------------------------\n")
# }
