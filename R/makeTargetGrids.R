
#' Makes a dataframe "grid" of estimands to estimate
#'
#' This function can be called prior to and then passed into
#' \code{\link{estimateTV_IPTW}} informally by specifying \code{target_grids =
#' makeTargetGrids(alphas = c(0.4,0.5))} via the \code{\link{...}} argument.
#'
#' @inheritParams estimateTV_IPTW
#' @param treatments defaults to \code{0:1} as in TV2012
#'
#' @seealso The \code{effect_grid} function in the  inferference package
#'
#' @return a list of pop_means and effects grids
#' @export
#'
#' @author Bradley Saul and Brian G. Barkley
makeTargetGrids <- function(alphas, treatments=c(0,1)){

  ##TODO: Add a column for `contrast_type` for the effects; then multiple
  ##contrast types can be estimated with the same call.
  ##contrast_types = c("difference", "negative risk ratio") as an argument
  stopifnot(all(sort(alphas)==alphas))
  stopifnot(all(alphas>=0))
  stopifnot(all(alphas<=1))

  seq_along_alphas <- seq_along(alphas)
  # # len3 <- 3*length(alphas)
  #
  # mu_grid <- data.frame(
  #   alpha1_num = rep(1:length(alphas), each = 3),
  #   alpha2_num = rep(NA, length(alphas)*3),
  #   alpha1 = rep(alphas, each = 3),
  #   trt1 = rep(c(0,1,NA), each = length(alphas))#,
  #   # alpha2 = rep(NA, len3),
  #   # trt2 = rep(NA, len3)
  # )
  # mu_grid$alpha2 <- NA
  # mu_grid$trt2 <- NA
  #
  # ce_grid <- expand.grid(
  #   alpha1_num = mu_grid$alpha1_num,
  #   alpha2_num = mu_grid$alpha1_num
  # )

  ## remove this later
  # full_grid <- inferference:::effect_grid(alphas)


  marginal <- c("TRUE", "FALSE")
  mu_cond <- expand.grid(alpha1_num = seq_along_alphas, trt1 = treatments,
                      alpha2_num = NA, trt2 = NA,
                      marginal = FALSE, effect_type = "outcome",
                      effect = "outcome", stringsAsFactors = FALSE)
  mu_marg <- expand.grid(alpha1_num = seq_along_alphas, trt1 = NA,
                      alpha2_num = NA, trt2 = NA,
                      marginal = TRUE, effect_type = "outcome",
                      effect = "outcome", stringsAsFactors = FALSE)

  ## effect grids
  dfx <- expand.grid(alpha1_num = seq_along_alphas, trt1 = treatments,
                    alpha2_num = seq_along_alphas, trt2 = treatments,
                    marginal = FALSE, effect_type = "contrast",
                    effect = "direct", stringsAsFactors = FALSE)
  dfx <- dfx[dfx$trt1 != dfx$trt2, ]
  dfx <- dfx[dfx$alpha1_num == dfx$alpha2_num, ]
  ifx <- expand.grid(alpha1_num = seq_along_alphas, trt1 = treatments,
                    alpha2_num = seq_along_alphas, trt2 = NA,
                    marginal = FALSE, effect_type = "contrast",
                    effect = "indirect", stringsAsFactors = FALSE)
  ifx$trt2 <- ifx$trt1
  tfx <- expand.grid(alpha1_num = seq_along_alphas, trt1 = treatments,
                    alpha2_num = seq_along_alphas, trt2 = treatments,
                    marginal = FALSE,
                    effect_type = "contrast", effect = "total",
                    stringsAsFactors = FALSE)
  tfx <- tfx[tfx$trt1 != tfx$trt2, ]
  ofx <- expand.grid(alpha1_num = seq_along_alphas, trt1 = NA,
                    alpha2_num = seq_along_alphas, trt2 = NA,
                    marginal = TRUE, effect_type = "contrast",
                    effect = "overall", stringsAsFactors = FALSE)



  ## combining

  mu_grid <- rbind(mu_cond, mu_marg)
  effects_grid <- rbind(dfx, ifx, tfx, ofx)

  ## adding the alphas back on

  mu_grid$alpha1 <- alphas[mu_grid$alpha1_num]
  mu_grid$alpha2 <- alphas[mu_grid$alpha2_num]

  effects_grid$alpha1 <- alphas[effects_grid$alpha1_num]
  effects_grid$alpha2 <- alphas[effects_grid$alpha2_num]

  # if (length(contrast_types==1)){
  #   effects_grid$contrast_type <- contrast_types
  # } else {
  #   multiple_fx_grids <- lapply(contrast_types, function(contrast_type_here){
  #     effects_grid$contast_type <- contrast_type_here
  #   }
  #   effects_grid <- do.call(rbind, multiple_fx_grids)
  # }

  list_out <- list(
    pop_means = mu_grid,
    effects = effects_grid
  )
  list_out
}
