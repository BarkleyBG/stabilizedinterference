
#' Helper function for \code{\link{tidyTargetsCombined}}
#'
#' @param trt either 1, 0, or NA. The value of \eqn{z} in \eqn{\mu(z,\alpha)}.
trtIdx <- function(trt){
  if (is.na(trt)){return(3)};   if (trt==0){return(2)};   if (trt==1){return(1)}
}

#' Take geex'd variance output and turn it into a tidy grid of target estimates
#'
#' @inheritParams getContrastVals
#' @param pop_mean_alphas_list looped geex output
#' @param target_grids output from \code{\link{makeTargetGrids}}
#' @param alpha_grid estimates and variance estimates for the pop means
#' @param geex_output from \code{\link[geex]{m_estimate}}.
#'
tidyTargetsCombined <- function(
  target_grids,
  alpha_grid,
  geex_output,
  contrast_type ## will eventually be phased into makeTargetGrids
){
  pop_means_grid <- alpha_grid
  # pop_means_grid <- target_grids$pop_means
  effects_grid <- target_grids$effects

  stopifnot( is.na( pop_means_grid$alpha2 ) )
  effects_grid$estimate <- NA
  effects_grid$variance <- NA

  Sigma_mat <- geex_output@vcov
  for (row_num in 1:(NROW(effects_grid)) ) {

    alpha1_num <- effects_grid[row_num,]$alpha1_num
    alpha2_num <- effects_grid[row_num,]$alpha2_num
    trt1 <- effects_grid[row_num,]$trt1
    trt2 <- effects_grid[row_num,]$trt2

    trt1_idx <- trtIdx(trt1)
    trt2_idx <- trtIdx(trt2)
    alpha1_idx <- (alpha1_num-1)*3 + 1:3
    alpha2_idx <- (alpha2_num-1)*3 + 1:3
    idx1 <- alpha1_idx[trt1_idx]
    idx2 <- alpha2_idx[trt2_idx]

    fx_estimate <- contrastThese(
      x1 = pop_means_grid$estimate[idx1],
      x2 = pop_means_grid$estimate[idx2],
      contrast_type = contrast_type
    )

    if (idx1 != idx2){
    contrast_mat <- matrix(rep(0, NCOL(Sigma_mat)), nrow = 1)
    contrast_vals <- getContrastVals(contrast_type = contrast_type)
    contrast_mat[ , idx1] <- contrast_vals[1]
    contrast_mat[ , idx2] <- contrast_vals[2]

    fx_variance <- contrast_mat %*% Sigma_mat %*% t(contrast_mat)
} else {fx_variance <- 0}
    effects_grid[row_num, ]$estimate <- fx_estimate
    effects_grid[row_num, ]$variance <- fx_variance

    # effects_grid$estimate[row_num] <- estimate
    # effects_grid$variance[row_num] <- variance
  }

  pop_means_grid$marginal <- is.na(pop_means_grid$trt1)
  pop_means_grid$effect_type <- "outcome"
  pop_means_grid$effect <- "outcome"

  pop_means_grid$std_error <- NA
  effects_grid$std_error <- NA
  pop_means_grid$lcl <- NA
  effects_grid$lcl <- NA
  pop_means_grid$ucl <- NA
  effects_grid$ucl <- NA

  names1 <- colnames(pop_means_grid)
  names2 <- colnames(effects_grid)
  names3 <- c(
    "alpha1", "trt1", "alpha2", "trt2",
    "estimate", "variance", "std_error", "lcl", "ucl",
    "marginal", "effect_type", "effect", "alpha1_num", "alpha2_num"
  )
  stopifnot(all(names1 %in% names2))
  stopifnot(all(names2 %in% names1))
  stopifnot(all(names3 %in% names1))
  stopifnot(all(names1 %in% names3))
  reorder1 <- match(names3, names1)
  reorder2 <- match(names3, names2)

  pop_means_grid <- pop_means_grid[ , reorder1]
  effects_grid <- effects_grid[ , reorder2]
  all(colnames(pop_means_grid)==colnames(effects_grid))


  full_grid <- rbind(pop_means_grid, effects_grid)
  full_grid$std_error <- sqrt(full_grid$variance)
  # full_grid$std_error
  full_grid$lcl <- full_grid$estimate - stats::qnorm(0.975)*full_grid$std_error
  full_grid$ucl <- full_grid$estimate + stats::qnorm(0.975)*full_grid$std_error

  # grid_names <- names(full_grid)
  # new_names <- c(grid_names[! grid_names %in% c("alpha1_num", "alpha2_num")],
  #                c("alpha1_num", "alpha2_num"))
  # new_order <- match(new_names,grid_names)
  # # names(full_grid) <- grid_names[new_order]
  # full_grid <- full_grid[,new_order]

  # class(full_grid) <- "ipw_interference_ests"
  full_grid
}

