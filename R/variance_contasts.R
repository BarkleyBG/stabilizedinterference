#' @param pop_mean_alphas_list looped geex output
#' @param target_grids output from \code{\link{makeTargetGrids}}
#' @param contrast_type "difference" or others (perhaps, later).
tidyTargets <- function(
  # target_grid,
  # pop_means_grid,
  # effects_grid,
  target_grids,
  pop_mean_alphas_list,
  # num_alphas,
  contrast_type ## will eventually be phased into makeTargetGrids
){
  pop_means_grid <- target_grids$pop_means
  effects_grid <- target_grids$effects

  stopifnot( is.na( pop_means_grid$alpha2 ) )
  # stopifnot( NROW(pop_means_grid)== 3*num_alphas )
  # stopifnot( NROW(effects_grid) == choose(num_alphas, 2)*????)
  pop_means_grid$estimate <- pop_means_grid$variance <- NA
  effects_grid$estimate <- effects_grid$variance <- NA
  # stopifnot( is.null( effects_grid$alpha2[1:(3*num_alphas)] ) )
  # target_grid$estimate <- target_grid$variance <- NA

  for (row_num in 1:(NROW(pop_means_grid)) ) {
    # mu_row <- pop_means_grid[row_num, ]
    geex_output <- pop_mean_alphas_list[[pop_means_grid$alpha1_num[row_num] ]]

    grab_geex_num <-
      getGeexNum(mu_row = pop_means_grid[row_num, ], trt_char= "trt1")

    # pop_means_mu_here <- pop_means_alphas[grab_geex_num]
    estimate <- geex_output$estimates[grab_geex_num]
    variance <- geex_output$vcov[grab_geex_num, grab_geex_num]

    pop_means_grid$estimate[row_num] <- estimate
    pop_means_grid$variance[row_num] <- variance

  }

  for (row_num in 1:(NROW(effects_grid)) ) {

    # if (effects_grid$effect[row_num]=="overall" &&
    #     effects_grid$alpha1[row_num] == 0.35  &&
    #     effects_grid$alpha2[row_num] == 0.35
    # ){
    #
    #   browser()
    #   "asdasd"
    #
    # }
    # mu_row <- pop_means_grid[row_num, ]
    # if ()
    geex_output1 <- pop_mean_alphas_list[[effects_grid$alpha1_num[row_num] ]]
    geex_output2 <- pop_mean_alphas_list[[effects_grid$alpha2_num[row_num] ]]

    grab_geex_num1 <-
      getGeexNum(mu_row = effects_grid[row_num, ], trt_char= "trt1")
    grab_geex_num2 <-
      getGeexNum(mu_row = effects_grid[row_num, ], trt_char= "trt2")
    if (effects_grid$effect[row_num] == "direct") {
      stopifnot(grab_geex_num1 %in% 1:2)
      stopifnot(grab_geex_num2 == (grab_geex_num1==1)+1)

      ## TODO: make direct effect have an alpha 2.
      # stopifnot(is.na(effects_grid$alpha2_num[row_num]))
      # stopifnot(
      if (effects_grid$alpha1_num[row_num] != effects_grid$alpha2_num[row_num]){
        browser()
      }
    } else
      if (effects_grid$effect[row_num] == "total") {
        stopifnot(grab_geex_num1 %in% 1:2)
        stopifnot(grab_geex_num2 == (grab_geex_num1==1)+1)
      } else{
        if ( grab_geex_num1!=grab_geex_num2 ){
          browser()
          effects_grid[row_num, ]
        }
        # stopifnot(grab_geex_num1==grab_geex_num2)
        # stopifnot( ### could be OE(0.4, 0.4), say.
        #   effects_grid$alpha1_num[row_num] != effects_grid$alpha2_num[row_num])
      }

    # pop_means_mu_here <- pop_means_alphas[grab_geex_num]
    estimate1 <- geex_output1$estimates[grab_geex_num1]
    estimate2 <- geex_output2$estimates[grab_geex_num2]
    vcov1 <- geex_output1$vcov#[grab_geex_num1, grab_geex_num1]
    vcov2 <- geex_output2$vcov#[grab_geex_num, grab_geex_num]

    estimate <- contrastThese(
      x1 = estimate1, x2 = estimate2, contrast_type = contrast_type)

    delta_args <- list(
      vcov1 = vcov1,
      vcov2 = vcov2,
      grab_geex_num1 = grab_geex_num1,
      grab_geex_num2 = grab_geex_num2,
      contrast_type = contrast_type
    )
    # saveRDS(delta_args, file = quickLookup("test_delta_args_glm.Rds"))


    variance <- do.call(calcDeltaMethodVariance, delta_args)
    # variance <- calcDeltaMethodVariance(
    #   vcov1,
    #   vcov2,
    #   grab_geex_num1,
    #   grab_geex_num2,
    #   contrast_type
    # )

    effects_grid$estimate[row_num] <- estimate
    effects_grid$variance[row_num] <- variance
  }


  full_grid <- rbind(pop_means_grid, effects_grid)
  full_grid$std_error <- sqrt(full_grid$variance)
  # full_grid$std_error
  full_grid$lcl <- full_grid$estimate - qnorm(0.975)*full_grid$std_error
  full_grid$ucl <- full_grid$estimate + qnorm(0.975)*full_grid$std_error

  grid_names <- names(full_grid)
  new_names <- c(grid_names[! grid_names %in% c("alpha1_num", "alpha2_num")],
                 c("alpha1_num", "alpha2_num"))
  new_order <- match(new_names,grid_names)
  # names(full_grid) <- grid_names[new_order]
  full_grid <- full_grid[,new_order]

  full_grid
}


#' Calculates the estimate and variance for causal effects
#'
#' @param vcov1 the vcov matrix for first population mean
#' @param vcov2 the vcov matrix for second population mean
#' @param grab_geex_num1 the correct from of vcov1
#' @param grab_geex_num2 the correct from of vcov2
#'
#' @return the estimated variance of the causal effect estimate
calcDeltaMethodVariance <- function(
  vcov1, vcov2, grab_geex_num1, grab_geex_num2, contrast_type
){
  dim_new <- NCOL(vcov1) + 1 ## this keeps the other 2 mu's - not necessary
  # dim_new <- NROW(vcov1) + 1
  new_Sigma_mat <- matrix(
    rep(NA, dim_new*dim_new),
    ncol = dim_new,
    nrow = dim_new
  )

  if (!isTRUE(all.equal(
    vcov1[-(1:3), -(1:3)],
    vcov1[-(1:3), -(1:3)],
    tolerance = 1e-3) )) {
    stop("why isn't this true?")
  }
  if (isTRUE(all.equal(vcov1, vcov2))) {

    if (grab_geex_num1==grab_geex_num2){
      ## overall effect
      delta_method_variance <- 0
    } else {
      ## direct effect

      contrast <- rep(0, NCOL(vcov1))

      contrast <- matrix(0, ncol = 1, nrow = dim_new-1)
      contrast_vals <- getContrastVals(contrast_type = contrast_type)
      contrast[grab_geex_num1, ] <- contrast_vals[1] ##aka 1 for diff
      contrast[grab_geex_num2, ] <- contrast_vals[2] ##aka -1 for diff

      delta_method_variance <- t(contrast) %*% vcov1 %*% contrast
    }
  } else {

  new_Sigma_mat[-dim_new, -dim_new] <- vcov1
  new_Sigma_mat[ dim_new,  dim_new] <- vcov2[grab_geex_num2, grab_geex_num2]

  new_Sigma_mat[4:(dim_new-1), dim_new] <- vcov2[-(1:3), grab_geex_num2]
  new_Sigma_mat[dim_new, 4:(dim_new-1)] <- vcov2[grab_geex_num2, -(1:3)]
  new_Sigma_mat[1:3, dim_new] <- 0
  new_Sigma_mat[dim_new, 1:3] <- 0

  contrast <- matrix(0, ncol = 1, nrow = dim_new)
  contrast_vals <- getContrastVals(contrast_type = contrast_type)
  contrast[grab_geex_num1, ] <- contrast_vals[1] ##aka 1 for diff
  # contrast[dim_new-1, ] <- contrast_vals[2] ##aka -1 for diff
  # contrast[dim_new-3 + grab_geex_num2, ] <- contrast_vals[2] ##aka -1 for diff
  contrast[dim_new, ] <- contrast_vals[2] ##aka -1 for diff

  delta_method_variance <- t(contrast) %*% new_Sigma_mat %*% contrast
  }
  delta_method_variance
}




#' Plumber to assist with geex-to-TV conversion
#'
#' @param mu_row one row of the pop_means_grid
#' @param trt_char "trt1" or "trt2"
#'
#' @return the row number of the geex output to grab
getGeexNum <- function(mu_row, trt_char) {
  stopifnot(trt_char %in% c("trt1", "trt2"))
  val <- mu_row[[trt_char]]
  if ( is.na(val) ) {
    grab_geex_num <- 3
  } else
    if (val == 0) {
      grab_geex_num <- 2
    } else if (val == 1) {
      grab_geex_num <- 1
    } else {stop()}
  grab_geex_num
}



getContrastVals <- function(contrast_type){
  # g(x_1, x_2)
  if (contrast_type=="difference"){
    out <- c(1,-1)
  } else
    if (contrast_type == "ratio"){
      stop("not implemented")
      out <- c(1/x_2, -x_1/(x_2^2)) ## not sure how this works
    } else
      if (contrast_type == "negative ratio"){
        stop("not implemented")
        out <- c(-1/x_2, x_1/(x_2^2)) ## not sure how this works
      } else
        if (contrast_type == "odds ratio") {
          stop("not implemented")
          # out <- deriv(  x1*(1-x2) / (x2 * ( 1-x1)))
        } else {stop("contrast_type not recognized")}
  out
}

contrastThese <- function(x1, x2, contrast_type){
  if (contrast_type=="difference"){
    out <- x1 - x2
  } else
    if (contrast_type == "risk ratio"){
      out <- x1/x2
    } else
      if (contrast_type == "negative risk ratio"){
        out <- (1 - x1/x2)
      } else
        if (contrast_type == "odds ratio") {
          out <-  x1*(1-x2) / (x2 * ( 1-x1))
        } else { stop("contrast_type not recognized")}
  out
}
