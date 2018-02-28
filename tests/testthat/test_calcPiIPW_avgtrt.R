
context("test calcPiIPW with avg trt")

test_that("calcPiIPW can be changed with non-alpha",{

  alpha <- 0.4
  ipw_args <- readRDS(file = quickLookup("calcIPW_avgtrt_args.Rds"))

  set.seed(333)
  pi_ipw_basic <- do.call(calcPiIPW, ipw_args)

  average_treatment <- ipw_args$alpha
  treatment <- ipw_args$treatment
  sum_trt <- sum(treatment)
  clust_size <- length(treatment)

  alpha_here <- alpha
  alpha_ratio <- alpha_here/average_treatment
  alpha_inv_ratio <- (1-alpha_here)/(1-average_treatment)
  pi_ratio <- (alpha_ratio)^(sum_trt)*(alpha_inv_ratio)^(clust_size - sum_trt)
  pi_ipw <- pi_ratio*pi_ipw_basic



  ipw_args$alpha <- alpha
  set.seed(333)
  pi_ipw_alpha <- do.call(calcPiIPW, ipw_args)
  # avg_trt <- ipw_args$alpha

  expect_equal(
    pi_ipw,
    pi_ipw_alpha,
    tol=1e-10
  )
})

