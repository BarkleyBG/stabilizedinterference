
context('estimateVariance')

# test_that("estimateVariance with 1 alpha", {
#   var_args <- readRDS(file = quickLookup("var_args.Rds"))
#
#   set.seed(233)
#   var_args$deriv_control <- geex::setup_deriv_control(method = "Richardson")
#   # var_args$compute_roote <- TRUE
#   my_var_orig <- do.call(estimateVarianceByAlpha, var_args)
#   # saveRDS(my_var_orig, file = quickLookup("var_output1_orig.Rds"))
#   my_var_orig_orig <- readRDS(file = quickLookup("var_output1_orig.Rds"))
#   expect_equal( ##scores
#     my_var_orig,
#     my_var_orig_orig,
#     tol=1e-7
#   )
#
#   expect_equal( ##scores
#     my_var_orig$pop_mean_alphas_list[[1]]$vcov[4:5, 4:5],
#     my_var_orig_orig$pop_mean_alphas_list[[1]]$vcov[4:5, 4:5],
#     tol=1e-7
#   )
#   expect_equal( ##mus
#     my_var_orig$pop_mean_alphas_list[[1]]$vcov[1:3, 1:3],
#     my_var_orig_orig$pop_mean_alphas_list[[1]]$vcov[1:3, 1:3],
#     tol=1e-7
#   )
#
#   set.seed(233)
#   var_args$alphas <- var_args$alphas[1]
#   var_args$num_alphas <- 1
#   var_args$mu_alphas_ests <-   var_args$mu_alphas_ests[[1]]
#   my_var_combined <- do.call(estimateVarianceCombined, var_args)
#
#   expect_equal( ##scores
#     my_var_orig$pop_mean_alphas_list[[1]]$vcov[4:5, 4:5],
#     my_var_combined@vcov[4:5, 4:5],#[7:8, 7:8],
#     tol=1e-7
#   )
#
#   # expect_failure(
#   expect_equal( ##scores
#     my_var_orig$pop_mean_alphas_list[[1]]$vcov,#[4:5, 4:5],
#     my_var_combined@vcov,#[4:5, 4:5],#[7:8, 7:8],
#     tol=1e-7
#   )
#   # )
#
#   expect_equal(
#     my_var_orig$pop_mean_alphas_list[[1]]$sandwich_components@.B_i[[1]],
#     my_var_combined@sandwich_components@.B_i[[1]],
#     tol=1e-7
#   )
#
# })


test_that("estimateVariance with 2 alphas", {
  var_args <- readRDS(file = quickLookup("var_args.Rds"))

  set.seed(233)
  var_args$deriv_control <- geex::setup_deriv_control(method = "Richardson")
  # var_args$compute_roote <- TRUE
  my_var_orig <- do.call(estimateVarianceByAlpha, var_args)
  # saveRDS(my_var_ori
  #
  # set.seed(233)
  # var_args$deriv_control <- geex::setup_deriv_control(method = "Richardson")
  # my_var_orig <- do.call(estimateVarianceByAlpha, var_args)
  set.seed(233)
  # var_args$alphas <- var_args$alphas[1]
  # var_args$num_alphas <- 1
  # var_args$mu_alphas_ests <-   var_args$mu_alphas_ests[[1]]

  var_args$x_levels <- geex::grab_design_levels(var_args$trt_model_obj)
  my_var_combined <- do.call(estimateVarianceCombined, var_args)
  #
  all.equal( ##scores
    my_var_orig$pop_mean_alphas_list[[1]]$vcov[4:5, 4:5],
    my_var_combined$geex_object@vcov[7:8, 7:8],
    tol=1e-7
  )

  #
  all.equal( ## mu alpha1's
    my_var_orig$pop_mean_alphas_list[[1]]$vcov[1:3, 1:3],
     my_var_combined$geex_object@vcov[1:3, 1:3],
    tol=1e-2
  )


  all.equal( ## mu alpha2's
    my_var_orig$pop_mean_alphas_list[[2]]$vcov[1:3, 1:3],
     my_var_combined$geex_object@vcov[4:6, 4:6],
    tol=1e-2
  )



}
)

