
context("test tidyTargets()")


test_that(
  "tidyTargets working well",{

    # saveRDS(tidy_args, file = quickLookup("test_tidyTargets_glm.Rds"))
    tidy_args <- readRDS(file = quickLookup("test_tidyTargets_glm.Rds"))

    # tidy_args$pop_mean_alphas_list[[1]]@GFUN <- NULL
    target_ests <- do.call(tidyTargets, tidy_args)
    #   target_grids,
    #   pop_mean_alphas_list,
    #   # num_alphas,
    #   contrast_type
    # )

    expect_true(
      target_ests$ucl[1] != target_ests$ucl[2]
    )
    expect_true(
      target_ests$std_error[1] != target_ests$std_error[2]
    )

    # expect_failure(
      expect_true(
        abs(target_ests$std_error[30] ) < 1e-1
      )
    # )
  }
)

test_that(
  "tidyTargets working well",{

    # saveRDS(tidy_args, file = quickLookup("test_tidyTargets_glm.Rds"))
    delta_args <- readRDS(file = quickLookup("test_delta_args_glm.Rds"))


    variance <- do.call(calcDeltaMethodVariance, delta_args)

    expect_equal(
      variance,
      0,
      # tol=1e-5,
      check.attributes= FALSE
    )


  }
)
