

context("tidyTargetsCombined")

test_that("tidyTargetsCombined works", {


  tidy_grid_args <- readRDS( file = quickLookup("test_tidyTargetsCombined_glm.Rds"))
  #
  te <- do.call(tidyTargetsCombined,tidy_grid_args)

  cte <- te[ (te$effect_type == "contrast"),]
  zdfm <- cte[(cte$alpha1==cte$alpha2) &
      (is.na(cte$trt1) | (cte$trt1==cte$trt2)) , ]

  expect_equal(
    zdfm$estimate,
    rep(0, NROW(zdfm))
  )

  expect_equal(
    zdfm$variance,
    rep(0, NROW(zdfm))
  )

  expect_equal(
    zdfm$std_error,
    rep(0, NROW(zdfm))
  )

  expect_equal(
    zdfm$lcl,
    rep(0, NROW(zdfm))
  )
  expect_equal(
    zdfm$ucl,
    rep(0, NROW(zdfm))
  )
})

