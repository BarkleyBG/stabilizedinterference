context("test_prepareOracle2.R")

test_that(
  "oracle with categorical predictors wont break",
  {


    break_args <- readRDS(file=quickLookup("breaking_hajek.Rds"))

    break_args$weight_type <- "Hajek2"
    break_args$data$gr <- rep(1:300, each = 10)

    break_args$data <- break_args$data[1:100,]
    break_args$alphas <- 4:5/10

    oracle_params <- structure(list(
      fixefs = structure(c(
        -2.39113592316612, -0.886830722442344, 0.938559317566697,
        -2.54026125893349, 0.0222751835707372, 3.0442752299863),
    .Names = c(
      "as.factor(node)3", "as.factor(node)5", "as.factor(node)6",
      "as.factor(node)9", "as.factor(node)10", "as.factor(node)11")),
    var_comp = structure(2.11198548859094, .Names = "gr.(Intercept)")),
    .Names = c("fixefs", "var_comp")
    )

    break_args2 <- break_args
    break_args2$model_method <- "oracle"


    oracle_prep <- prepareOracle(
      data = break_args2$data,
      model_options = oracle_params,
      modeling_formula = Y~ -1 + as.factor(node)
    )

    expect_true(oracle_prep$will_it_work)

    break_args2$model_options <- oracle_prep$model_options
    break_args2$formula <- infection | Y ~ -1 + as.factor(node)   | gr

    test_run2 <- do.call(estimateEffects, break_args2)

    expect_true(
      all(!is.na(test_run2$estimates$std_error ))
    )
    expect_true(
      all(!is.na(test_run2$estimates$estimate ))
    )


  }
)
