
context('prepareOracle')

test_that("prepareOracle works", {

  data <- data.frame(resp = rbinom(20,1,0.5),
                     x1 = rnorm(20), x2= rbinom(20, 1, 0.5),
                     x3 = letters[sample.int(3, 20, replace = TRUE)])

  model_formula <- resp ~ x1*x2 + as.factor(x3)

  fixefs_true <- 1:6/10 ## 'true' params
  names(fixefs_true) <- c(
    "(Intercept)",
    "x1",
    "x2",
    "as.factor(x3) == 'b'",
    "as.factor(x3) == 'c'",
    "x1:x2"
    )
  oracle_parms_list <- list(
    fixefs = fixefs_true,
    var_comp = NA
  )

  expect_warning(
    try_me <- prepareOracle(
      data = data,
      model_options = oracle_parms_list,
      modeling_formula = model_formula
    )
  )
  expect_equal(try_me$will_it_work, FALSE)


  names(oracle_parms_list$fixefs) <- c(
    "(Intercept)",
    "x1",
    "x2",
    "as.factor(x3)2",
    "as.factor(x3)3",
    "x1:x2"
  )
  expect_warning(
    try_me <- prepareOracle(
      data = data,
      model_options = oracle_parms_list,
      modeling_formula = model_formula
    )
  )
  expect_equal(try_me$will_it_work, FALSE)


  names(oracle_parms_list$fixefs) <- c(
    "(Intercept)",
    "x1",
    "x2",
    "as.factor(x3)b",
    "as.factor(x3)c",
    "x1:x2"
  )
  expect_silent(
    try_me <- prepareOracle(
      data = data,
      model_options = oracle_parms_list,
      modeling_formula = model_formula
    )
  )
  expect_equal(try_me$will_it_work, TRUE)



})


test_that("prepOracle can avoid some errors",{
  data <- data.frame(resp = rbinom(20,1,0.5),
                     x1 = rnorm(20), x2= rbinom(20, 1, 0.5),
                     x3 = letters[sample.int(3, 20, replace = TRUE)])

  model_formula <- resp ~ x1*x2 + as.factor(x3)

  fixefs_true <- 1:6/10 ## 'true' params
  names(fixefs_true) <- c(
    "(Intercept)",
    "x1",
    "x2",
    "as.factor(x3) == 'b'",
    "as.factor(x3) == 'c'",
    "x1:x2"
  )



  oracle_parms_list <- list(
    fixefs = fixefs_true,
    var_comp = FALSE
  )
  expect_error(
    try_me <- prepareOracle(
      data = data,
      model_options = oracle_parms_list,
      modeling_formula = model_formula
    )
  )

  oracle_parms_list <- list(
    fixefs = fixefs_true,
    var_comp = TRUE
  )
  expect_error(
    try_me <- prepareOracle(
      data = data,
      model_options = oracle_parms_list,
      modeling_formula = model_formula
    )
  )


  oracle_parms_list <- list(
    fixefs = fixefs_true,
    var_comp = -3
  )
  expect_error(
    try_me <- prepareOracle(
      data = data,
      model_options = oracle_parms_list,
      modeling_formula = model_formula
    )
  )

  oracle_parms_list <- list(
    fixefs = fixefs_true
  )
  expect_error(
    try_me <- prepareOracle(
      data = data,
      model_options = oracle_parms_list,
      modeling_formula = model_formula
    )
  )

  names(fixefs_true) <- NULL
  oracle_parms_list <- list(
    fixefs = fixefs_true,
    var_comp = NA
  )
  expect_error(
    try_me <- prepareOracle(
      data = data,
      model_options = oracle_parms_list,
      modeling_formula = model_formula
    )
  )

})
