
context("test_glmer8")


test_that(
  "I can reproduce geex vignette",
  {

    glmer8_args <- readRDS(quickLookup("glmer8_killer.Rds"))

    glmer8_args$data <- glmer8_args$data[1:100,]
    glmer8_args$verbose <- FALSE
    expect_silent(
      foo <- do.call(estimateTV_IPTW, glmer8_args)
    )
    # glmer8_args$data <- glmer8_args$data[1:100,]
    #
    glmer8_args$weight_type <- "HT"
    # glmer8_args$compute_roots <- TRUE
#
    set.seed(33)
    expect_silent(
      foo <- do.call(estimateTV_IPTW, glmer8_args)
    )
    library(inferference)
    set.seed(33)
    bar <- interference(
      data = glmer8_args$data[1:100,],
      formula = glmer8_args$formula,
      model_options = list(family = "binomial", nAGQ = 5),
      allocations = glmer8_args$alphas,
      method = "simple"
    )

    bar$estimates[c(3,1,5,4,2,6), c("alpha1", "trt1", "alpha2", "estimate", "std.error")]
    foo$estimates[1:6,  c("alpha1", "trt1", "alpha2", "estimate", "std_error")]

  }
)
