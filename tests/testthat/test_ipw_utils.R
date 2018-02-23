
context('ipw utils')


test_that(
  "makeTargetGrids returns appropriate estimands",{

    zz <- makeTargetGrids(alphas = c(0.3333, .77,.9))
    fx <- zz$effects
    dfx <- fx[fx$effect=="direct",]

    expect_equal(
      dfx$alpha1_num,
      dfx$alpha2_num
    )

  }
)


