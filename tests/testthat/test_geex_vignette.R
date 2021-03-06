
context("test_geex_vignette")




test_that(
  "I can reproduce geex vignette, with extra DE",
  {
    eefun <- function(data, model, alpha1, alpha2){
      X <- model.matrix(model, data = data)
      A <- model.response(model.frame(model, data = data))
      Y <- data$Y

      function(theta){
        p  <- length(theta)
        p1 <- length(coef(model))
        lp  <- X %*% theta[1:p1]
        rho <- plogis(lp)

        hh1  <- ((rho/alpha1)^A * ((1-rho)/(1-alpha1))^(1 - A))
        IPW1 <- 1/(exp(sum(log(hh1))))
        hh2  <- ((rho/alpha2)^A * ((1-rho)/(1-alpha2))^(1 - A))
        IPW2 <- 1/(exp(sum(log(hh2))))

        score_eqns <- apply(X, 2, function(x) sum((A - rho) * x))
        ce01 <- mean(Y * (A == 0)) * IPW1 / (1 - alpha1)
        ce11 <- mean(Y * (A == 1)) * IPW1 / (alpha1)
        ce02 <- mean(Y * (A == 0)) * IPW2 / (1 - alpha2)
        ce12 <- mean(Y * (A == 1)) * IPW2 / (alpha2)
        c(score_eqns,
          ce01 - theta[p - 5],
          ce11 - theta[p-4],
          (ce11-ce01) - theta[p-3],

          ce02 - theta[p - 2],
          ce12 - theta[p-1],
          (ce12-ce02) - theta[p]
          )
      }
    }

    deriv_method <- "simple"

    library(geex)
    library(inferference)
    my_glm_formula <-  Y | A ~ X1 | group
    data  <-  vaccinesim[1:100,]
    data$group <- rep(1:10, each = 10)


    model_method <-  'glm'
    allocations <-  c(.35, .4)
    my_seed <- 222

    set.seed(my_seed)

    # tv <- interference(
    #   formula = my_glm_formula,
    #   data   = data,
    #   model_method = model_method,
    #   allocations = allocations,
    #   method = deriv_method,
    #   runSilent = !interactive()
    # )
    set.seed(my_seed)

    mglm <- glm(A ~ X1, data = data, family = binomial)

    bs <- m_estimate(
      estFUN = eefun,
      data  = data,
      units = 'group',
      root_control = setup_root_control(start = c(
        coef(mglm), .4,  .13,.120 ,  .35,  .10,.1)),
      outer_args = list(alpha1 = allocations[1],
                        alpha2 = allocations[2], model = mglm),
      deriv_control = geex::setup_deriv_control(method = deriv_method)

    )


    set.seed(my_seed)
    glm_fit1 <- estimateEffects(
      data = data,
      formula = my_glm_formula,
      alphas = allocations[1:2],
      weight_type = "HT",
      # weight_type = c("HT", "Hajek1", "Hajek2")[3],
      # model_method = c("glm", "glmer")[1],
      model_method =  "glm" ,
      deriv_control = geex::setup_deriv_control(method = deriv_method),
          keep_components = TRUE,
      verbose = interactive(),
      model_options = NULL#list(nAGQ = 5, family = "binomial"),
      # ...
    )

    roots(bs)[c(3,6,4,7,5,8)]
    (vcov(bs) %>% diag())[c(3,6,4,7,5,8)] %>% sqrt
    # direct_effect(tv, allocation = .35)$estimate
    zz1 <- glm_fit1$estimates
    zz1[9,]$estimate
    # zz1[zz1$alpha1 ==allocations[1] & zz1$trt1 ==1 & zz1$trt2==0, ]

    # me_idx <- c(1,2,7,8)
    me_idx <- c(1,2, 4,5, 7,8)
    bs_idx <- c(4,3, 7,6,1,2)
    expect_equal(
      glm_fit1$geex_obj@estimates[me_idx],#,, me_idx],
      bs@estimates[bs_idx],#, bs_idx],
      tol=1e-8,
      check.attributes = FALSE
    )
    expect_equal(
      glm_fit1$geex_obj@vcov[me_idx, me_idx],
      bs@vcov[bs_idx, bs_idx],
      tol=1e-7
    )

    l11 <- c(0,0,1,-1,0, 0,0,0)
    expect_equal(
      t(l11) %*% vcov(bs) %*% l11,
      vcov(bs)[5,5],
      # glm_fit1$geex_obj$vcov[me_idx, me_idx],
      # bs@vcov[bs_idx, bs_idx],
      tol=1e-8,
      check.attributes = FALSE
    )

    l11 <- c(0,0,1,-1,0)
    expect_equal(
      glm_fit1$estimates$std_error[7],
      sqrt(vcov(bs)[5,5]),
      # glm_fit1$geex_obj$vcov[me_idx, me_idx],
      # bs@vcov[bs_idx, bs_idx],
      tol=1e-8,
      check.attributes = FALSE
    )

    sigma1 <- glm_fit1$geex_obj@vcov#[me_idx, me_idx]
    lll1 <- c(1,-1,0, 0,0,0, 0,0)
    myvar1 <- t(lll1) %*% sigma1 %*% lll1
    expect_equal(
      myvar1,
      vcov(bs)[5,5],
      check.attributes = FALSE,
      tol=1e-8
    )
    expect_equal(
      sqrt(myvar1),
      glm_fit1$estimates$std_error[7],
      check.attributes = FALSE,
      tol=1e-8
    )

#
#     ####### Y1
#     ### ests
#     tvy1 <- tv$estimates[1,]
#     bsy1 <- roots(bs)[3]
#     mey1 <- zz1[1,]
#
#     expect_equal( tvy1$estimate, bsy1, tol=1e-7 , check.attributes = FALSE)
#     expect_equal( tvy1$estimate, mey1$estimate, tol=1e-7 , check.attributes = FALSE)
#
#     ### Var
#
#     L1 <- c(0, 0, 1, 0)
#     Sigma <- vcov(bs)
#     bssd1 <- sqrt(t(L1) %*% Sigma %*% L1)  # from GEEX
#     tvsd1 <- tv$estimates[1,]$std.error # from inferference
#     mesd1 <- zz1[1,]$std_error
#
#     expect_equal( tvsd1,bssd1, tol=1e-2 , check.attributes = FALSE)
#     # expect_failure(
#       expect_equal(bssd1, mesd1, tol=1e-7 , check.attributes = FALSE)
#     # )
#
#     ####### DE
#
#     ## [1] 0.2667871
#     tvde <- direct_effect(tv, allocation = .35)$estimate
#     ## 0.2667872
#     bsde <- roots(bs)[3] - roots(bs)[4]
#     mede <- zz1[9, ]$estimate
#
#     expect_equal( tvde, bsde, tol=1e-7 , check.attributes = FALSE)
#     expect_equal( bsde, mede, tol=1e-7 , check.attributes = FALSE)
#
#
#     ### Var
#
#     LDE <- c(0, 0, 1, -1)
#     Sigma <- vcov(bs)
#     bssdde <- sqrt(t(LDE) %*% Sigma %*% LDE)  # from GEEX
#     tvsdde <- direct_effect(tv, allocation = .35)$std.error # from inferference
#     mesdde <- zz1[9,]$std_error
#
#     expect_equal( tvsdde,bssdde, tol=1e-1 , check.attributes = FALSE)
#     # expect_failure(
#       expect_equal(bssdde, mesdde, tol=1e-2, check.attributes = FALSE)
#     # )


      # LDE1 <- c(0, 0, -1, 1,0, 0,0,0)
      # LDE1%*%roots(bs)
      # roots(bs)[5]
      # Sigma <- vcov(bs)
      # # Sigma[4,5] <- 1
      # # Sigma[5,4] <- 1
      # mmat1 <- t(LDE1) %*% Sigma %*% LDE1
      # bssdde1 <- sqrt(mmat1)  # from GEEX
      # bssdde1
      # glm_fit1$estimates[7, c("estimate", "std_error")]
      # s1 <- glm_fit1$geex_obj$sandwich_components
      # sbs <- bs@sandwich_components


      # me_idx <- c(1,2,4,5)
      # bs_idx <- c(4,3,1,2)

      # s1@.A_i[[1]]
      # sbs@.A_i[[1]]
      # s1@.B_i[[1]][me_idx,me_idx]
      # sbs@.B_i[[1]][bs_idx, bs_idx]
      # s1@.A_i[[1]][me_idx,me_idx]
      # sbs@.A_i[[1]][bs_idx, bs_idx]
      # s1@.A[me_idx,me_idx]
      # sbs@.A[bs_idx, bs_idx]

      # solve(s1@.A)[me_idx,me_idx]
      # solve(sbs@.A)[bs_idx, bs_idx]
      # s1@.B[me_idx,me_idx]
      # sbs@.B[bs_idx, bs_idx]
      #
      #
      # (solve(s1@.A) %*% s1@.B)[me_idx,me_idx]
      # (solve(sbs@.A) %*% sbs@.B)[bs_idx, bs_idx]

      # Ss1 <- solve(s1@.A) %*% s1@.B %*% t(solve(s1@.A))
      # Ssbs <- solve(sbs@.A) %*% sbs@.B %*% t(solve(sbs@.A))

      # glm_fit1$geex_obj@vcov[me_idx, me_idx]
      # bs@vcov[bs_idx, bs_idx]


    #   var1 <- t(c(1,-1,0,0)) %*%
    #      glm_fit1$geex_obj@vcov[me_idx, me_idx] %*% c(1,-1,0,0)
    #
    #   sd1 <- sqrt(var1)
    #
    #   Ss1[me_idx,me_idx]
    #   Ssbs[bs_idx, bs_idx]
    #   LDE2 <- c(0, 0, 0,0,1)
    #
    #
    #
    #   # Sigma <- vcov(bs)
    #   # Sigma[4,5] <- 1
    #   # Sigma[5,4] <- 1
    #   mmat2 <- t(LDE2) %*% Sigma %*% LDE2
    #   bssdde2 <- sqrt(mmat2)
    # # glm_fit1$estimates[7:9,]
    #
    #
    #
    # # bsb1 <- bs@sandwich_components@.B_i[[1]]
    # # meb1 <- glm_fit1$geex_obj@sandwich_components@.B_i[[1]]
    # #
    # # bsa1 <- bs@sandwich_components@.A_i[[1]]
    # # mea1 <- glm_fit1$geex_obj@sandwich_components@.A_i[[1]]
    # #
    # # ssme1 <- t(mea1) %*% meb1 %*% mea1
    # # ssbs1 <- t(bsa1) %*% bsb1 %*% bsa1
    # #
    # # me_idx <- c(1,2,4,5)
    # # bs_idx <- c(4,3,1,2)
    # #
    # # diag(bsb1)[bs_idx]
    # # diag(meb1)[me_idx]
    # #
    # # diag(bsa1)[bs_idx]
    # # diag(mea1)[me_idx]
  }
)

test_that(
  "I can reproduce geex vignette",
  {
    eefun <- function(data, model, alpha){
      X <- model.matrix(model, data = data)
      A <- model.response(model.frame(model, data = data))
      Y <- data$Y

      function(theta){
        p  <- length(theta)
        p1 <- length(coef(model))
        lp  <- X %*% theta[1:p1]
        rho <- plogis(lp)

        hh  <- ((rho/alpha)^A * ((1-rho)/(1-alpha))^(1 - A))
        IPW <- 1/(exp(sum(log(hh))))

        score_eqns <- apply(X, 2, function(x) sum((A - rho) * x))
        ce0 <- mean(Y * (A == 0)) * IPW / (1 - alpha)
        ce1 <- mean(Y * (A == 1)) * IPW / (alpha)

        c(score_eqns,
          ce0 - theta[p - 1],
          ce1 - theta[p]
        )
      }
    }

    deriv_method <- "simple"

    library(geex)
    library(inferference)
    my_glm_formula <-  Y | A ~ X1 | group
    data  <-  vaccinesim
    data$group <- rep(1:300, each = 10)


    model_method <-  'glm'
    allocations <-  c(.35, .4)
    my_seed <- 222

    set.seed(my_seed)

    tv <- interference(
      formula = my_glm_formula,
      data   = data,
      model_method = model_method,
      allocations = allocations,
      method = deriv_method,
      runSilent = !interactive()
    )
    set.seed(my_seed)

    mglm <- glm(A ~ X1, data = data, family = binomial)

    bs <- m_estimate(
      estFUN = eefun,
      data  = data,
      units = 'group',
      root_control = setup_root_control(start = c(coef(mglm), .4,  .13)),
      outer_args = list(alpha = allocations[1], model = mglm),
      deriv_control = geex::setup_deriv_control(method = deriv_method)

    )


    set.seed(my_seed)
    glm_fit1 <- estimateEffects(
      data = data,
      formula = my_glm_formula,
      alphas = allocations,
      weight_type = "HT",
      # weight_type = c("HT", "Hajek1", "Hajek2")[3],
      # model_method = c("glm", "glmer")[1],
      model_method =  "glm" ,
      deriv_control = geex::setup_deriv_control(method = deriv_method),
      verbose = interactive(),
      model_options = NULL#list(nAGQ = 5, family = "binomial"),
      # ...
    )

    roots(bs)
    direct_effect(tv, allocation = .35)$estimate
    zz1 <- glm_fit1$estimates
    zz1[9,]$estimate
    # zz1[zz1$alpha1 ==allocations[1] & zz1$trt1 ==1 & zz1$trt2==0, ]



    ####### Y1
    ### ests
    tvy1 <- tv$estimates[1,]
    bsy1 <- roots(bs)[3]
    mey1 <- zz1[2,]

    expect_equal( tvy1$estimate, bsy1, tol=1e-7 , check.attributes = FALSE)
    expect_equal( tvy1$estimate, mey1$estimate, tol=1e-7 , check.attributes = FALSE)

    ### Var

    L1 <- c(0, 0, 1, 0)
    Sigma <- vcov(bs)
    bssd1 <- sqrt(t(L1) %*% Sigma %*% L1)  # from GEEX
    tvsd1 <- tv$estimates[1,]$std.error # from inferference
    mesd1 <- zz1[2,]$std_error

    expect_equal( tvsd1,bssd1, tol=1e-2 , check.attributes = FALSE)
    # expect_failure(
    expect_equal(bssd1, mesd1, tol=1e-7 , check.attributes = FALSE)
    # )

    ####### DE

    ## [1] 0.2667871
    tvde <- direct_effect(tv, allocation = .35)$estimate
    ## 0.2667872
    bsde <- roots(bs)[3] - roots(bs)[4]
    mede <- zz1[9, ]$estimate

    expect_equal( tvde, bsde, tol=1e-7 , check.attributes = FALSE)
    expect_equal( bsde, mede, tol=1e-7 , check.attributes = FALSE)


    ### Var

    LDE <- c(0, 0, 1, -1)
    Sigma <- vcov(bs)
    bssdde <- sqrt(t(LDE) %*% Sigma %*% LDE)  # from GEEX
    tvsdde <- direct_effect(tv, allocation = .35)$std.error # from inferference
    mesdde <- zz1[9,]$std_error

    expect_equal( tvsdde,bssdde, tol=1e-1 , check.attributes = FALSE)
    # expect_failure(
    expect_equal(bssdde, mesdde, tol=1e-7, check.attributes = FALSE)
    # )

  }
)
