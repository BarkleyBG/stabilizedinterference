
context('test_TV; integration')
library(inferference)

tv_names <- c("alpha1", "trt1", "alpha2", "trt2", "estimate", "std.error",
              "conf.low", "conf.high")
liu_names <- c("alpha1", "trt1", "alpha2", "trt2", "estimate", "std_error",
               "lcl", "ucl")


# test_that("integration Hajek2 works",{
#
#   data <- inferference::vaccinesim
#   data <- data[1:150,]
#   #
#   # data,
#   # formula,
#   # alphas,
#   # weight_type = c("HT", "Hajek1", "Hajek2")[3],
#   # model_method = c("glm", "glmer")[2],
#   # model_options = list(nAGQ = 5, family = "binomial"),
#   alphas <- 4:5/10
#   my_formula <- Y | A ~ X1*X2 + (1|group) | group
#
#   glmer_fit <- estimateTV_IPTW(
#     data = data,
#     formula = my_formula,
#     alphas = alphas,
#
#     verbose = interactive(),
#     deriv_control = geex::setup_deriv_control(method = "simple"),
#     # weight_type = "HT"
#     weight_type = c("HT", "Hajek1", "Hajek2")[3]
#     # model_method = c("glm", "glmer")[2],
#     # model_options = list(nAGQ = 5, family = "binomial"),
#     # ...
#   )
#
#   # debug(
#   #
#   glm_fit <- estimateTV_IPTW(
#     data = data,
#     formula = Y | A ~ X1*X2  | group,
#     alphas = alphas,
#     # weight_type = "HT",
#     weight_type = c("HT", "Hajek1", "Hajek2")[3],
#     # model_method = c("glm", "glmer")[1],
#     model_method =  "glm" ,
#     deriv_control = geex::setup_deriv_control(method = "simple"),
#     verbose = interactive(),
#     model_options = NULL#list(nAGQ = 5, family = "binomial"),
#     # ...
#   )
#   # )
#
# }
# )


test_that("integration HT-glm works",{


  data <- inferference::vaccinesim
  data <- data[1:1200,]
  gsize <- 5
  data$group <- rep(seq.int(NROW(data)/gsize), each =gsize)
  # data$group <- rep(1:5, each = 10)
  alphas <- c(3,6)/10


  my_glm_formula <- Y | A ~ X1*X2  | group
  # debug(
  my_seed <- 213

  set.seed(my_seed)
  glm_fit1 <- estimateTV_IPTW(
    data = data,
    formula = my_glm_formula,
    alphas = alphas,
    weight_type = "HT",
    # weight_type = c("HT", "Hajek1", "Hajek2")[3],
    # model_method = c("glm", "glmer")[1],
    model_method =  "glm" ,
    deriv_control = geex::setup_deriv_control(method = "simple"),
    # verbose = TRUE,
    model_options = NULL#list(nAGQ = 5, family = "binomial"),
    # ...
  )
  # )
  #

  set.seed(my_seed)
  glm_fit2 <- interference(
    data = data,
    formula = my_glm_formula,
    model_method  = "glm",
    method = "simple",
    # integrand = inferference::logit_integrand
    allocations = alphas
  )

  zz1 <-   glm_fit1$estimates
  zz1 <- zz1[,liu_names]
  zz2 <- glm_fit2$estimates
  zz2 <- zz2[, tv_names]
  colnames(zz2) <- liu_names
  head(zz1)
  head(zz2)
  tail(zz1)
  tail(zz2)
  expect_equal(
    glm_fit1$estimates$estimate,
    glm_fit2$estimates$estimate,
    tol = 1e-7
  )
  # expect_failure(
  expect_equal(
    glm_fit1$estimates$std_error,
    glm_fit2$estimates$std.error,
    tol = 1e-1 ### different methods
  )

  adj_se <- (glm_fit2$estimates$std.error - glm_fit1$estimates$std_error) /
    glm_fit1$estimates$std_error

  adj_se <- na.omit(adj_se)
  expect_true( all( abs(adj_se)< .25) )
  # )
}
)


# test_that("integration HT-glmer works",{
#
#   data <- inferference::vaccinesim
#   data <- data[1:200,]
#   data$group <- rep(1:10, each = 20)
#
#   alphas <- 4:5/10
#
#
#   my_glmer_formula <- Y | A ~ X1*X2 + (1|group) | group
#
#   glmer_fit1 <- estimateTV_IPTW(
#     data = data,
#     formula = my_glmer_formula,
#     alphas = alphas,
#
#     verbose = TRUE,
#     deriv_control = geex::setup_deriv_control(method = "simple"),
#     weight_type = "HT"
#     # weight_type = c("HT", "Hajek1", "Hajek2")[3],
#     # model_method = c("glm", "glmer")[2],
#     # model_options = list(nAGQ = 5, family = "binomial"),
#     # ...
#   )
#   # )
#   #
#   library(inferference)
#   glmer_fit2 <- interference(
#     data = data,
#     formula = my_glmer_formula,
#     # model_method  = "glm",
#     method = "simple",
#     model_options = list(nAGQ = 5, family = "binomial"),
#     # integrand = inferference::logit_integrand
#     allocations = alphas
#   )
#
#   expect_equal(
#     glmer_fit1$estimates$estimate,
#     glmer_fit2$estimates$estimate
#   )
# }
# )


