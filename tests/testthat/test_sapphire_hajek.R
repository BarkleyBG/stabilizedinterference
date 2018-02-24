context("test_sapphire_hajek-r.R")

test_that(
  "sapphire-GMERT wont break",
  {


    break_args <- readRDS(file=quickLookup("breaking_hajek.Rds"))

    break_args$weight_type <- "Hajek2"
    break_args$data$gr <- rep(1:300, each = 10)

    break_args$data <- break_args$data[1:200,]
    break_args$alphas <- 4:5/10
    # expect_failure(
    zz <- do.call(estimateTV_IPTW,break_args)


    # )
  }
)
