
context("geex::m_estimate() geex_args_alpha")


test_that(
  "glmer eeFunTV_IPTW with geex::m_estimate",{

    geex_args_alpha <- readRDS(file = quickLookup("geex_args_alpha.Rds"))
    # testthat::source_test_helpers()
    # devtools::load_all()
    # if (!exists("eeFunTV_IPTW")){ devtools::load_all()}
    geex_args_alpha$deriv_control <-
      geex::setup_deriv_control(method = "simple")
    # geex_args_alpha$estFUN <- eeFunTV_IPTW
    variance_from_geex <- do.call(geex::m_estimate, geex_args_alpha)


  }
)

test_that(
  "glmer eeFunTV_IPTW with geex::m_estimate",{


    # saveRDS(geex_args_alpha, file = quickLookup("geex_args_alpha_glm_hajek2.Rds"))

    geex_args_alpha <- readRDS(file = quickLookup("geex_args_alpha_glm_hajek2.Rds"))
    # testthat::source_test_helpers()
    # devtools::load_all()
    # if (!exists("eeFunTV_IPTW")){ devtools::load_all()}
    geex_args_alpha$deriv_control <-
      geex::setup_deriv_control(method = "simple")
    geex_args_alpha$estFUN <- eeFunTV_IPTW
    variance_from_geex <- do.call(geex::m_estimate, geex_args_alpha)


  }
)
