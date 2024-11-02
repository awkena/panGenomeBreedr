test_that("calc_rpp_exp works", {

  # Example: Calculate RPP up to BC5
  rpp_exp <- calc_rpp_exp(bc_gen = 10, rpp2n = TRUE)

  rpp_bc10 <- calc_rpp_exp(bc_gen = 10, rpp2n = FALSE)


  expect_equal(length(rpp_exp), 11)
  expect_equal(rpp_bc10, c(bc10 = 0.9995))

})
