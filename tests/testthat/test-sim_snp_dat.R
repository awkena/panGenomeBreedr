test_that("sim_snp_dat works", {

  # Simulate data for 20 snp and 100 individuals on chr 1
  geno_data <- sim_snp_dat(nsnp = 20,
                           nobs = 100,
                           chr = 1,
                           start = 1000,
                           end = 20000,
                           add_LD = TRUE,
                           LD_range = c(0.2, 1))

  expect_equal(dim(geno_data), c(100, 20))
})
