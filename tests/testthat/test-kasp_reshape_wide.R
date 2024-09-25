test_that("kasp_reshape_wide works", {
  dat1 <- panGenomeBreedr::beta_carotene
  plate_wide <- kasp_reshape_wide(x = dat1,
                                  subset = 'MasterPlate',
                                  snp_id = 'SNPID',
                                  geno_call = 'Call',
                                  idvar = "SubjectID",
                                  blank = 'NTC')

  # Get Master Plate 1
  plate1 <- plate_wide$`SE-24-1088_P01_d1`
  expect_equal(dim(plate1), c(94, 5))

})
