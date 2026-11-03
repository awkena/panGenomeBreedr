test_that("kasp_call_rate() summarises call rates and returns plots", {
  path1 <- system.file("extdata", "Genotyping_141.010_01.csv",
                       package = "panGenomeBreedr", mustWork = TRUE)
  dat <- read_kasp_csv(file = path1, data_type = 'raw')$Data

  cr <- kasp_call_rate(x = dat, threshold = 0.95)

  expect_named(cr, c('marker_summary', 'plate_summary', 'marker_call_rate',
                     'no_call_reason', 'plate_call_rate'))
  expect_s3_class(cr$marker_call_rate, 'ggplot')
  expect_s3_class(cr$no_call_reason, 'ggplot')
  expect_s3_class(cr$plate_call_rate, 'ggplot')

  ms <- cr$marker_summary
  expect_equal(nrow(ms), 4)
  # NTC wells are excluded: 2 plates x 94 non-NTC wells
  expect_true(all(ms$n_wells == 188))
  expect_equal(ms$called + ms$no_call + ms$uncallable + ms$other, ms$n_wells)
  expect_equal(ms$SNPID[1], 'snpSB00800')       # lowest call rate first
  expect_equal(ms$called[ms$SNPID == 'snpSB00800'], 54)
  expect_true(all(ms$ci_lower <= ms$call_rate & ms$call_rate <= ms$ci_upper))
  expect_equal(ms$pass, ms$call_rate >= 0.95)

  expect_equal(nrow(cr$plate_summary), 8)
  expect_true(all(c('SNPID', 'MasterPlate') %in% colnames(cr$plate_summary)))
})

test_that("kasp_call_rate() works without a plate column", {
  dat <- data.frame(SNPID = rep(c('m1', 'm2'), each = 4),
                    Call = c('A:A', 'A:G', '?', 'NTC',
                             'C:C', 'Uncallable', 'Bad', 'C:T'))
  cr <- kasp_call_rate(dat, plate_id = NULL, threshold = 0.6)

  expect_null(cr$plate_summary)
  expect_null(cr$plate_call_rate)
  expect_equal(cr$marker_summary$n_wells, c(4, 3))
  expect_equal(cr$marker_summary$other, c(1, 0))
  expect_equal(cr$marker_summary$pass, c(FALSE, TRUE))
})

test_that("kasp_call_rate() validates its inputs", {
  dat <- data.frame(SNPID = 'm1', MasterPlate = 'p1', Call = 'A:A')
  expect_error(kasp_call_rate(dat, geno_call = 'Genotype'), 'not found')
  expect_error(kasp_call_rate(dat, threshold = 2), 'threshold')
  expect_error(kasp_call_rate(data.frame(SNPID = 'm1', MasterPlate = 'p1',
                                         Call = 'NTC')), 'non-NTC')
})
