test_that("pred_summary_plot works", {

  #skip_on_cran()
  path <- tempdir()
  setwd(path)
  dat1 <- panGenomeBreedr::beta_carotene
  my_sum <- kasp_color(x = dat1,
                          subset = 'plates',
                          sep = ':',
                          geno_call = 'Call',
                          uncallable = 'Uncallable',
                          unused = '?',
                          blank = 'NTC') |>

  pred_summary(snp_id = 'SNPID',
                 Group_id = 'Group',
                 Group_unknown = '?',
                 geno_call = 'Call',
                 rate_out = TRUE)

  expect_invisible(
  pred_summary_plot(x = my_sum$summ,
                    pdf = TRUE,
                    pred_cols = c('false' = 'red', 'true' = 'blue',
                                  'unverified' = 'orange2'),
                    alpha = 1,
                    text_size = 12,
                    width = 6,
                    height = 6,
                    angle = 45)
  )

  expect_true(length(list.files(path = ".", pattern = "\\.pdf$")) > 0)
  on.exit(unlink(path))
})
