test_that("rpp_barplot works", {

  #skip_on_cran()
  path <- tempdir()
  setwd(path)

  # Observed RPP values
  rpp_df <- data.frame(sample_id = c('rp', 'dp', paste0('bc3_', 1:8)),
                       rpp = c(1, 0, round(seq(0.75, 0.97, len = 8), 2)))

  # Generate bar plot for RPP values
  expect_invisible(
  rpp_barplot(rpp_df,
              rpp_threshold = 0.85,
              text_size = 18,
              text_scale_fct = 0.1,
              alpha = 0.9,
              bar_width = 0.5,
              aspect_ratio = 0.5,
              pdf = TRUE,
              width = 8,
              height = 6))

  expect_true(length(list.files(path = ".", pattern = "\\.pdf$")) > 0)
  on.exit(unlink(path))



})
