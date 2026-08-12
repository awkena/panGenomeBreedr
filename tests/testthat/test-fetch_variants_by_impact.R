test_that("fetch_variants_by_impact() works with a single impact level", {
  high_impact_variants <- fetch_variants_by_impact(
    con = con_test,
    impact_level = "HIGH",
    chrom = "Chr05",
    start = 75104537,
    end = 75106403
  )

  expect_s3_class(high_impact_variants, "data.frame")
  expect_gt(nrow(high_impact_variants), 0)
  expect_true(all(high_impact_variants$impact == "HIGH"))
})

test_that("fetch_variants_by_impact() works with multiple impact levels", {
  variants <- fetch_variants_by_impact(
    con = con_test,
    impact_level = c("HIGH", "MODERATE"),
    chrom = "Chr05" # Test without start/end
  )

  expect_s3_class(variants, "data.frame")
  expect_gt(nrow(variants), 0)
  expect_true(all(variants$impact %in% c("HIGH", "MODERATE")))
})

test_that("fetch_variants_by_impact() handles invalid inputs gracefully", {
  # Error on invalid connection
  expect_error(
    fetch_variants_by_impact(con = NULL, impact_level = "HIGH"),
    "A valid 'con' object is required for 'local' mode. Please run connect_local_db() first.",
    fixed = TRUE
  )

  # Error on invalid impact level
  expect_error(
    fetch_variants_by_impact(con = con_test, impact_level = "INVALID_IMPACT"),
    "No valid impact levels provided.",
    fixed = TRUE
  )
})
