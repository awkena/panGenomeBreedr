test_that("query_by_impact() filters variants by impact and region", {
  # Skip check if the global test connection wasn't initialized cleanly
  skip_if_not(
    exists("con_test"),
    message = "Global test connection 'con_test' not found."
  )

  # Query by impact level globally across the active mini dataset
  result_high <- query_by_impact(con = con_test, impact_level = "HIGH")

  expect_true(is.data.frame(result_high))
  if (nrow(result_high) > 0) {
    expect_true(all(result_high$impact == "HIGH"))
    expect_true(all(
      c("variant_id", "chrom", "pos", "annotation", "impact") %in%
        colnames(result_high)
    ))
  }

  # Query by impact confirming case-insensitivity (lowercase input handling)
  result_case <- query_by_impact(con = con_test, impact_level = "high")
  expect_true(is.data.frame(result_case))
  expect_equal(nrow(result_case), nrow(result_high))

  # Query with explicit coordinate window region filtering on Chr05
  result_region <- query_by_impact(
    con = con_test,
    impact_level = "HIGH",
    chrom = "Chr05",
    start = 75104537,
    end = 75106403
  )

  expect_true(is.data.frame(result_region))
  if (nrow(result_region) > 0) {
    expect_true(all(result_region$chrom == "Chr05"))
    expect_true(all(
      result_region$pos >= 75104537 & result_region$pos <= 75106403
    ))
  }

  # Query with an invalid impact parameter value to trigger validation errors
  expect_error(
    query_by_impact(con = con_test, impact_level = "nonsense"),
    "No valid impact levels provided."
  )
})