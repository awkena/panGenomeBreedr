test_that("pg_query_by_af enforces required chromosome argument", {
  # Verify the protective stop condition triggers when chrom is omitted
  expect_error(
    pg_query_by_af(start = 1, end = 100),
    "Chromosome must be specified to prevent memory overflow"
  )
})

test_that("pg_query_by_af safely handles empty genotype returns", {
  # Mock the database query to simulate a region containing no variants
  local_mocked_bindings(
    pg_query_db = function(...) data.frame()
  )

  # Verify the function outputs the correct message and halts gracefully
  expect_message(
    res <- pg_query_by_af(chrom = "Chr05", start = 1, end = 100),
    "No genotype data found for the specified region."
  )

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 0)
})

test_that("pg_query_by_af correctly fetches, calculates, and filters by thresholds", {
  # Create a dummy data frame representing the output from calc_af
  dummy_af_res <- data.frame(
    variant_id = c("v1", "v2", "v3"),
    chrom = rep("Chr05", 3),
    pos = c(10, 20, 30),
    alt_af = c(0.01, 0.50, 0.99),
    stringsAsFactors = FALSE
  )

  # Mock both internal dependencies to strictly control the data flow
  local_mocked_bindings(
    pg_query_db = function(...) data.frame(dummy_col = "raw_genotypes"),
    calc_af = function(...) dummy_af_res
  )

  # Execute the function using thresholds that isolate the middle variant
  res <- pg_query_by_af(
    min_af = 0.05,
    max_af = 0.95,
    chrom = "Chr05"
  )

  # Validate that the extreme allele frequencies were successfully filtered out
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 1)
  expect_equal(res$variant_id, "v2")
  expect_equal(res$alt_af, 0.50)
})