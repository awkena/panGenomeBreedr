test_that("pgsql_query_genotypes retrieves and unpacks specific variant IDs", {
  skip_if_not_installed("dittodb")
  skip_if_not_installed("RPostgres")

  dittodb::with_mock_path("fixtures", {
    dittodb::with_mock_db({
      con <- DBI::dbConnect(
        RPostgres::Postgres(),
        dbname = 'sorghum_pangenome_db',
        host = 'localhost',
        port = 5432,
        user = 'israeltawiahtetteh'
      )

      # IDs used during recording
      my_ids <- c("INDEL_Chr03_79037889", "SNP_Chr03_79037855")

      # Execute
      res <- pgsql_query_genotypes(con, variant_ids = my_ids)

      # Structural Check
      expect_s3_class(res, "data.frame")
      expect_true("variant_id" %in% colnames(res))

      # Unpacking Check
      expect_false("calls" %in% colnames(res))
      expect_true(ncol(res) > 1000)

      # Content Check
      expect_equal(nrow(res), length(my_ids))

      DBI::dbDisconnect(con)
    })
  })
})

test_that("pgsql_query_genotypes handles missing IDs and warnings correctly", {
  # Logic check using SQLite
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")

  # close database
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  # Create necessary table structures
  DBI::dbExecute(
    con,
    "CREATE TABLE variants (variant_id TEXT, chrom TEXT, pos INTEGER)"
  )
  DBI::dbExecute(con, "CREATE TABLE genotypes (variant_id TEXT, calls TEXT)")
  DBI::dbExecute(con, "CREATE TABLE metadata (lib TEXT, array_index INTEGER)")

  # Test empty ID warning (Updated from expect_error to expect_warning)
  expect_warning(
    res_empty_input <- pgsql_query_genotypes(con, variant_ids = character(0)),
    "The 'variant_ids' vector is empty"
  )
  expect_equal(nrow(res_empty_input), 0)

  # Test warning for IDs not in database
  expect_warning(
    res_not_found <- pgsql_query_genotypes(con, variant_ids = c("non_existent_id")),
    "No data found"
  )
  expect_equal(nrow(res_not_found), 0)
})
