test_that("pgsql_query_by_af filtering logic works independently of DB", {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")

  # Setup dummy tables
  DBI::dbExecute(
    con,
    "CREATE TABLE genotypes (variant_id TEXT, chrom TEXT, pos INTEGER, calls TEXT)"
  )
  DBI::dbExecute(con, "CREATE TABLE metadata (lib TEXT, array_index INTEGER)")
  DBI::dbExecute(
    con,
    "CREATE TABLE variants (variant_id TEXT, variant_type TEXT, ref TEXT, alt TEXT)"
  )

  # Insert data that will NOT pass a 0.1 to 0.9 filter
  DBI::dbExecute(
    con,
    "INSERT INTO genotypes VALUES ('v1', 'Chr01', 10, '0/0'), ('v2', 'Chr01', 20, '1/1')"
  )
  DBI::dbExecute(con, "INSERT INTO metadata VALUES ('Sample1', 0)")
  DBI::dbExecute(
    con,
    "INSERT INTO variants VALUES ('v1', 'SNP', 'A', 'T'), ('v2', 'SNP', 'G', 'C')"
  )

  # Wrap the call in expect_warning to catch the expected "No variants passed" message
  expect_warning(
    {
      res_none <- pgsql_query_by_af(
        con = con,
        chrom = "Chr01",
        start = 1,
        end = 100,
        min_af = 0.1,
        max_af = 0.9
      )
    },
    "No variants passed"
  )

  # Ensure the result is still a 0-row data frame
  expect_equal(nrow(res_none), 0)

  DBI::dbDisconnect(con)
})