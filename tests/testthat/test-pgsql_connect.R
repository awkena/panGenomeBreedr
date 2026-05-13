test_that("pgsql_connect creates a valid connection object (Mocked)", {
  # simulate a database environment
  dittodb::with_mock_db({
    # Test that the function returns a connection object
    con <- pgsql_connect(password = "test_pass")

    expect_s4_class(con, "DBIConnection")

    # Clean up the mock connection
    DBI::dbDisconnect(con)
  })
})

test_that("pgsql_connect throws a helpful error on failure", {
  # We force an error by providing a non-existent host
  expect_error(
    pgsql_connect(host = "non-existent-link", password = "pass"),
    "Failed to connect to the database"
  )
})