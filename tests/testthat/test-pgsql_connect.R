# Retrieve the internal unexported function cleanly into the test workspace
pgsql_connect_internal <- getFromNamespace("pgsql_connect", "panGenomeBreedr")

test_that("pgsql_connect creates a valid connection object (Mocked)", {
  # Set temporary dummy environment variables for this test scope
  withr::local_envvar(c(
    "PGSQL_HOST" = "localhost",
    "PGSQL_USER" = "postgres",
    "PGSQL_DB" = "test_db"
  ))

  # Simulate the database environment using dittodb
  dittodb::with_mock_db({
    con <- pgsql_connect_internal(password = "test_pass")

    expect_s4_class(con, "DBIConnection")

    # Clean up the mock connection
    DBI::dbDisconnect(con)
  })
})

test_that("pgsql_connect throws a helpful error on failure", {
  # Set the environment variables so it passes the pre-flight guard check
  # but fails during the actual connection attempt
  withr::local_envvar(c(
    "PGSQL_HOST" = "non-existent-link",
    "PGSQL_USER" = "postgres",
    "PGSQL_DB" = "test_db"
  ))

  # Test that the connection failure error string is handled gracefully
  expect_error(
    pgsql_connect_internal(password = "pass"),
    "Database host not found|Failed to connect to the database"
  )
})