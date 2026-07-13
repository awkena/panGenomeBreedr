# Helper function to ensure a clean state for each test, preventing tests
# from interfering with each other or the user's global environment.
with_clean_api_state <- function(code) {
  # Capture original state
  orig_opt <- getOption("panGenomeBreedr.api_url")
  orig_env <- Sys.getenv("PANGENOME_API_URL", unset = NA)

  # Guarantee cleanup after the test code runs
  on.exit({
    options(panGenomeBreedr.api_url = orig_opt)
    if (is.na(orig_env)) {
      Sys.unsetenv("PANGENOME_API_URL")
    } else {
      Sys.setenv(PANGENOME_API_URL = orig_env)
    }
  }, add = TRUE)

  # Unset state before running the test code
  options(panGenomeBreedr.api_url = NULL)
  Sys.unsetenv("PANGENOME_API_URL")

  # Execute the test code
  force(code)
}

test_that("set_api_url() correctly configures the API endpoint", {
  with_clean_api_state({
    test_url <- "http://test-server.edu:8000"

    # It should set both the R option and the system environment variable
    expect_message(
      set_api_url(test_url),
      "panGenomeBreedr API endpoint successfully set to: http://test-server.edu:8000"
    )
    expect_equal(getOption("panGenomeBreedr.api_url"), test_url)
    expect_equal(Sys.getenv("PANGENOME_API_URL"), test_url)
  })
})

test_that("set_api_url() automatically strips trailing slashes", {
  with_clean_api_state({
    test_url_with_slash <- "http://another-server.com:9000/"
    expected_stripped_url <- "http://another-server.com:9000"

    expect_message(
      set_api_url(test_url_with_slash),
      "successfully set to: http://another-server.com:9000"
    )
    expect_equal(getOption("panGenomeBreedr.api_url"), expected_stripped_url)
    expect_equal(Sys.getenv("PANGENOME_API_URL"), expected_stripped_url)
  })
})


test_that("get_api_url() retrieves the URL based on the correct hierarchy", {
  # Test Case 1: Should return the default URL when nothing is set.
  with_clean_api_state({
    expect_equal(get_api_url(), "http://79.72.72.212:8000")
  })

  # Test Case 2: Should read from the environment variable if the option is not set.
  with_clean_api_state({
    Sys.setenv(PANGENOME_API_URL = "http://env-var-url.com")
    expect_equal(get_api_url(), "http://env-var-url.com")
  })

  # Test Case 3: Should prioritize the R option over the environment variable.
  with_clean_api_state({
    options(panGenomeBreedr.api_url = "http://option-url.com")
    Sys.setenv(PANGENOME_API_URL = "http://env-var-url.com")
    expect_equal(get_api_url(), "http://option-url.com")
  })

  # Test Case 4: Should correctly retrieve the URL set by set_api_url().
  with_clean_api_state({
    set_api_url("http://integration-test.com")
    expect_equal(get_api_url(), "http://integration-test.com")
  })
})