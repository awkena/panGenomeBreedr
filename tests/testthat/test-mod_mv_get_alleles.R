
# Test the UI function
test_that("mod_mv_get_alleles_ui creates expected UI elements", {
  # Create UI with test ID
  ui <- mod_mv_get_alleles_ui("test")

  # Test that function returns a tagList
  expect_s3_class(ui, "shiny.tag.list")

  # Convert to HTML to check structure
  html <- as.character(ui)

  # Test for presence of key UI elements
  expect_true(grepl('id="test-data_type_id"', html))
  expect_true(grepl('id="test-sep_id"', html))
  expect_true(grepl('id="test-plates_pres"', html))
  expect_true(grepl('id="test-alleles_output"', html))
  expect_true(grepl('id="test-genotypes_output"', html))

  # Test selectInput choices
  expect_true(grepl('kasp', html))
  expect_true(grepl('agriplex', html))

  # Test default separator value
  expect_true(grepl('value=":"', html))

  # Test labels
  expect_true(grepl('Select Data Source', html))
  expect_true(grepl('Genotype Call Separator', html))
  expect_true(grepl('Choose plate', html))
  expect_true(grepl('Identified Alleles', html))
  expect_true(grepl('Genotype Classifications', html))
})

test_that("mod_mv_get_alleles_ui handles different IDs correctly", {
  # Test with different IDs
  ui1 <- mod_mv_get_alleles_ui("module1")
  ui2 <- mod_mv_get_alleles_ui("module2")

  html1 <- as.character(ui1)
  html2 <- as.character(ui2)

  # Check that IDs are properly namespaced
  expect_true(grepl('id="module1-data_type_id"', html1))
  expect_true(grepl('id="module2-data_type_id"', html2))
  expect_false(grepl('id="module2-data_type_id"', html1))
  expect_false(grepl('id="module1-data_type_id"', html2))
})

# Test the server function structure
test_that("mod_mv_get_alleles_server function has correct structure", {
  # Check that function exists and has correct parameters
  expect_true(exists("mod_mv_get_alleles_server"))

  # Get function arguments
  args <- formals(mod_mv_get_alleles_server)
  expect_named(args, c("id", "kasp_data"))
})

# Test server function with mock data (basic structure test)
test_that("mod_mv_get_alleles_server initializes without error", {
  # Create mock reactive data
  mock_kasp_data <- reactive({
    data.frame(
      plate = c("plate1", "plate1", "plate2"),
      sample = c("A1", "A2", "B1"),
      call = c("A:T", "G:C", "A:A")
    )
  })

  # Mock the helper functions that would be imported
  mock_uniq_plates <- function(data) {
    c("plate1", "plate2")
  }

  mock_get_calls <- function(x, a) {
    data.frame(
      sample = c("A1", "A2"),
      call = c("A:T", "G:C")
    )
  }

  mock_get_alleles <- function(x, data_type, sep) {
    list(
      alleles = data.frame(allele1 = c("A", "G"), allele2 = c("T", "C")),
      genotypes = data.frame(sample = c("A1", "A2"), genotype = c("AT", "GC"))
    )
  }

  mock_alleles_df <- function(x) x
  mock_genotypes <- function(x) x

  # Temporarily assign mock functions to global environment for testing
  assign("uniq_plates", mock_uniq_plates, envir = .GlobalEnv)
  assign("get_calls", mock_get_calls, envir = .GlobalEnv)
  assign("get_alleles", mock_get_alleles, envir = .GlobalEnv)
  assign("alleles_df", mock_alleles_df, envir = .GlobalEnv)
  assign("genotypes", mock_genotypes, envir = .GlobalEnv)

  # Test that server function can be called without immediate error
  expect_no_error({
    testServer(mod_mv_get_alleles_server, args = list(kasp_data = mock_kasp_data), {
      # Basic test that server initializes
      expect_true(TRUE)
    })
  })

  # Clean up mock functions
  rm("uniq_plates", "get_calls", "get_alleles", "alleles_df", "genotypes",
     envir = .GlobalEnv)
})

# Test server reactivity with testServer
test_that("mod_mv_get_alleles_server reactive behavior", {
  skip_if_not_installed("shinytest2")

  # Mock helper functions
  assign("uniq_plates", function(data) c("plate1", "plate2"), envir = .GlobalEnv)
  assign("get_calls", function(x, a) data.frame(sample = "A1", call = "A:T"), envir = .GlobalEnv)
  assign("get_alleles", function(x, data_type, sep) {
    list(
      alleles = data.frame(allele1 = "A", allele2 = "T"),
      genotypes = data.frame(sample = "A1", genotype = "AT")
    )
  }, envir = .GlobalEnv)
  assign("alleles_df", function(x) x, envir = .GlobalEnv)
  assign("genotypes", function(x) x, envir = .GlobalEnv)

  mock_kasp_data <- reactive({
    data.frame(plate = "plate1", sample = "A1", call = "A:T")
  })

  testServer(mod_mv_get_alleles_server, args = list(kasp_data = mock_kasp_data), {
    # Test that distinct_plates reactive works
    expect_equal(distinct_plates(), c("plate1", "plate2"))

    # Simulate user selecting a plate
    session$setInputs(plates_pres = "plate1")
    session$setInputs(data_type_id = "kasp")
    session$setInputs(sep_id = ":")

    # Test that Get_calls reactive responds to input
    expect_true(is.data.frame(Get_calls()))

    # Test that obtain_alleles reactive works
    expect_true(is.list(obtain_alleles()))
    expect_named(obtain_alleles(), c("alleles", "genotypes"))
  })

  # Clean up
  rm("uniq_plates", "get_calls", "get_alleles", "alleles_df", "genotypes",
     envir = .GlobalEnv)
})

# Test edge cases
test_that("mod_mv_get_alleles_ui handles edge cases", {
  # Test with empty ID
  expect_no_error(mod_mv_get_alleles_ui(""))

  # Test with special characters in ID
  ui_special <- mod_mv_get_alleles_ui("test-module_1")
  html_special <- as.character(ui_special)
  expect_true(grepl('id="test-module_1-data_type_id"', html_special))
})

test_that("UI component structure is correct", {
  ui <- mod_mv_get_alleles_ui("test")

  # Test that it contains sidebar layout components
  html <- as.character(ui)
  expect_true(grepl('col-sm-4', html) || grepl('sidebar', html))

  # Test for accordion components
  expect_true(grepl('accordion', html))

  # Test for DT output components by checking the actual output IDs
  expect_true(grepl('test-alleles_output', html))
  expect_true(grepl('test-genotypes_output', html))
})
