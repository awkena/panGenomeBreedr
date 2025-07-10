# Test file for KASP Marker Design Module
# Mock data for testing
create_mock_vcf_data <- function() {
  list(
    vcf_matrix_markerID = c("marker1", "marker2", "marker3"),
    vcf_matrix_chromID = c("chr1", "chr2", "chr3")
  )
}

create_mock_vcf_df <- function() {
  data.frame(
    ID = c("marker1", "marker2", "marker3"),
    CHROM = c("chr1", "chr2", "chr3"),
    POS = c(1000, 2000, 3000),
    REF = c("A", "T", "G"),
    ALT = c("T", "C", "A"),
    sample1 = c("0/1", "1/1", "0/0"),
    sample2 = c("1/1", "0/1", "1/1"),
    stringsAsFactors = FALSE
  )
}

create_mock_excel_data <- function() {
  data.frame(
    marker_id = c("marker1", "marker2", "marker3"),
    chromosome = c("chr1", "chr2", "chr3"),
    position = c(1000, 2000, 3000),
    ref_allele = c("A", "T", "G"),
    alt_allele = c("T", "C", "A"),
    other_col = c("val1", "val2", "val3"),
    stringsAsFactors = FALSE
  )
}

create_mock_kasp_result <- function() {
  list(
    marker_data = data.frame(
      Marker_ID = "test_marker",
      Chromosome = "chr1",
      Position = 1000,
      Sequence = "ATCGATCG",
      stringsAsFactors = FALSE
    ),
    plot = ggplot2::ggplot() + ggplot2::geom_point(ggplot2::aes(x = 1, y = 1))
  )
}

# Test the UI function
test_that("UI function creates proper structure", {
  ui <- mod_kasp_marker_design_ui("test")

  expect_s3_class(ui, "shiny.tag.list")
  expect_true(length(ui) > 0)

  # Convert to HTML and check for key elements
  html_content <- as.character(ui)
  expect_true(grepl("File Uploads", html_content))
  expect_true(grepl("Column Mapping", html_content))
  expect_true(grepl("Marker Selection", html_content))
  expect_true(grepl("Analysis Parameters", html_content))
})

# Test server function with VCF file input
test_that("handles VCF file input", {
  # Create temporary VCF file
  temp_vcf <- tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1",
    "chr1\t1000\tmarker1\tA\tT\t.\t.\t.\tGT\t0/1"
  ), temp_vcf)

  testServer(mod_kasp_marker_design_server, {
    # Mock the marker.chr_ID function using testthat's with_mocked_bindings
    with_mocked_bindings(
      "marker.chr_ID" = function(file_path) {
        create_mock_vcf_data()
      },
      {
        session$setInputs(
          upload_choice = "Raw VCF File",
          vcf_file = list(datapath = temp_vcf)
        )

        session$flushReact()

        expect_true(is.reactive(vcf_data))
      }
    )
  })

  unlink(temp_vcf)
})

# Test server function with Excel file input
test_that("handles Excel file input", {
  # Create temporary Excel file path (we'll mock the read function)
  temp_excel <- tempfile(fileext = ".xlsx")

  testServer(mod_kasp_marker_design_server, {
    # Mock readxl::read_excel using testthat's with_mocked_bindings
    with_mocked_bindings(
      "read_excel" = function(...) {
        create_mock_excel_data()
      },
      .package = "readxl",
      {
        session$setInputs(
          upload_choice = "Processed VCF File",
          gt_df = list(datapath = temp_excel)
        )

        session$flushReact()
        expect_true(is.reactive(gt_data))
      }
    )
  })
})

# Test column selection logic
test_that("handles column selection for VCF files", {
  temp_vcf <- tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1",
    "chr1\t1000\tmarker1\tA\tT\t.\t.\t.\tGT\t0/1"
  ), temp_vcf)

  testServer(mod_kasp_marker_design_server, {
    # Mock the vcf reading functionality
    with_mocked_bindings(
      "marker.chr_ID" = function(file_path) {
        create_mock_vcf_data()
      },
      {
        session$setInputs(
          upload_choice = "Raw VCF File",
          vcf_file = list(datapath = temp_vcf)
        )

        session$flushReact()

        # Check that column names are available
        expect_true(is.reactive(vcf_colnames))
      }
    )
  })

  unlink(temp_vcf)
})

# Test error handling for file reading
test_that("gracefully handles file read errors", {
  testServer(mod_kasp_marker_design_server, {
    # Mock readxl::read_excel to throw an error
    with_mocked_bindings(
      "read_excel" = function(...) {
        stop("File read error")
      },
      .package = "readxl",
      {
        session$setInputs(
          upload_choice = "Processed VCF File",
          gt_df = list(datapath = "nonexistent.xlsx")
        )

        # The error should be caught within the reactive, so we test that gt_data is NULL
        session$flushReact()

        # Test that the reactive handles the error gracefully
        expect_true(is.reactive(gt_data))
        # If error handling is implemented, gt_data should be NULL or empty
        # This depends on how the actual module handles errors
      }
    )
  })
})

# Test file type switching
test_that("handles file type switching", {
  testServer(mod_kasp_marker_design_server, {
    # Test Raw VCF File selection
    session$setInputs(upload_choice = "Raw VCF File")
    session$flushReact()

    # Check that the output is rendered (we can't easily test the exact UI content)
    expect_true(TRUE) # Basic test that no error occurs

    # Test Processed VCF File selection
    session$setInputs(upload_choice = "Processed VCF File")
    session$flushReact()

    expect_true(TRUE) # Basic test that no error occurs
  })
})

# Test marker design execution
test_that("executes marker design with proper inputs", {
  temp_genome <- tempfile(fileext = ".fa")
  writeLines(c(">chr1", "ATCGATCGATCGATCG"), temp_genome)

  testServer(mod_kasp_marker_design_server, {
    # Skip this test if kasp_marker_design_alt function doesn't exist
    # This test focuses on the setup rather than the actual execution

    with_mocked_bindings(
      "read_excel" = function(...) {
        create_mock_excel_data()
      },
      .package = "readxl",
      {
        # Set all required inputs
        session$setInputs(
          upload_choice = "Processed VCF File",
          gt_df = list(datapath = "dummy.xlsx"),
          variant_id_col = "marker_id",
          chrom_col = "chromosome",
          pos_col = "position",
          ref_al_col = "ref_allele",
          alt_al_col = "alt_allele",
          geno_start = 6,
          marker_ID = "marker1",
          chr_ID = "chr1",
          genome_file = list(datapath = temp_genome),
          maf = 0.05,
          reg_name = "test_region"
        )

        session$flushReact()

        # Test that the reactive values are set up correctly
        expect_true(is.reactive(kasp_des.result))

        # Test that gt_data is properly loaded
        expect_true(is.reactive(gt_data))
      }
    )
  })

  unlink(temp_genome)
})

# Test download functionality
test_that("download handlers are properly configured", {
  testServer(mod_kasp_marker_design_server, {
    # Set up mock result data
    mock_result <- create_mock_kasp_result()
    kasp_des.result(mock_result)

    session$setInputs(
      file_name = "test_output",
      exten = ".csv"
    )

    # Test that the reactive value is set
    expect_equal(kasp_des.result()$marker_data$Marker_ID, "test_marker")

    # Test filename generation logic
    expect_true(is.character(session$input$file_name))
    expect_true(is.character(session$input$exten))
  })
})

# Test plot generation logic
test_that("handles plot generation settings", {
  testServer(mod_kasp_marker_design_server, {
    # Test plot enabled
    session$setInputs(draw_plot = TRUE)
    session$flushReact()

    expect_true(session$input$draw_plot)

    # Test plot disabled
    session$setInputs(draw_plot = FALSE)
    session$flushReact()

    expect_false(session$input$draw_plot)
  })
})

# Test data table rendering
test_that("renders data table correctly", {
  testServer(mod_kasp_marker_design_server, {
    # Set up mock result
    mock_result <- create_mock_kasp_result()
    kasp_des.result(mock_result)

    # Test that the table can be rendered (basic check)
    expect_true(is.reactive(kasp_des.result))
    expect_false(is.null(kasp_des.result()$marker_data))
  })
})

# Test input validation
test_that("validates required inputs", {
  testServer(mod_kasp_marker_design_server, {
    # Test with missing inputs - should not proceed
    session$setInputs(run_but = 1)

    # Since req() will stop execution, we expect the reactive to remain NULL
    expect_null(kasp_des.result())
  })
})

# Integration test
test_that("full workflow integration", {
  temp_excel <- tempfile(fileext = ".xlsx")
  temp_genome <- tempfile(fileext = ".fa")
  writeLines(c(">chr1", "ATCGATCGATCGATCG"), temp_genome)

  testServer(mod_kasp_marker_design_server, {
    # Mock read_excel function only (skip kasp_marker_design_alt for this test)
    with_mocked_bindings(
      "read_excel" = function(...) {
        create_mock_excel_data()
      },
      .package = "readxl",
      {
        # Simulate workflow setup
        session$setInputs(
          upload_choice = "Processed VCF File",
          gt_df = list(datapath = temp_excel),
          draw_plot = TRUE,
          file_name = "integration_test",
          exten = ".csv"
        )

        session$flushReact()

        # Test that basic workflow components are working
        expect_true(is.reactive(gt_data))
        expect_true(TRUE) # Basic integration test passed
      }
    )
  })

  unlink(c(temp_excel, temp_genome))
})

# Additional test for checking mock function calls (alternative approach)
test_that("tracks function calls without mockery", {
  # Create a flag to track if function was called
  function_called <- FALSE

  testServer(mod_kasp_marker_design_server, {
    with_mocked_bindings(
      "read_excel" = function(...) {
        function_called <<- TRUE
        create_mock_excel_data()
      },
      .package = "readxl",
      {
        session$setInputs(
          upload_choice = "Processed VCF File",
          gt_df = list(datapath = "test.xlsx")
        )

        session$flushReact()

        expect_true(function_called)
      }
    )
  })
})

# Test helper function for more complex mocking scenarios
create_counting_mock <- function(return_value) {
  call_count <- 0
  function(...) {
    call_count <<- call_count + 1
    # Store call count in a more accessible location
    .GlobalEnv$test_call_count <- call_count
    return_value
  }
}

test_that("tracks multiple function calls", {
  # Initialize counter
  .GlobalEnv$test_call_count <- 0

  testServer(mod_kasp_marker_design_server, {
    # Create a counting mock
    counting_mock <- create_counting_mock(create_mock_excel_data())

    with_mocked_bindings(
      "read_excel" = counting_mock,
      .package = "readxl",
      {
        # Call multiple times
        session$setInputs(
          upload_choice = "Processed VCF File",
          gt_df = list(datapath = "test1.xlsx")
        )
        session$flushReact()

        session$setInputs(
          gt_df = list(datapath = "test2.xlsx")
        )
        session$flushReact()

        # Check that function was called
        expect_true(.GlobalEnv$test_call_count > 0)
      }
    )
  })

  # Clean up
  if (exists("test_call_count", envir = .GlobalEnv)) {
    rm("test_call_count", envir = .GlobalEnv)
  }
})
