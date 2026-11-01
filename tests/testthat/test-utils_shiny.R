test_that("makeTraitInput generates consistent Shiny tag structures", {
  # Mock a Shiny namespace function
  ns <- function(x) x

  # Generate row for id 1 (should not have a trash button)
  res_1 <- panGenomeBreedr:::makeTraitInput(1, ns)

  # Generate row for id 2 (should have a trash button)
  res_2 <- panGenomeBreedr:::makeTraitInput(2, ns)

  # Verify structure
  expect_s3_class(res_1, "shiny.tag")
  expect_true(grepl("trait_1", as.character(res_1)))

  # Ensure only ID > 1 has the actionButton (trash icon)
  expect_false(grepl("remove_1", as.character(res_1)))
  expect_true(grepl("remove_2", as.character(res_2)))
})

test_that("render_ghost_plate_error produces a valid ggplot object", {
  plot <- panGenomeBreedr:::render_ghost_plate_error("Invalid file format")

  # Verify ggplot integrity
  expect_s3_class(plot, "ggplot")
  expect_equal(plot$labels$title, "PLATE LAYOUT ERROR")
})


