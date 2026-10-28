file_path <- testthat::test_path(
  "fixtures",
  "filtered_meth_raw_chr1_chrRDNAm_small.rds"
)

loaded_data <- readRDS(file_path)

treatment = c(0,0,0,0,0,0,
              1,1,1,1)

num_samples <- length(treatment)

test_that("test extract unique chromosome names", {

  file_path <- testthat::test_path("fixtures", "chromosome_names.rds")

  expect_true(nzchar(file_path))
  expect_true(file.exists(file_path))

  # Load the expected character vector
  expected_output <- readRDS(file_path)

  actual_output <- get_unique_chr_names(
    loaded_data,
    num_samples
  )

  expect_type(actual_output, "character")
  expect_equal(actual_output, expected_output)
})