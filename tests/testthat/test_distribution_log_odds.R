file_path <- testthat::test_path(
  "fixtures",
  "filtered_meth_raw_chr1_chrRDNAm_small.rds"
)

loaded_data <- readRDS(file_path)

# Samples
samples <- list(
  'OLD' = c(
    "SRR28431752", "SRR28431753", "SRR28431754", "SRR28431755", 
    "SRR28431756", "SRR28431757"
  ),
  'YOUNG' = c(
    "SRR28431748", "SRR28431749", "SRR28431750", "SRR28431751"
  )
)

treatment = c(0,0,0,0,0,0,
              1,1,1,1)

test_that("test log odds for numCs and numTs", {

  file_path <- testthat::test_path("fixtures", "log_odds.rds")
  
  expected_output <- S4Vectors::DataFrame(readRDS(file_path))
  
  actual_output <- distribution_log_odd(
    loaded_data,
    treatment,
    chromo = "chrRDNAm"
  )
  
  expect_s4_class(actual_output, "DataFrame")
  expect_equal(actual_output, expected_output)
})