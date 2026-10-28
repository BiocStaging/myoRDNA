file_path <- testthat::test_path(
  "fixtures",
  "filtered_meth_raw_chr1_chrRDNAm_small.rds"
)

loaded_data <- readRDS(file_path)

treatment = c(0,0,0,0,0,0,
              1,1,1,1)

num_samples <- length(treatment)

test_that("test calculate the relative rDNA copy number", {

  file_path <- testthat::test_path("fixtures", "rel_out.rds")
  
  expected_output <- S4Vectors::DataFrame(readRDS(file_path))
  
  actual_output <- relative_rDNA_CN(
    Meth_object = loaded_data,
    num_samples = num_samples,
    autosomes = TRUE,
    rDNA_chr = "chrRDNAm",
    diploid = TRUE
  )
  
  expect_s4_class(actual_output, "DataFrame")
  expect_equal(actual_output, expected_output)
})