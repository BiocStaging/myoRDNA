file_path <- testthat::test_path(
  "fixtures",
  "filtered_meth_raw_chr1_chrRDNAm_small.rds"
)

loaded_data <- readRDS(file_path)

treatment = c(0,0,0,0,0,0,
              1,1,1,1)

num_samples <- length(treatment)

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

loaded_data <- new(
  "methylRawList",
  loaded_data,
  treatment = treatment
)

meth_rDNA <- get_rDNA_data(
  loaded_data,
  num_samples,
  rDNA_chr_name = "chrRDNAm")

cutoff = 35000
meth_rDNA_shifted <- shifting(
  meth_rdna_raw=meth_rDNA,
  cutoff=cutoff,
  treatment=treatment
)

# Test
test_that("test site by site differential methylation", {

  file_path <- testthat::test_path("fixtures", "differential_meth_site.rds")

  expect_true(nzchar(file_path))
  expect_true(file.exists(file_path))

  expected_output <- S4Vectors::DataFrame(readRDS(file_path))
  
  # Run the function
  actual_output <- diff_meth(
    meth_data=meth_rDNA_shifted,
    samples=samples,
    species="mm39",
    num_samples=num_samples,
    min_cov=10,
    region="18S",
    model_type="binomial",
    p_adjust_method="BH",
    site = TRUE)
  
  rownames(actual_output) <- NULL
  rownames(expected_output) <- NULL
  
  expect_s4_class(actual_output, "DataFrame")
  expect_equal(actual_output, expected_output)
})

# Test
test_that("test region by region differential methylation", {

  file_path <- testthat::test_path("fixtures", "differential_meth_region.rds")

  expect_true(nzchar(file_path))
  expect_true(file.exists(file_path))

  expected_output <- S4Vectors::DataFrame(readRDS(file_path))
  
  # Run the function
  actual_output <- diff_meth(
    meth_data=meth_rDNA_shifted,
    samples=samples,
    species="mm39",
    num_samples=num_samples,
    min_cov=10,
    region="18S",
    model_type="binomial", #beta-binomial, quasi-binomial
    p_adjust_method="BH",
    site = FALSE)
  
  expect_s4_class(actual_output, "DataFrame")
  
  # Check that both outputs contain the same columns
  expect_setequal(names(actual_output), names(expected_output))
  
  # Align columns by name
  expected_output <- expected_output[
    ,
    names(actual_output),
    drop = FALSE
  ]
  
  rownames(actual_output) <- NULL
  rownames(expected_output) <- NULL
  
  expect_s4_class(actual_output, "DataFrame")
  expect_equal(actual_output, expected_output)
})