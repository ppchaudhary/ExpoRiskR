test_that("align_omics runs and returns aligned objects", {
  d <- generate_dummy_exporisk()

  aligned <- align_omics(
    microbiome = d$microbiome,
    metabolome = d$metabolome,
    exposures  = d$exposures,
    meta       = d$meta,
    id_col     = "sample_id",
    strict     = TRUE
  )

  expect_true(is.list(aligned) || inherits(aligned, "SummarizedExperiment"))
})
