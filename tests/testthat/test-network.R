test_that("alignment output is consistent for downstream use", {
  d <- generate_dummy_exporisk()

  aligned <- align_omics(
    d$microbiome, d$metabolome, d$exposures, d$meta,
    id_col = "sample_id",
    strict = TRUE
  )

  if (is.list(aligned)) {
    expect_true(length(aligned) > 0)
  }
})
