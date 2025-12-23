test_that("align_omics aligns samples and prep_omics returns matrices", {
  set.seed(2)
  d <- generate_dummy_exporisk(n = 25)
  al <- align_omics(d$microbiome, d$metabolome, d$exposures, d$meta, id_col = "sample_id", strict = TRUE)
  pr <- prep_omics(al$microbiome, al$metabolome, al$exposures, na_action = "error")

  expect_true(is.matrix(pr$X))
  expect_true(is.matrix(pr$Y))
  expect_true(is.matrix(pr$E))
  expect_equal(nrow(pr$X), 25)
})
