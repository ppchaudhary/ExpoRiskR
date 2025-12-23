test_that("generate_dummy_exporisk returns expected components", {
  set.seed(1)
  d <- generate_dummy_exporisk(n = 20)
  expect_true(is.matrix(d$microbiome))
  expect_true(is.matrix(d$metabolome))
  expect_true(is.matrix(d$exposures))
  expect_true(is.data.frame(d$meta))
  expect_true(all(c("sample_id","outcome") %in% colnames(d$meta)))
})
