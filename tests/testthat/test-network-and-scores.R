test_that("build_exposure_network returns edges and igraph", {
  set.seed(1)
  d <- generate_dummy_exporisk(n = 40)
  al <- align_omics(d$microbiome, d$metabolome, d$exposures, d$meta, id_col = "sample_id", strict = TRUE)
  pr <- prep_omics(al$microbiome, al$metabolome, al$exposures)
  
  net <- build_exposure_network(pr$X, pr$Y, pr$E, fdr = 0.5, max_pairs = 500, seed = 1)
  expect_true(is.list(net))
  expect_true(all(c("edges","graph","meta") %in% names(net)))
  expect_true(inherits(net$graph, "igraph") || length(net$graph) == 0)
})

test_that("exposure_perturbation_score returns ranked data.frame", {
  set.seed(2)
  d <- generate_dummy_exporisk(n = 40)
  al <- align_omics(d$microbiome, d$metabolome, d$exposures, d$meta, id_col = "sample_id", strict = TRUE)
  pr <- prep_omics(al$microbiome, al$metabolome, al$exposures)
  
  sc <- exposure_perturbation_score(pr$X, pr$Y, pr$E, fdr = 0.5, max_pairs = 400, seed = 1)
  expect_true(is.data.frame(sc))
  expect_true(all(c("exposure","perturbation_score") %in% colnames(sc)))
})
