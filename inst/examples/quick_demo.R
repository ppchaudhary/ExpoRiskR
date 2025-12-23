library(ExpoRiskR)

d <- generate_dummy_exporisk()
al <- align_omics(d$microbiome, d$metabolome, d$exposures, d$meta, id_col = "sample_id")
pr <- prep_omics(al$microbiome, al$metabolome, al$exposures)

print(dim(pr$X))
