context("running ddMEGENA")

set.seed(42)

data(darmanis)
data(design_mat)
darmanis_test = darmanis[1:30, ]

test_that("full ddcorAll process yields correct results on this test data set", {

	testthat::skip(message = "skipping MEGENA integration test by default to save compute time and avoid output")

	ddcor_res = ddcorAll(inputMat = darmanis_test, design = design_mat,
		compare = c("oligodendrocyte", "neuron"),
		adjust = "none", heatmapPlot = FALSE, nPerm = 0, nPairs = "all")

	megena_res = suppressWarnings(ddMEGENA(ddcor_res, adjusted = FALSE, evalCompactness = FALSE, nPerm = 10, minModSize = 10))

	expect_equal(megena_res$modules[1, "Genes"], "ACSL3")

	megena_res = suppressWarnings(ddMEGENA(ddcor_res, adjusted = FALSE, evalCompactness = TRUE, nPerm = 10, minModSize = 10, pval_gene_thresh = 0.5))

	expect_equal(megena_res$modules[2, "Genes"], "AASDHPPT")

})
