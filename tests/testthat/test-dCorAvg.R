context("running dCorAvg test")

set.seed(42)

data(darmanis)
data(design_mat)
darmanis_test = darmanis[1:30, ]

test_that("differential correlation average works as expected", {

  ddcor_res_avg_gene = ddcorAll(inputMat = darmanis_test, design = design_mat, compare = c("oligodendrocyte", "neuron"), adjust = "perm", heatmapPlot = FALSE, nPerm = 5, nPairs = "all", getDCorAvg = TRUE, dCorAvgType = "gene_average", dCorAvgMethod = "median")

  expect_equal(round(ddcor_res_avg_gene[[2]]$avgZDiff[1], 1), 0.9)

  ddcor_res_avg_gene = ddcorAll(inputMat = darmanis_test, design = design_mat, compare = c("oligodendrocyte", "neuron"), adjust = "perm", heatmapPlot = FALSE, nPerm = 5, nPairs = "all", getDCorAvg = TRUE, dCorAvgType = "total_average", dCorAvgMethod = "median")

  expect_equal(round(ddcor_res_avg_gene[[2]][["total_zdiff"]], 1), 0.4)

  ddcor_res_avg_gene_mean = ddcorAll(inputMat = darmanis_test, design = design_mat, compare = c("oligodendrocyte", "neuron"), adjust = "perm", heatmapPlot = FALSE, nPerm = 5, nPairs = "all", getDCorAvg = TRUE, dCorAvgType = "total_average", dCorAvgMethod = "mean")

  expect_equal(round(ddcor_res_avg_gene_mean[[2]][["total_zdiff"]], 1), 0.3)

})
