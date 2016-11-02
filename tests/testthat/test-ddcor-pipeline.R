context("running full ddcorAll pipeline")

set.seed(42)

data(darmanis)
data(design_mat)
darmanis_test = darmanis[1:30, ]

test_that("full ddcorAll process yields correct results on this test data set", {

  ddcor_res = ddcorAll(inputMat = darmanis_test, design = design_mat,
    compare = c("oligodendrocyte", "neuron"),
    adjust = "none", heatmapPlot = FALSE, nPerm = 0, nPairs = "all")

  ddcor_res_spearman = ddcorAll(inputMat = darmanis_test, design = design_mat,
    compare = c("oligodendrocyte", "neuron"), corrType = "spearman",
    adjust = "none", heatmapPlot = FALSE, nPerm = 0, nPairs = "all")

  expect_equal(ddcor_res[1, ]$Gene1, "ACSL3")

})

test_that("p-value adjustments work", {

  ddcor_res = ddcorAll(inputMat = darmanis_test, design = design_mat, compare = c("oligodendrocyte", "neuron"), adjust = "none", heatmapPlot = FALSE, nPerm = 0)

  expect_equal(sum("pValDiff_adj" %in% colnames(ddcor_res)), 0)

  ddcor_res = ddcorAll(inputMat = darmanis_test, design = design_mat, compare = c("oligodendrocyte", "neuron"),  adjust = "BH", nPerm = 0)

  expect_equal(sum("pValDiff_adj" %in% colnames(ddcor_res)), 1)

  expect_equal(round(ddcor_res$pValDiff_adj[1], 1), 0)

  ddcor_res = ddcorAll(inputMat = darmanis_test, design = design_mat, compare = c("oligodendrocyte", "neuron"), adjust = "perm", heatmapPlot = FALSE, nPerm = 5)

  expect_equal(sum("empPVals" %in% colnames(ddcor_res)), 1)

})

test_that("split set works", {

  ddcor_res = ddcorAll(inputMat = darmanis_test, design = design_mat, compare = c("oligodendrocyte", "neuron"), adjust = "none", heatmapPlot = FALSE, nPerm = 0, splitSet = c("ACSL3", "ARHGAP5"), nPairs = "all")

	expect_equal(nrow(ddcor_res), (nrow(darmanis_test)-2) * 2)

	expect_true(all(ddcor_res$Gene2 %in% c("ACSL3", "ARHGAP5")))

  darmanis_test_a = darmanis_test[!rownames(darmanis_test) %in% c("ACSL3", "ARHGAP5"), ]
  darmanis_test_b = darmanis_test[c("ACSL3", "ARHGAP5"), ]

  ddcor_res_b = ddcorAll(inputMat = darmanis_test_a, design = design_mat, inputMatB = darmanis_test_b, compare = c("oligodendrocyte", "neuron"), adjust = "none", heatmapPlot = FALSE, nPerm = 0, nPairs = "all")

  expect_equal(ddcor_res, ddcor_res_b)

})

test_that("corr cutoff works", {

  ddcor_res = ddcorAll(inputMat = darmanis_test, design = design_mat,
    compare = c("oligodendrocyte", "neuron"), adjust = "none",
    heatmapPlot = FALSE, nPerm = 0, corr_cutoff = 0.7)

  #usually the zScore diff is greater than 5, but with the corr cutoff to 0.7, this is affected.
	expect_lt(ddcor_res[1, ]$zScoreDiff, 5)

})

test_that("making design matrices works", {

  n_oligo_samples = 38; n_neuron_samples = 120
  cell_type = c(rep("oligodendrocyte", n_oligo_samples), rep("neuron", n_neuron_samples))
  design_mat = model.matrix(~ cell_type + 0)
  colnames(design_mat) = c("neuron", "oligodendrocyte")

  design_mat2 = makeDesign(cell_type)

  expect_equal(design_mat, design_mat2)

})
