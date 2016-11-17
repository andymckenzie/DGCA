context("running correlation matrix and associated matrices")

nums = matrix(c(
	c(0.1, 0.2, 0.3, 0.4, 0.5, 0.1, 0.2, 0.4, 0.3, 100),
	c(0.3, 0.1, 0.5, 0.3, 0.5, 0.1, 0.2, 0.4, 0.3, 100)),
	nrow = 2)

fac = factor(c(rep("a", 5), rep( "b", 5)))
design_mat = model.matrix(~ fac + 0)

gold = cor.test(nums[1,c(1:5)], nums[2,c(1:5)])

nums_mis = matrix(c(
	c(0.1, 0.2, 0.3, 0.4, 0.5, 0.1, 0.2, 0.4, 0.3, 100),
	c(0.3, 0.1, 0.5, 0.3, NA, 0.1, 0.2, 0.4, 0.3, 100)),
	nrow = 2)

gold_mis = cor.test(nums_mis[1,c(1:5)], nums_mis[2,c(1:5)],
	use = "pairwise.complete.obs")

test_that("getCor works", {

	expect_equal(slot(getCors(nums, design_mat), "corMatList")[["faca"]]$pvals[1,2], gold$p.value)

	#on data.frame inputs
	expect_equal(slot(getCors(data.frame(nums), design_mat), "corMatList")[["faca"]]$pvals[1,2], gold$p.value)

	gold_mis_spearman = suppressWarnings(cor.test(nums_mis[1,c(1:5)], nums_mis[2,c(1:5)],
		use = "pairwise.complete.obs", method = "spearman"))

	#can handle missing sample data?
	expect_equal(slot(getCors(nums_mis, design_mat, corrType = "spearman"),
		"corMatList")[["faca"]]$pvals[1,2], gold_mis_spearman$p.value)

})

test_that("design matrix works", {

	expect_equal(getGroupsFromDesign(nums, design_mat)[[2]], list("faca",
		"facb"))

	colnames(design_mat) = c("group1", "group2")

	expect_equal(getGroupsFromDesign(nums, design_mat)[[2]], list("group1",
		"group2"))

	expect_equal(names(slot(getCors(nums_mis, design_mat, corrType = "spearman"),
		"corMatList")), c("group1", "group2"))

})
