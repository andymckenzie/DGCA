context("running correlation matrix and associated matrices")

nums = matrix(c(0.1, 0.2, 0.3, 0.4, 0.5, 0.1, 0.2, 0.4, 0.3, 100), nrow = 5)
cors = cor(nums)
nsamps = matrix(c(5, 5, 5, 5), nrow = 2)
gold = cor.test(nums[,1], nums[,2])

nums_mis = matrix(c(0.1, 0.2, 0.5, NA, 0.3, 0.1, 0.2, 0.3, 0.4, 100), nrow = 5)
nsamps_mis = matrix(c(4, 4, 4, 5), nrow = 2)
cors_mis = cor(nums_mis, use = "pairwise.complete.obs")
gold_mis = cor.test(nums_mis[,1], nums_mis[,2], use = "pairwise.complete.obs")

test_that("matCorrSig works", {

	expect_equal(matCorSig(cors, nsamps)[1,2], gold$p.value)

	#can handle missing sample data?
	expect_equal(matCorSig(cors_mis, nsamps_mis)[1,2], gold_mis$p.value)

})


test_that("matCorr works", {

	expect_equal(matCorr(nums, corrType = "pearson"), cor(nums))

	expect_equal(matCorr(nums_mis, corrType = "spearman"), cor(nums_mis,
		method = "spearman", use = "pairwise.complete.obs"))

})

test_that("matNSamp works", {

	expect_equal(matNSamp(nums), nsamps)

	expect_equal(matNSamp(nums_mis), nsamps_mis)

	expect_equal(matNSamp(nums_mis, impute = TRUE), nsamps)

})
