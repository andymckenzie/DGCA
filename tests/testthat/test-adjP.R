context("running adjustPVals")

set.seed(42)

pvals = runif(0, 1, n = 100)

test_that("adjPVals works as expected on improper input", {

	expect_error(adjustPVals(pvals, adjust = "foo"),
		"Adjust method is one of the available methods")

})

test_that("adjPVals works as expected on proper input", {

	expect_message(adjustPVals(pvals, adjust = "none", verbose = TRUE),
		"No adjustment of p-values")

	expect_message(adjustPVals(pvals, adjust = "bonferroni", verbose = TRUE),
		"p.adjust")

	expect_gt(adjustPVals(pvals, adjust = "bonferroni")[50], pvals[50])

})

test_that("bigEmpPVals handles normal test stats as well as thresholded test stats", {

	test_stat = rnorm(100, 0, 1)
	test0_stat = rnorm(1000, 0, 1)
	emp_pvals = bigEmpPVals(test_stat, test0_stat)

	expect_gt(round(max(emp_pvals), 1), 0.8)

	expect_lt(round(min(emp_pvals), 1), 0.2)

	test_stat = rnorm(100, 0.75, 1)
	test_stat[test_stat < 0] = 0
	test0_stat = rnorm(1000, 0.4, 1)
	test0_stat[test0_stat < 0] = 0

	emp_pvals_ties = bigEmpPVals(test_stat, test0_stat)

	expect_gt(round(max(emp_pvals_ties), 1), 0.8)

	expect_lt(round(min(emp_pvals_ties), 1), 0.2)

})
