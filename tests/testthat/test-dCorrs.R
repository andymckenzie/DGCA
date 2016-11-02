context("running dCorrs")

ncols = 10
nrows = 8
length_names = 5 #needs to remain 5 to not conflict with fake id "foo" below
set.seed(42)
corrType = "pearson"

create_string <- function(length){
	string = paste(sample(letters, length, replace = TRUE),
		collapse = "")
	return(string)
}

row_names = make.unique(replicate(nrows, create_string(length_names)))
col_names = make.unique(replicate(ncols, create_string(length_names)))

simData_A = matrix(rnorm(ncols * nrows), nrow = nrows, ncol = ncols)
rownames(simData_A) = row_names
colnames(simData_A) = col_names

simData_B = matrix(rnorm(ncols * nrows), nrow = nrows, ncol = ncols)
rownames(simData_B) = row_names
colnames(simData_B) = col_names

test_that("pairwiseDCor works", {

	#create simulation data with a strong differential correlation signal
	simData_A2 = simData_A
	simData_B2 = simData_B
	simData_A2[ , 1] = c(2, 1.1, rep(1, times = (nrows - 2)))
	simData_A2[ , 2] = c(0, rep(1, times = (nrows - 1)))
	simData_B2[ , 1] = c(2, 1.1, rep(1, times = (nrows - 2)))
	simData_B2[ , 2] = c(2, rep(1, times = (nrows - 1)))

	#cor(c(1,1,1,1,1,1,1,1.1,2), c(1,1,1,1,1,1,1,1,2), method = "pearson")
	# [1] 0.9999507
	expect_equal(round(matCorr(simData_B2,
		corrType = "pearson")[1,2], digits = 2), 0.99)

	# cor(c(1,1,1,1,1,1,1,1.01,2), c(1,1,1,1,1,1,1,1,2), method = "spearman")
	# [1] 0.75592895
	expect_equal(round(matCorr(simData_B2,
		corrType = "spearman")[1,2], digits = 2), 0.76)

	fac = factor(c(rep("A", times = nrow(simData_A2)),
		rep("B", times = nrow(simData_B2))))
	design_mat = model.matrix(~ fac + 0)

	cor_res = getCors(cbind(t(simData_A2), t(simData_B2)), design_mat)

	PDC = pairwiseDCor(cor_res, compare = c("facA", "facB"))

	#expect strong differential correlation to be found here
	expect_equal(round(slot(PDC, "ZDiff")[1,2], 0), 9)
	expect_equal(round(log(slot(PDC, "PValDiff")[1,2], 10), 0), -18)

	expect_error(pairwiseDCor(cor_res, compare = c("foo", "bar")),
		"names are not in")

	#extract the resulting object with a pair slice
	dc_slice = dcTopPairs(PDC, 40, classify = FALSE)

	#gene name #2 in first row should correspond to the second colname
	expect_equal(dc_slice[1, "Gene2"], col_names[2])

	#top gene pair should be just as significant as before
	expect_equal(round(dc_slice[1, "zScoreDiff"], 0), 9)
	expect_equal(round(log(dc_slice[1, "pValDiff"],  10), 0), -18)

	dc_slice_class = dcTopPairs(PDC, 40, classify = TRUE)

	#test that classifying of the differential correlations works
	expect_equal(as.character(dc_slice_class[1, "Classes"]), "-/+")

	#can extract the resulting object with top pairs function
	dc_top = dcTopPairs(PDC, nPairs = 30)

	#top differentially correlated genes can be predicted
	expect_equal(as.character(dc_top[1, c("Gene1", "Gene2")]),
		c(col_names[1], col_names[2]))

})
