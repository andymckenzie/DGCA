context("running filtering of matrix")

set.seed(42)

high_var = data.frame(matrix(rnorm(140, mean = 10, sd = 10), nrow = 7))
low_var = data.frame(matrix(rnorm(140, mean = 10, sd = 1), nrow = 7))
low_mean = data.frame(matrix(rnorm(120, mean = 1, sd = 10), nrow = 6))
total_mat = rbind(high_var, low_var, low_mean)
rownames(total_mat) = letters[1:20]

groups = c(rep(c("A", "B", "C"), each = 6), "A")
design = makeDesign(groups)

test_that("filterGenes works", {

  library(matrixStats)

  filtered_mean = filterGenes(total_mat, filterTypes = c("central"),
    keepRows = NULL, filterCentralType = "median", filterCentralPercentile = 0.3)

  #filtering by mean works as expected
  expect_equal(rownames(filtered_mean), letters[1:14])

  filtered_mean = filterGenes(total_mat, filterTypes = c("central"),
    keepRows = NULL, filterCentralType = "mean", filterCentralPercentile = 0.3)

  #filtering by mean works as expected
  expect_equal(rownames(filtered_mean), letters[1:14])

  filtered_variance = filterGenes(total_mat, filterTypes = c("dispersion"),
    keepRows = NULL, filterDispersionType = "cv", filterDispersionPercentile = 0.35)

  #filtering by variance works as expected
  expect_equal(rownames(filtered_variance), c(letters[1:7], letters[15:20]))

  filtered_variance = filterGenes(total_mat, filterTypes = c("dispersion"),
    keepRows = NULL, filterDispersionType = "variance", filterDispersionPercentile = 0.35)

  #filtering by variance works as expected
  expect_equal(rownames(filtered_variance), c(letters[1:7], letters[15:20]))

  filtered_mat = filterGenes(total_mat, filterTypes = c("central", "dispersion"),
    keepRows = NULL, filterCentralType = "median",
    filterDispersionType = "cv", filterCentralPercentile = 0.3,
    filterDispersionPercentile = 0.5, sequential = TRUE)

  #filtering by mean and variance sequentially works as expected
  expect_equal(rownames(filtered_mat), c(letters[1:7]))

  filtered_mat = filterGenes(total_mat, filterTypes = c("central", "dispersion"),
    keepRows = NULL, filterCentralType = "median",
    filterDispersionType = "cv", filterCentralPercentile = 0.3,
    filterDispersionPercentile = 0.35, sequential = FALSE)

  #filtering by mean and variance independently works as expected
  expect_equal(rownames(filtered_mat), c(letters[1:7]))

  filtered_mat = filterGenes(total_mat, filterTypes = c("central", "dispersion"),
    keepRows = c(letters[10], letters[20]), filterCentralType = "median",
    filterDispersionType = "cv", filterCentralPercentile = 0.3,
    filterDispersionPercentile = 0.45, sequential = TRUE)

  #keeprows works as expected for sequential filtering
  expect_equal(rownames(filtered_mat), c(letters[1:7], letters[10], letters[20]))

  filtered_mat = filterGenes(total_mat, filterTypes = c("central", "dispersion"),
    keepRows = c(letters[10], letters[20]), filterCentralType = "median",
    filterDispersionType = "cv", filterCentralPercentile = 0.3,
    filterDispersionPercentile = 0.35, sequential = FALSE)

  #keeprows works as expected for non-sequential filtering
  expect_equal(rownames(filtered_mat), c(letters[1:7], letters[10], letters[20]))

  filtered_mat = filterGenes(total_mat, filterTypes = c("central"), keepRows = NULL, filterCentralType = "mean", filterCentralPercentile = 0.35, allGroups = TRUE, design = design)

  #keeprows works as expected given multiple groups
  expect_equal(rownames(filtered_mat), c("b", "d", "e", "f", "h", "i", "j", "k", "l", "m", "n"))

})
