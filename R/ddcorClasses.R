#' An S4 class to store correlation matrices and associated info.
#'
#' @slot corMatList List of correlation, number of sample, and correlation significance matrices in each condition.
#' @slot design Design matrix inputted into the function.
#' @slot options Character vector of options used in the function.
#' @export
setClass ("corMats",
	representation(
	corMatList = "list",
	design = "matrix",
	options = "vector"
	)
)

#' S4 class for pairwise differential correlation matrices and associated info.
#'
#' @slot corA Correlation matrix for identifiers in condition A.
#' @slot corPvalA Matrix of correlation significances in condition A.
#' @slot corB Correlation matrix for identifiers in condition B.
#' @slot corPvalB Matrix of correlation significances in condition B.
#' @slot ZDiff Matrix of differences in z-values between conditions.
#' @slot PValDiff Matrix of p-values for differences in z-values between conditions.
#' @slot options Character vector of options used in the function.
#' @export
setClass ("dcPair",
	representation(
	corA = "matrix",
	corPvalA = "matrix",
	corB = "matrix",
	corPvalB = "matrix",
	ZDiff = "matrix",
	PValDiff = "matrix",
	options = "vector"
	)
)
