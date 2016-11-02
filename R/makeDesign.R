#' @title Create a design matrix from a character vector. 
#' @description This function wraps around model.matrix to create a design matrix based on the different levels in model.matrix. Note that the order of the column names of the design matrix will be in lexicographic order according to the locale in use. For more, see ?Comparison
#' @param x Character vector to be used to create the design matrix.
#' @return A design matrix as the same type as returned using ?model.matrix.
#' @export
makeDesign <- function(x){

  col_names = sort(unique(x))
  design_mat = model.matrix(~ x + 0)
  colnames(design_mat) = col_names
  names(attr(design_mat, "contrasts")) = deparse(substitute(x))
  return(design_mat)

}
