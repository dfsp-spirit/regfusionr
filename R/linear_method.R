

#' @title Transform MNI305 coords (FreeSurfer fsaverage surface) to MNI152 coordinates using the linear method.
#'
#' @description This uses the 4x4 matrix from the FreeSurfer Coordinate Systems documentation.
#'
#' @param vertex_coords nx3 matrix of coordinates, e.g., typically from fsaverage surface vertices.
#'
#' @note This is the opposite of using the \cite{Wu et al.} approach: a linear transformation matrix is used. This approach is mainly implemented in this package to allow you to easily check the difference between the methods.
#'
#' @return nx3 numerical matrix if MNI152 coords.
#'
#' @export
linear_fsaverage_coords_to_MNI152_coords <- function(vertex_coords) {
  return(freesurferformats::doapply.transform.mtx(vertex_coords, mni152reg_mtx()));
}


#' @title Get fsaverage (MNI305) to MNI152 transformation matrix.
#'
#' @description This returns the 4x4 matrix from the FreeSurfer Coordinate Systems documentation.
#'
#' @note This is the opposite of using the \cite{Wu et al.} approach. It is mainly implemented in this package to allow you to easily check the difference between the methods.
#'
#' @keywords internal
mni152reg_mtx <- function() {
  return(matrix(c(0.9975, -0.0073, 0.0176, -0.0429, 0.0146, 1.0009, -0.0024, 1.5496, -0.0130, -0.0093, 0.9971, 1.1840, 0, 0, 0, 1), ncol = 4, byrow = TRUE));
}
