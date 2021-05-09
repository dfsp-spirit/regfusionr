
#' @title Apply spatial transformation matrix to input coords/voxel indices.
#'
#' @note The current implementation is copied from the develop version of freesurferformats. It will be called from there once the function is in the release.
#'
#' @keywords internal
ras_to_vox <- function(ras, affine) {
  #freesurferformats::doapply.transform.mtx(ras, affine)
  doapply.transform.mtx(ras, affine)
}

#' @title Apply a spatial transformation matrix to the given coordinates.
#'
#' @param coords nx3 (cartesian) or nx4 (homogeneous) numerical matrix, the input coordinates. If nx4, left as is for homogeneous notation, if nx3 (cartesian) a 1 will be appended as the 4th position.
#'
#' @param mtx a 4x4 numerical transformation matrix
#'
#' @return the coords after applying the transformation. If coords was nx3, nx3 is returned, otherwise nx4.
#'
#' @note This function is copied from the develop version of freesurferformats. It will be called from there once the function is in the release.
#'
#' @examples
#'     coords_tf = doapply.transform.mtx(c(1.0, 1.0, 1.0), mni152reg());
#'     coords_tf;
#'     doapply.transform.mtx(coords_tf, solve(mni152reg()));
#'
#' @export
doapply.transform.mtx <- function(coords, mtx) {
  if(is.vector(coords)) {
    coords = matrix(coords, nrow = 1L);
  }
  was_cartesian = FALSE;
  if(ncol(coords) == 3L) {
    was_cartesian = TRUE;
    coords = cbind(coords, 1.0);
  }

  if(ncol(coords) != 4L) {
    stop(sprintf("Parameter coords must have 3 or 4 colums (or 3 or 4 entries for a vector), found %d.\n", ncol(coords)));
  }

  if(any(coords[,4] == 0)) {
    warning("Invalid 'coords': last column must not be zero.");
  }

  transformed_coords = t(mtx %*% t(coords));
  if(was_cartesian) {
    transformed_coords = transformed_coords[,1:3] / transformed_coords[,4]; # convert from homogeneous to cartesian
  }
  return(transformed_coords);
}


project_data <- function() {

}


