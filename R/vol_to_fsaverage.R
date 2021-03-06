# Main functions for vol to surf mapping for an input 3D/4D image.


#' @title Project or map values from MNI volume to fsaverage surface.
#'
#' @description Applies the Wu et al. regfusion method to obtain surface coords, then interpolates values.
#'
#' @param input_img 3D or 4D NIFTI or MGZ image instance of type \code{fs.volume}, or a character string that will be interpreted as a file system path to such a volume that should be loaded with \code{freesurferformats::read.fs.volume}. If 4D, the 4th dimension is considered the time/subject dimension.
#'
#' @param template_type character string, the source template or the space that your input image is in. One of 'MNI152_orig', 'Colin27_orig', 'MNI152_norm', 'Colin27_norm'.
#'
#' @param rf_type the regfusion type to use, one of 'RF_ANTs' or 'RF_M3Z'.
#'
#' @param interp interpolation method, currently only 'linear' is supported.
#'
#' @param out_type character string, the format of the output files. One of the following: 'curv' for FreeSurfer curv format, 'mgz' for FreeSurfer MGZ format, 'gii' for GIFTI format.
#'
#' @param out_dir character string, the path to a writable output directory. If \code{NULL}, the returned named list contains the projected data (instead of the path of the file it was written to) at the keys 'lh' and 'rh', and the parameter 'out_type' is ignored.
#'
#' @return named list of 2 character strings, the output files (for the 2 hemispheres) at keys 'lh' and 'rh'. See the documentation for parameter 'out_dir' if you want the data in R instead.
#'
#' @author Tim Schäfer for the R version, Wu Jianxiao and CBIG for the original Matlab version.
#'
#' @importFrom data.table fread
#'
#' @export
vol_to_fsaverage <- function(input_img, template_type, rf_type='RF_ANTs', interp='linear', out_type='curv', out_dir=".") {

  if(! freesurferformats::is.fs.volume(input_img)) {
    if(is.character(input_img)) {
      input_img = freesurferformats::read.fs.volume(input_img, with_header = TRUE);
    } else {
      stop("Parameter 'input_img' must be fs.volume instance or a valid path to a volume file.");
    }
  }

  check_rf_and_template(template_type = template_type, rf_type = rf_type);

  valid_out_types = c('curv', 'mgz', 'gii');
  if(! (out_type %in% valid_out_types)) {
    stop(sprintf("Parameter 'out_type' must be one of: %s.", paste(valid_out_types(), collapse=", ")));
  }

  mapping = sprintf(".avgMapping_allSub_%s_%s_to_fsaverage.txt", rf_type, template_type);

  if(! is.null(out_dir)) {
    if(! dir.exists(out_dir)) {
      dir.create(out_dir, recursive = FALSE);
    }
  }

  out = list('lh' = NULL, 'rh' = NULL);
  for (hemi in c('lh', 'rh')) {
    mapping_file = system.file("extdata", "mappings", sprintf("%s%s", hemi, mapping), package = "regfusionr", mustWork = TRUE);
    ras = t(as.matrix(data.table::fread(mapping_file, nrows = 3, header = FALSE)));
    affine = freesurferformats::mghheader.ras2vox(input_img$header);
    projected = project_data(input_img$data, affine, ras, interp);

    if(is.null(out_dir)) {
      out[[hemi]] = projected;
    } else {
      if(is.array(projected) && out_type %in% c("curv", "gii")) {
        stop("The 'curv' and 'gii' output formats are not supported for 4D input data, try 'mgz' instead.");
        # We could write one output curv/gifti file per frame, but I guess simply using MGZ is better.
      }
      out_file = file.path(out_dir, sprintf("%s%s.%s", hemi, mapping, out_type));
      freesurferformats::write.fs.morph(out_file, projected);
      out[[hemi]] = out_file;
    }
  }

  return(out);
}


#' @title Perform coord transform and trilinear data interpolation.
#'
#' @param data 3D or 3D source image values
#'
#' @param affine 4x4 double matrix, the ras2vox transformation matrix
#'
#' @param ras nx3 double matrix, the source RAS coordinates
#'
#' @param interp character string, the interpolation mode. Only 'linear' is currently supported. We should support 'nearest' as well.
#'
#' @return a numerical vector for 3D input data or a numerical 4D array for 4D input data, the data values interpolated at the RAS coordinates. For the 4D input case, the 2D output is encoded in a 4D array with two dimensions of length 1 as follows: The first dimension contains the projected and interpolated per-vertex data values for a single frame, dimensions 2 and 3 are both of length 1, and the fourth dimensions is the timepoint/subject dimension (the frames).
#'
#' @keywords internal
project_data <- function(data, affine, ras, interp='linear') {

  if(requireNamespace("oce", quietly = TRUE)) {

    supported_interp = c("linear");
    if(!(interp %in% supported_interp)) {
      stop("Parameter 'interp' must be 'linear', others not supported yet.");
    }

    if(is.null(data)) {
      stop("Parameter 'data' must not be NULL (expected numerical array).");
    }

    check_affine(affine);

    coords = freesurferformats::doapply.transform.mtx(ras, affine);

    data = drop(data); # remove empty dimensions (typically the last one)

    if(length(dim(data)) == 3L) {
      nvols = 1L;
      x = seq_len(dim(data)[1]);
      y = seq_len(dim(data)[2]);
      z = seq_len(dim(data)[3]);
      proj_data = oce::approx3d(x, y, z, data, coords[,1], coords[,2], coords[,3]);
    } else if (length(dim(data)) == 4L) {
      nvols = dim(data)[4];
      proj_data = array(data = rep(NA, (nrow(ras)*nvols)), dim = c(nrow(ras), 1L, 1L, nvols));
      for(vol_idx in seq_len(nvols)) {
        x = seq_len(dim(data)[1]);
        y = seq_len(dim(data)[2]);
        z = seq_len(dim(data)[3]);
        approx_dta = oce::approx3d(x, y, z, data[,,,vol_idx], coords[,1], coords[,2], coords[,3]);
        proj_data[,1L,1L,vol_idx] = approx_dta;
      }
    } else {
      stop("Only 3D or 4D input data supported.");
    }
    return(proj_data);
  } else {
    stop("The 'oce' package must be installed to use this functionality. See https://github.com/dfsp-spirit/regfusionr for installation instructions.");
  }
}

