# Main functions for vol to surf mapping for an input 3D/4D image.


#' @title Map values from MNI volume to fsaverage surface.
#'
#' @description Applies the Wu et al. regfusion method to obtain surface coords, then interpolates values.
#'
#' @param input_img 3D or 4D NIFTI or MGZ image instance of type \code{fs.volume}. If 4D, the 4th dimension is considered the time/subject dimension.
#'
#' @param template_type character string, the source template or the space that your input image is in. One of 'MNI152_orig', 'Colin27_orig', 'MNI152_norm', 'Colin27_norm'.
#'
#' @param rf_type the regfusion type to use, one of 'RF_ANTs' or 'RF_M3Z'.
#'
#' @param interp interpolation method, currenlty only 'linear' is supported.
#'
#' @param out_type character string, the format of the output files. One of the following: 'curv' for FreeSurfer curv format, 'mgz' for FreeSurfer MGZ format.
#'
#' @param out_dir character string, the path to a writable output directory. If \code{NULL}, the returned named list contains the projected data (instead of the path of the file it was written to), and the parameter 'out_type' is ignored.
#'
#' @return named list of 2 character strings, the output files (for the 2 hemispheres) at keys 'lh' and 'rh'. See 'out_dir' if you want the data in R instead.
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

  valid_out_types = c('curv', 'mgz');
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
#' @return data values interpolated at the RAS coordinates.
#'
#' @keywords internal
project_data <- function(data, affine, ras, interp='linear') {

  supported_interp = c("linear");
  if(!(interp %in% supported_interp)) {
    stop("Parameter 'interp' must be 'linear', others not supported yet.");
  }

  if(is.null(data)) {
    stop("Parameter 'data' must not be NULL (expected numerical array).");
  }

  check_affine(affine);

  coords = doapply.transform.mtx(ras, affine);

  if(length(dim(data)) == 3L) {
    nvols = 1L;
    proj_data = array(data = rep(NA, (nrow(ras)*nvols)), dim = c(nvols, nrow(ras)));
    x = seq_len(dim(data)[1]);
    y = seq_len(dim(data)[2]);
    z = seq_len(dim(data)[3]);
    approx_dta = oce::approx3d(x, y, z, data, coords[,1], coords[,2], coords[,3]);
    proj_data[1,] = approx_dta;
  } else if (length(dim(data)) == 4L) {
    nvols = dim(data)[4];
    proj_data = array(data = rep(NA, (nrow(ras)*nvols)), dim = c(nvols, nrow(ras)));
    for(vol_idx in seq_len(nvols)) {
      x = seq_len(dim(data)[1]);
      y = seq_len(dim(data)[2]);
      z = seq_len(dim(data)[3]);
      approx_dta = oce::approx3d(x, y, z, data[,,,vol_idx], coords[,1], coords[,2], coords[,3]);
      proj_data[vol_idx,] = approx_dta;
    }
  } else {
    stop("Only 3D and 4D data supported.");
  }

  return(approx_dta);
}

