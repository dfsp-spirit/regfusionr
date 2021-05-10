# Main functions for vol to surf mapping.

#' @title Get valid brain template spaces.
#' @keywords internal
valid_template_types <- function() c('MNI152_orig', 'Colin27_orig', 'MNI152_norm', 'Colin27_norm')


#' @title Get valid registration fusion types.
#' @keywords internal
valid_rf_types <- function() c('RF_ANTs', 'RF_M3Z')


#' @keywords internal
check_coords <- function(coords) {
  if(! is.matrix(coords)) {
    stop("Parameter 'coords' must be a numeric matrix.");
  } else {
    if(ncol(coords) != 3) {
      stop(sprintf("Parameter 'coords' must be a numeric matrix with 3 columns, has %d.", ncol(coords)));
    }
  }
}

#' @title Ensure the template and rf types and valid and are an allowed combination.
#'
#' @keywords internal
check_rf_and_template <- function(template_type, rf_type) {
  if(! (template_type %in% valid_template_types())) { stop(sprintf("Parameter template_type must be one of: '%s' but is '%s'.", paste(valid_template_types(), collapse=", "), template_type)); }
  if(! (rf_type %in% valid_rf_types())) { stop(sprintf("Parameter rf_type must be one of: '%s' but is '%s'.", paste(valid_rf_types(), collapse=", "), rf_type)); }

  if (rf_type == 'RF_ANTs') {
    accepted_template_types = c('MNI152_orig', 'Colin27_orig');
  } else {
    accepted_template_types = c('MNI152_norm', 'Colin27_norm');
  }
  if(! (template_type %in% accepted_template_types)) {
    stop(sprintf("When using rf_type '%s', template_type must be one of: '%s' but is '%s'.", rf_type,  paste(accepted_template_types, collapse=", "), template_type));
  }
}


#' @title Map MNI coords to fsaverage coords.
#'
#' @inheritParams vol_to_fsaverage
#'
#' @param coords nx3 numeric matrix, the source coordinates in the input image which must be in space 'template_type'.
#'
#' @return nx3 numeric matrix of target coordinates.
#'
#' @export
vol_coords_to_fsaverage_coords <- function(coords, template_type='MNI152_orig', rf_type='RF_ANTs') {
  check_coords(coords);
  check_rf_and_template(rf_type, template_type);
  stop("not implemented yet")

}

#' @title Map values from MNI volume to fsaverage surface.
#'
#' @description Applies the Wu et al. regfusion method to obtain surface coords, then interpolates values.
#'
#' @param input_img 3D or 4D NIFTI or MGZ image instance of type \code{fs.volume}. If 4D, the 4th dimension is considered the time/subject dimension.
#'
#' @param out_dir character string, the path to a writable output directory.
#'
#' @param template_type character string, the source template
#'
#' @param rf_type the regfusion type to use.
#'
#' @param interp interpolation method
#'
#' @param out_type character string, the format of the output files. One of the following: 'curv' for FreeSurfer curv format, 'mgz' for FreeSurfer MGZ format.
#'
#' @return named list of 2 character strings, the output files (for the 2 hemispheres) at keys 'lh' and 'rh'.
#'
#' @export
vol_to_fsaverage <- function(input_img, out_dir=".", template_type='MNI152_orig', rf_type='RF_ANTs', interp='linear', out_type='curv') {

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

  if(! dir.exists(out_dir)) {
    dir.create(out_dir, recursive = FALSE);
  }

  out_files = list('lh' = NULL, 'rh' = NULL);
  for (hemi in c('lh', 'rh')) {
    mapping_file = system.file("extdata", "mappings", sprintf("%s%s", hemi, mapping), package = "regfusionr", mustWork = TRUE);
    ras = t(as.matrix(data.table::fread(mapping_file, nrows = 3, header = FALSE)));
    affine = freesurferformats::mghheader.ras2vox(input_img$header);
    projected = project_data(input_img$data, affine, ras, interp);

    out_file = file.path(out_dir, sprintf("%s%s.%s", hemi, mapping, out_type));
    freesurferformats::write.fs.morph(out_file, projected);
    out_files[[hemi]] = out_file;
  }
  # mapping_file = file.path("~/develop/regfusionr/inst/extdata/lh.avgMapping_allSub_RF_ANTs_Colin27_orig_to_fsaverage.txt");
  return(out_files);
}


#' @title Perform trilinear data interpolation.
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
  coords = doapply.transform.mtx(ras, affine);

  if(is.null(data)) {
    stop("Parameter 'data' must not be NULL (expected numerical array).");
  }

  if(length(dim(data)) == 4L) {
    data = data[,,,1];
    warning("Using first slice of 4D data only."); # TODO: handle ALL slices (independently).
  };

  if(length(dim(data)) == 3L) {
    # https://rdrr.io/cran/oce/man/approx3d.html
    x = seq_len(dim(data)[1]);
    y = seq_len(dim(data)[2]);
    z = seq_len(dim(data)[3]);
    approx_dta = oce::approx3d(x, y, z, data, coords[,1], coords[,2], coords[,3]);
  } else {
    stop("Only 3D and 4D data supported.");
  }

  return(approx_dta);
}

