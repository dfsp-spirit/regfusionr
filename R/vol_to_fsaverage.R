# Main functions for vol to surf mapping.

#' @keywords internal
valid_template_types <- function() c('MNI152_orig', 'Colin27_orig', 'MNI152_norm', 'Colin27_norm')

#' @keywords internal
valid_rf_types <- function() c('RF_ANTs', 'RF_M3Z')

#' @keywords internal
check_rf_and_template <- function(template_type, rf_type) {
  if(! (template_type %in% valid_template_types())) { stop(sprintf("Parameter template_type must be one of: %s.", paste(valid_template_types(), collapse=", "))); }
  if(! (rf_type %in% valid_rf_types())) { stop(sprintf("Parameter rf_type must be one of: %s.", paste(valid_rf_types(), collapse=", "))); }

  if (rf_type == 'RF_ANTs') {
    accepted_template_types = c('MNI152_orig', 'Colin27_orig');
  } else {
    accepted_template_types = c('MNI152_norm', 'Colin27_norm');
  }
  if(! (template_type %in% accepted_template_types)) {
    stop(sprintf("When using rf_type '%s', template_type must be one of: %s.", rf_type,  paste(accepted_template_types, collapse=", ")));
  }
}

#' @title Map MNI coords to fsaverage coords.
#'
#'
#'
#' @export
vol_coords_to_fsaverage_coords <- function(mni_coords, template_type='MNI152_orig', rf_type='RF_ANTs') {
  check_rf_and_template(rf_type, template_type);

}

#' @title Map values from MNI volume to fsaverage surface.
#'
#' Applies the Wu et al. regfusion method to obtain surface coords, then interpolates values.
#'
#' @param input_img 3D or 4D NIFTI or MGZ image. If 4D, the 4th dimension is considered the time/subject dimension.
#'
#' @param out_dir character string, the path to a writeable output directory.
#'
#' @param template_type character string, the source template
#'
#' @param rf_type the regfusion type to use.
#'
#' @param interp interpolation method
#'
#' @param out_type character string, the format of the output files. One of the following: 'curv' for FreeSurfer curv format, 'mgz' for FreeSurfer MGZ format.
#'
#' @export
vol_to_fsaverage <- function(input_img, out_dir=".", template_type='MNI152_orig', rf_type='RF_ANTs', interp='linear', out_type='curv') {
  check_rf_and_template(rf_type, template_type);
  mapping = sprintf(".avgMapping_allSub_%s_%s_to_fsaverage.txt", rf_type, template_type);
  for (hemi in c('lh', 'rh')) {
    mapping_file = system.file("extdata", sprintf("%s%s", hemi, mapping), package = "regfusionr", mustWork = TRUE);
  }
  # mapping_file = file.path("~/develop/regfusionr/inst/extdata/lh.avgMapping_allSub_RF_ANTs_Colin27_orig_to_fsaverage.txt");
  ras = as.matrix(data.table::fread(mapping_file, nrows = 3, header = FALSE));

}
