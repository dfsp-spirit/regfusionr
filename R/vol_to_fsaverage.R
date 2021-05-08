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
#' @export
vol_to_fsaverage <- function(input_img, out_dir, template_type='MNI152_orig', rf_type='RF_ANTs', interp='linear', out_type='nii.gz') {
  check_rf_and_template(rf_type, template_type);

}
