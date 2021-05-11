# Internal utility functions.

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

#' @keywords internal
check_affine <- function(affine) {
  if(! is.matrix(affine)) {
    stop("Parameter 'affine' must be a numeric matrix.");
  } else {
    if(ncol(affine) != 4) {
      stop(sprintf("Parameter 'affine' must be a numeric matrix with 4 columns, has %d.", ncol(affine)));
    }
    if(nrow(affine) != 4) {
      stop(sprintf("Parameter 'affine' must be a numeric matrix with 4 rows, has %d.", nrow(affine)));
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

#' @keywords internal
get_data_file <- function(filename, subdir=NULL) {
  if(is.null(subdir)) {
    data_file = system.file("extdata", filename, package = "regfusionr", mustWork = TRUE);
  } else {
    data_file = system.file("extdata", subdir, filename, package = "regfusionr", mustWork = TRUE);
  }
  return(data_file);
}
