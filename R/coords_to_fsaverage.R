# Standalone function to covert MNI/Colin coordinates to fsaverage coords.

#' @title Map MNI coords to fsaverage coords.
#'
#' @inheritParams vol_to_fsaverage
#'
#' @param coords nx3 numeric matrix, the source coordinates in the input image which must be in space 'template_type'.
#'
#' @return nx3 numeric matrix of target coordinates.
#'
#' @export
vol_coords_to_fsaverage_coords <- function(coords, template_type='MNI152', rf_type='RF_ANTs') {
  check_coords(coords);
  check_rf_and_template(rf_type, template_type);

  # Load mappings
  lh_map_file = get_data_file("FSL_MNI152_FS4.5.0_RF_ANTs_avgMapping.vertex.lh.mgz", subdir = "coordmap");
  rh_map_file = get_data_file("FSL_MNI152_FS4.5.0_RF_ANTs_avgMapping.vertex.rh.mgz", subdir = "coordmap");
  lh_vertex = freesurferformats::read.fs.mgh(lh_map_file, with_header = FALSE, drop_empty_dims = TRUE);
  rh_vertex = freesurferformats::read.fs.mgh(rh_map_file, with_header = FALSE, drop_empty_dims = TRUE);

  # Load cortex mask.
  cortex_mask_filename = "FSL_MNI152_FS4.5.0_cortex_estimate.nii.gz";
  #cortex_mask_filename_colin27 = "SPM_Colin27_FS4.5.0_cortex_estimate.nii.gz";
  cortex_mask_file = get_data_file(cortex_mask_filename, subdir = "coordmap");
  cortex_mask = freesurferformats::read.fs.volume(cortex_mask_file);

  #

}
