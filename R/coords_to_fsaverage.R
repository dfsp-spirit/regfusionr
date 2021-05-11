# Standalone function to covert MNI/Colin coordinates to fsaverage coords.

#' @title Map MNI152 coords to fsaverage coords.
#'
#' @inheritParams vol_to_fsaverage
#'
#' @param coords nx3 numeric matrix, the source RAS coordinates in the input image which must be in MNI152 space.
#'
#' @return nx3 numeric matrix of target coordinates.
#'
#' @export
mni152_coords_to_fsaverage_coords <- function(coords, template_type='MNI152', surface='white', fs_home=Sys.getenv("FS_HOME")) {
  check_coords(coords);

  if(nchar(fs_home) == 0) {
    stop("Parameter 'fs_home' must not be empty. Make sure that the environment variable FS_HOME is set or pass a valid path.");
  }
  if(! dir.exists(fs_home)) {
    stop(sprintf("Parameter 'fs_home' points to '%d', but that directory does not exist (or is not readable).", fs_home));
  }

  # Load surface
  lh_surf = freesurferformats::read.fs.surface(file.path(fs_home, 'subjects', 'fsaverage', 'surf', sprintf("lh.%s", surface)));
  rh_surf = freesurferformats::read.fs.surface(file.path(fs_home, 'subjects', 'fsaverage', 'surf', sprintf("rh.%s", surface)));

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
