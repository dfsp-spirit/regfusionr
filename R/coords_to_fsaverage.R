# Standalone function to covert MNI/Colin coordinates to fsaverage coords.

#' @title Map MNI152 coords to fsaverage coords.
#'
#' @param coords nx3 numeric matrix, the source RAS coordinates in the input image which must be in MNI152 space.
#'
#' @param surface character string, the fsaverage surface (brain mesh) to load. Must be a valid FreeSurfer surface name like 'white', 'pial', 'orig, 'inflated'.
#'
#' @param fs_home character string, path to the FreeSurfer installation. Alternatively, a hemilist of \code{freesurferformats::fs.surface} instances like \code{surface = list("lh"=mysurflh, "rh"=mysurfrh)}. Used to find the surfaces, at \code{<fs_home>/subjects/fsaverage/surf/<hemi>.<surface>}, where hemi is 'lh' and 'rh'. Can be NULL if 'surface' is a hemilist of fs.surface instances.
#'
#' @return named list with entries 'fsaverage_vertices': integer vector of fsaverage surface vertex indices, 'hemi': vector of hemi strings for the vertices, 'fsaverage_coords': nx3 numeric matrix of target coordinates.
#'
#' @note see standalone_scripts_for_MNI_fsaverage_coordinates_conversion/CBIG_RF_MNICoord2fsaverageVertex.m
#'
#' @export
mni152_coords_to_fsaverage <- function(coords, surface='white', fs_home=Sys.getenv("FS_HOME")) {
  if(is.vector(coords)) {
    if(length(coords) %% 3L == 0L) {
      coords = matrix(coords, ncol = 3L);
    }
  }
  check_coords(coords);

  # Load surfaces
  if(is.list(surface)) {
    if(freesurferformats::is.fs.surface(surface$lh) & freesurferformats::is.fs.surface(surface$rh)) {
      lh_surf = surface$lh;
      rh_surf = surface$rh;
    } else {
      stop("Parameter 'surface' must be a character string like 'white' or a hemilist of fs.surface instances.");
    }
  } else if (is.character(surface)) {
    if(nchar(fs_home) == 0) {
      stop("Parameter 'fs_home' must not be empty. Make sure that the environment variable FS_HOME is set or pass a valid path.");
    }
    if(! dir.exists(fs_home)) {
      stop(sprintf("Parameter 'fs_home' points to '%s', but that directory does not exist (or is not readable).", fs_home));
    }
    fsavg_path = file.path(fs_home, 'subjects', 'fsaverage');
    if(! dir.exists(fsavg_path)) {
      stop(sprintf("Parameter 'fs_home' points to '%s', but expected fsaverage sub directory '%s' does not exist.", fs_home, fsavg_path));
    }

    lh_surf = freesurferformats::read.fs.surface(file.path(fs_home, 'subjects', 'fsaverage', 'surf', sprintf("lh.%s", surface)));
    rh_surf = freesurferformats::read.fs.surface(file.path(fs_home, 'subjects', 'fsaverage', 'surf', sprintf("rh.%s", surface)));
  } else {
    stop("Parameter 'surface' must be a character string like 'white' or a hemilist of fs.surface instances.");
  }

  # Load mappings
  lh_map_file = get_data_file("FSL_MNI152_FS4.5.0_RF_ANTs_avgMapping.vertex.lh.mgz", subdir = "coordmap");
  rh_map_file = get_data_file("FSL_MNI152_FS4.5.0_RF_ANTs_avgMapping.vertex.rh.mgz", subdir = "coordmap");
  lh_vertex = freesurferformats::read.fs.mgh(lh_map_file, with_header = FALSE, drop_empty_dims = TRUE);
  rh_vertex = freesurferformats::read.fs.mgh(rh_map_file, with_header = FALSE, drop_empty_dims = TRUE);

  # Load cortex mask.
  cortex_mask_file = get_data_file("FSL_MNI152_FS4.5.0_cortex_estimate.nii.gz", subdir = "coordmap");
  cortex_mask = freesurferformats::read.fs.volume(cortex_mask_file, with_header = TRUE);

  # Do the masking.
  lh_vertex[which(cortex_mask$data == 0)] = 0;
  rh_vertex[which(cortex_mask$data == 0)] = 0;

  # Convert input RAS coords to voxel indices (IJK) for the matrix.
  mni_voxels = doapply.transform.mtx(coords, freesurferformats::mghheader.ras2vox(cortex_mask)) + 1L;
  mni_array = array(data = c(mni_voxels[2,], mni_voxels[1,], mni_voxels[3,]) , dim = c(256, 256, 256)); # TODO: get dim from image

  num_coords = nrow(coords);
  verts = rep(0L, num_coords);
  fs_coords = matrix(rep(0.0, (num_coords * 3L)), ncol = 3L);
  hemi = rep(NULL, num_coords);
  for(coord_idx in seq_len(num_coords)) {
    lh_corr = lh_vertex[mni_array[1, coord_idx], mni_array[2, coord_idx], mni_array[3, coord_idx]];
    rh_corr = rh_vertex[mni_array[1, coord_idx], mni_array[2, coord_idx], mni_array[3, coord_idx]];
    if(lh_corr != 0) { # vertex is from left hemi
      verts[coord_idx] = lh_corr;
      hemi[coord_idx] = 'lh';
      fs_coords[coord_idx, ] = lh_surf$vertices[coord_idx, ];
    } else if(rh_corr != 0) {
      verts[coord_idx] = rh_corr;
      hemi[coord_idx] = 'rh';
      fs_coords[coord_idx, ] = rh_surf$vertices[coord_idx, ];
    } else {
      message(sprintf("Input coord set %d not within MNI cortex mask, returning NaNs.", coord_idx));
      verts[coord_idx] = NaN;
      fs_coords[coord_idx, ] = rep(NaN, 3L);
    }
  }
  return(list("fsaverage_vertices"=verts, "hemi"=hemi, "fsaverage_coords"=fs_coords));
}
