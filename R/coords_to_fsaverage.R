# Standalone function to covert MNI/Colin coordinates to fsaverage coords.

#' @title Map MNI152 coords to fsaverage coords and vertices.
#'
#' @param coords nx3 numeric matrix, the source RAS coordinates in the input image which must be in MNI152 space. The coords must be within the cortex, otherwise the mapping makes no sense and \code{NaN} values are returned for the respective coords.
#'
#' @param surface character string, the fsaverage surface (brain mesh) to load. Must be a valid FreeSurfer surface name like 'white', 'pial', 'orig, 'inflated'.
#'
#' @param fs_home character string, path to the FreeSurfer installation. Alternatively, a hemilist of \code{freesurferformats::fs.surface} instances like \code{surface = list("lh"=mysurflh, "rh"=mysurfrh)}. Used to find the surfaces, at \code{<fs_home>/subjects/fsaverage/surf/<hemi>.<surface>}, where hemi is 'lh' and 'rh'. Can be NULL if 'surface' is a hemilist of fs.surface instances.
#'
#' @param silent logical, whether to suppress output messages in case of coords outside of cortex.
#'
#' @return named list with entries 'fsaverage_vertices': integer vector of fsaverage surface vertex indices, 'hemi': vector of hemi strings for the vertices, 'fsaverage_coords': nx3 numeric matrix of target coordinates, 'query_mni_coords': copy of input parameter coords, 'query_mni_voxels': the voxel indices at the query RAS coords.
#'
#' @author Tim Schäfer for the R version, Wu Jianxiao and CBIG for the original Matlab version.
#'
#' @note See the more general function \code{\link{vol_coords_to_fsaverage}} for more options.
#'
#' @examples
#' \dontrun{
#'   mni_ras = c(60.0, 0.0, 10.0)
#'   res = mni152_coords_to_fsaverage(mni_ras, surface = "white");
#'   res$fsaverage_vertices;   # 9092
#' }
#'
#' @export
mni152_coords_to_fsaverage <- function(coords, surface='white', fs_home=Sys.getenv("FS_HOME"), silent = TRUE) {
  return(vol_coords_to_fsaverage(coords, surface=surface, fs_home = fs_home, silent = silent, rf_type="RF_ANTs", template_type="FSL_MNI152"));
}


#' @title Map Colin27 coords to fsaverage coords and vertices.
#'
#' @param coords nx3 numeric matrix, the source RAS coordinates in the input image which must be in Colin27 space. The coords must be within the cortex, otherwise the mapping makes no sense and \code{NaN} values are returned for the respective coords.
#'
#' @param surface character string, the fsaverage surface (brain mesh) to load. Must be a valid FreeSurfer surface name like 'white', 'pial', 'orig, 'inflated'.
#'
#' @param fs_home character string, path to the FreeSurfer installation. Alternatively, a hemilist of \code{freesurferformats::fs.surface} instances like \code{surface = list("lh"=mysurflh, "rh"=mysurfrh)}. Used to find the surfaces, at \code{<fs_home>/subjects/fsaverage/surf/<hemi>.<surface>}, where hemi is 'lh' and 'rh'. Can be NULL if 'surface' is a hemilist of fs.surface instances.
#'
#' @param silent logical, whether to suppress output messages in case of coords outside of cortex.
#'
#' @return named list with entries 'fsaverage_vertices': integer vector of fsaverage surface vertex indices, 'hemi': vector of hemi strings for the vertices, 'fsaverage_coords': nx3 numeric matrix of target coordinates, 'query_coords': copy of input parameter coords, 'query_voxels': the voxel indices at the query RAS coords.
#'
#' @author Tim Schäfer for the R version, Wu Jianxiao and CBIG for the original Matlab version.
#'
#' @note See the more general function \code{\link{vol_coords_to_fsaverage}} for more options.
#'
#' @examples
#' \dontrun{
#'   c27_ras = c(60.0, 0.0, 10.0)
#'   res = colin27_coords_to_fsaverage(mni_ras, surface = "white");
#'   res$fsaverage_vertices;
#' }
#'
#' @keywords internal
colin27_coords_to_fsaverage <- function(coords, surface='white', fs_home=Sys.getenv("FS_HOME"), silent = TRUE) {
  res = vol_coords_to_fsaverage(coords, surface=surface, fs_home = fs_home, silent = silent, rf_type="RF_ANTs", template_type="SPM_Colin27");
  # rename fields.
  res$query_coords = res$query_mni_coords;
  res$query_voxels = res$query_mni_voxels;
  res$query_mni_coords = NULL;
  res$query_mni_voxels = NULL;
  return(res);
}


#' @title Map MNI152 or Colin27 volume coords to fsaverage coords and vertices.
#'
#' @inheritParams vol_to_fsaverage
#'
#' @param coords nx3 numeric matrix, the source RAS coordinates in the input image which must be in MNI152/Colin27 space. The coords must be within the cortex, otherwise the mapping makes no sense and \code{NaN} values are returned for the respective coords.
#'
#' @param surface character string, the fsaverage surface (brain mesh) to load. Must be a valid FreeSurfer surface name like 'white', 'pial', 'orig, 'inflated'.
#'
#' @param fs_home character string, path to the FreeSurfer installation. Alternatively, a hemilist of \code{freesurferformats::fs.surface} instances like \code{surface = list("lh"=mysurflh, "rh"=mysurfrh)}. Used to find the surfaces, at \code{<fs_home>/subjects/fsaverage/surf/<hemi>.<surface>}, where hemi is 'lh' and 'rh'. Can be NULL if 'surface' is a hemilist of fs.surface instances.
#'
#' @param silent logical, whether to suppress output messages in case of coords outside of cortex.
#'
#' @return named list with entries 'fsaverage_vertices': integer vector of fsaverage surface vertex indices, 'hemi': vector of hemi strings for the vertices, 'fsaverage_coords': nx3 numeric matrix of target coordinates, 'query_mni_coords': copy of input parameter coords, 'query_mni_voxels': the voxel indices at the query RAS coords.
#'
#' @author Tim Schäfer for the R version, Wu Jianxiao and CBIG for the original Matlab version.
#'
#' @examples
#' \dontrun{
#'   mni_ras = c(60.0, 0.0, 10.0)
#'   res = vol_coords_to_fsaverage(mni_ras, surface = "white");
#'   res$fsaverage_vertices;   # 9092
#' }
#'
#' @export
vol_coords_to_fsaverage <- function(coords, surface='white', fs_home=Sys.getenv("FS_HOME"), silent = TRUE, rf_type="RF_ANTs", template_type="FSL_MNI152") {
  if(is.vector(coords)) {
    if(length(coords) %% 3L == 0L) {
      coords = matrix(coords, ncol = 3L);
    }
  }
  check_coords(coords);

  if(! (rf_type %in% c('RF_ANTs', 'RF_M3Z'))) {
    stop(sprintf("Parameter 'rf_type' must be one of 'RF_ANTs' or 'RF_M3Z', but is '%s'.\n", rf_type));
  }

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
  lh_map_file = get_data_file(sprintf("%s_FS4.5.0_%s_avgMapping.vertex.lh.mgz", template_type, rf_type), subdir = "coordmap");
  rh_map_file = get_data_file(sprintf("%s_FS4.5.0_%s_avgMapping.vertex.rh.mgz", template_type, rf_type), subdir = "coordmap");
  lh_vertex = freesurferformats::read.fs.mgh(lh_map_file, with_header = FALSE, drop_empty_dims = TRUE); # 256x256x256 array
  rh_vertex = freesurferformats::read.fs.mgh(rh_map_file, with_header = FALSE, drop_empty_dims = TRUE); # 256x256x256 array

  # Load cortex mask volume file.
  cortex_mask_file = get_data_file(sprintf("%s_FS4.5.0_cortex_estimate.nii.gz", template_type), subdir = "coordmap");
  # Note: cortex mask is an fs.mgh instance, the cortex_mask$data is a 256x256x256 matrix due to drop_empty_dims (256x256x256x1 originally).
  cortex_mask = freesurferformats::read.fs.volume(cortex_mask_file, with_header = TRUE, drop_empty_dims = TRUE);

  # Do the masking in the 3D arrays.
  lh_vertex[which(cortex_mask$data == 0)] = 0;
  rh_vertex[which(cortex_mask$data == 0)] = 0;

  # Convert input RAS coords to voxel indices (IJK) for the matrix.
  mni_voxels = freesurferformats::doapply.transform.mtx(coords, freesurferformats::mghheader.ras2vox(cortex_mask)) + 1L;
  if(is.vector(mni_voxels)) {
    mni_voxels = matrix(mni_voxels, ncol = 3, byrow = TRUE);
  }
  # Reorder columns.
  mni_array = matrix(data = c(mni_voxels[,2], mni_voxels[,1], mni_voxels[,3]) , ncol = 3, byrow = FALSE);

  num_coords = nrow(coords);
  verts = rep(0L, num_coords);
  fs_coords = matrix(rep(0.0, (num_coords * 3L)), ncol = 3L);
  hemi = rep(NA, num_coords);
  for(coord_idx in seq_len(num_coords)) {
    lh_corr = lh_vertex[mni_array[coord_idx, 1], mni_array[coord_idx, 2], mni_array[coord_idx, 3]];
    rh_corr = rh_vertex[mni_array[coord_idx, 1], mni_array[coord_idx, 2], mni_array[coord_idx, 3]];
    if(lh_corr != 0) { # vertex is from left hemi
      verts[coord_idx] = lh_corr;
      hemi[coord_idx] = 'lh';
      fs_coords[coord_idx, ] = lh_surf$vertices[coord_idx, ];
    } else if(rh_corr != 0) {
      verts[coord_idx] = rh_corr;
      hemi[coord_idx] = 'rh';
      fs_coords[coord_idx, ] = rh_surf$vertices[coord_idx, ];
    } else {
      if(! silent) {
        message(sprintf("Input coord set %d not within volume cortex mask, returning NaNs.", coord_idx));
      }
      verts[coord_idx] = NaN;
      fs_coords[coord_idx, ] = rep(NaN, 3L);
    }
  }
  return(list("fsaverage_vertices"=verts, "hemi"=hemi, "fsaverage_coords"=fs_coords, "query_mni_coords"=coords, "query_mni_voxels"=mni_voxels));
}


#' @title Map MNI152 voxels of reference file to fsaverage coords and vertices.
#'
#' @description The voxel indices are specific to the reference volume and thus this function is of limited use in general, but it serves as an example on how to achieve voxel mapping.
#'
#' @param voxels integer nx3 matrix, the IJK voxel indices in the 256x256x256 cortex mask reference file, 1-based.
#'
#' @inheritParams mni152_coords_to_fsaverage
#'
#' @return see \code{mni152_coords_to_fsaverage}
#'
#' @note The voxel indices used by this function are specific to the cortex mask reference volume, so it is preferable to use the \code{} function. However, the code of this function illustrates how to get the fsaverage coords based on the IJK voxel indices for your own image, so take it as a demo.
#'
#' @keywords internal
mni152_voxels_to_fsaverage <- function(voxels, surface='white', fs_home=Sys.getenv("FS_HOME"), silent = TRUE) {
  if(is.vector(voxels)) {
    if(length(voxels) %% 3L == 0L) {
      voxels = matrix(voxels, ncol = 3L);
    }
  }
  check_voxels(voxels);
  voxels = voxels - 1L; # required 0-based for freesurferformats::doapply.transform.mtx() below.

  cortex_mask_file = get_data_file("FSL_MNI152_FS4.5.0_cortex_estimate.nii.gz", subdir = "coordmap");
  # Note: cortex mask is an fs.mgh instance, the cortex_mask$data is a 256x256x256 matrix due to drop_empty_dims (256x256x256x1 originally).
  cortex_mask = freesurferformats::read.fs.volume(cortex_mask_file, with_header = TRUE, drop_empty_dims = TRUE);

  # Convert input RAS coords to voxel indices (IJK) for the matrix.
  mni_coords = freesurferformats::doapply.transform.mtx(voxels, freesurferformats::mghheader.vox2ras(cortex_mask));
  return(mni152_coords_to_fsaverage(mni_coords, surface = surface, fs_home = fs_home, silent = silent));
}

#' @title Map Colin27 voxels of reference file to fsaverage coords and vertices.
#'
#' @description The voxel indices are specific to the reference volume and thus this function is of limited use in general, but it serves as an example on how to achieve voxel mapping.
#'
#' @param voxels integer nx3 matrix, the IJK voxel indices in the 256x256x256 cortex mask reference file, 1-based.
#'
#' @inheritParams mni152_voxels_to_fsaverage
#'
#' @return see \code{mni152_coords_to_fsaverage}
#'
#' @note The voxel indices used by this function are specific to the cortex mask reference volume, so it is preferable to use the \code{} function. However, the code of this function illustrates how to get the fsaverage coords based on the IJK voxel indices for your own image, so take it as a demo.
#'
#' @keywords internal
colin27_voxels_to_fsaverage <- function(voxels, surface='white', fs_home=Sys.getenv("FS_HOME"), silent = TRUE) {
  if(is.vector(voxels)) {
    if(length(voxels) %% 3L == 0L) {
      voxels = matrix(voxels, ncol = 3L);
    }
  }
  check_voxels(voxels);
  voxels = voxels - 1L; # required 0-based for freesurferformats::doapply.transform.mtx() below.

  cortex_mask_file = get_data_file("SPM_Colin27_FS4.5.0_cortex_estimate.nii.gz", subdir = "coordmap");
  # Note: cortex mask is an fs.mgh instance, the cortex_mask$data is a 256x256x256 matrix due to drop_empty_dims (256x256x256x1 originally).
  cortex_mask = freesurferformats::read.fs.volume(cortex_mask_file, with_header = TRUE, drop_empty_dims = TRUE);

  # Convert input RAS coords to voxel indices (IJK) for the matrix.
  mni_coords = freesurferformats::doapply.transform.mtx(voxels, freesurferformats::mghheader.vox2ras(cortex_mask));
  return(colin27_coords_to_fsaverage(mni_coords, surface = surface, fs_home = fs_home, silent = silent));
}


