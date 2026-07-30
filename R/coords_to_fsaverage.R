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
#' @export
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

  # Check whether vertex .mgz mapping files exist. If not, fall back to .txt-based approach.
  lh_map_path = system.file("extdata", "coordmap", sprintf("%s_FS4.5.0_%s_avgMapping.vertex.lh.mgz", template_type, rf_type), package = "regfusionr");
  if(lh_map_path == "") {
    if(! silent) {
      message(sprintf("Vertex map .mgz file not found for template_type='%s', rf_type='%s'. Using .txt-based nearest-vertex mapping instead.", template_type, rf_type));
    }
    return(vol_coords_to_fsaverage_txt(coords, lh_surf, rh_surf, silent, rf_type, template_type));
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
      fs_coords[coord_idx, ] = lh_surf$vertices[verts[coord_idx], ];
    } else if(rh_corr != 0) {
      verts[coord_idx] = rh_corr;
      hemi[coord_idx] = 'rh';
      fs_coords[coord_idx, ] = rh_surf$vertices[verts[coord_idx], ];
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


#' @title Map volume coords to fsaverage using .txt mapping files (fallback for missing .mgz files).
#'
#' @description Uses the per-vertex RAS coordinate .txt mapping files and nearest-neighbor Euclidean distance to find the closest fsaverage vertex for each query coordinate. Available for all template_type/rf_type combinations, but slower than the .mgz-based approach.
#'
#' @param coords nx3 numeric matrix of RAS coordinates in the source volume space.
#'
#' @param lh_surf fs.surface instance for the left hemisphere.
#'
#' @param rh_surf fs.surface instance for the right hemisphere.
#'
#' @param silent logical, whether to suppress messages for coordinates outside cortex.
#'
#' @param rf_type character string, the registration fusion type ('RF_ANTs' or 'RF_M3Z').
#'
#' @param template_type character string, the template type ('FSL_MNI152' or 'SPM_Colin27').
#'
#' @return named list, same structure as \code{\link{vol_coords_to_fsaverage}}.
#'
#' @keywords internal
vol_coords_to_fsaverage_txt <- function(coords, lh_surf, rh_surf, silent, rf_type, template_type) {

  # Map the template_type (FSL_MNI152/SPM_Colin27) + rf_type to the mapping file naming convention.
  if(template_type == "FSL_MNI152") {
    if(rf_type == "RF_ANTs") {
      mapping_template = "MNI152_orig";
    } else {
      mapping_template = "MNI152_norm";
    }
  } else if(template_type == "SPM_Colin27") {
    if(rf_type == "RF_ANTs") {
      mapping_template = "Colin27_orig";
    } else {
      mapping_template = "Colin27_norm";
    }
  } else {
    stop(sprintf("Invalid template_type '%s'.", template_type));
  }

  # Load the per-vertex RAS coordinate mapping files.
  # These files have 3 rows (x, y, z) and 163842 columns (one per fsaverage vertex).
  # Transpose to get 163842 x 3 matrices.
  lh_mapping_file = get_data_file(sprintf("lh.avgMapping_allSub_%s_%s_to_fsaverage.txt", rf_type, mapping_template), subdir = "mappings");
  rh_mapping_file = get_data_file(sprintf("rh.avgMapping_allSub_%s_%s_to_fsaverage.txt", rf_type, mapping_template), subdir = "mappings");

  ras_lh = t(as.matrix(data.table::fread(lh_mapping_file, nrows = 3, header = FALSE)));
  ras_rh = t(as.matrix(data.table::fread(rh_mapping_file, nrows = 3, header = FALSE)));

  num_coords = nrow(coords);
  num_verts_lh = nrow(ras_lh);
  num_verts_rh = nrow(ras_rh);

  verts = rep(NaN, num_coords);
  fs_coords = matrix(rep(NaN, (num_coords * 3L)), ncol = 3L);
  hemi = rep(NA_character_, num_coords);

  # Distance threshold for "outside cortex" detection (in mm).
  # A query coordinate more than this distance from any cortical vertex is considered outside cortex.
  dist_threshold = 15.0;

  for(coord_idx in seq_len(num_coords)) {
    query_coord = coords[coord_idx, ];

    # Compute Euclidean distances to all fsaverage vertices in each hemisphere.
    lh_dists = sqrt(rowSums((ras_lh - matrix(query_coord, nrow = num_verts_lh, ncol = 3L, byrow = TRUE))^2));
    rh_dists = sqrt(rowSums((ras_rh - matrix(query_coord, nrow = num_verts_rh, ncol = 3L, byrow = TRUE))^2));

    lh_min_dist = min(lh_dists);
    rh_min_dist = min(rh_dists);

    if(lh_min_dist <= rh_min_dist && lh_min_dist < dist_threshold) {
      verts[coord_idx] = which.min(lh_dists);
      hemi[coord_idx] = 'lh';
      fs_coords[coord_idx, ] = lh_surf$vertices[verts[coord_idx], ];
    } else if(rh_min_dist < dist_threshold) {
      verts[coord_idx] = which.min(rh_dists);
      hemi[coord_idx] = 'rh';
      fs_coords[coord_idx, ] = rh_surf$vertices[verts[coord_idx], ];
    } else {
      if(! silent) {
        message(sprintf("Input coord set %d too far from any cortical vertex (min dist: %.2f mm), returning NaNs.", coord_idx, min(lh_min_dist, rh_min_dist)));
      }
    }
  }

  # Compute voxel indices for the query coords using the cortex mask (for consistency with the .mgz code path).
  cortex_mask_file = get_data_file(sprintf("%s_FS4.5.0_cortex_estimate.nii.gz", template_type), subdir = "coordmap");
  cortex_mask = freesurferformats::read.fs.volume(cortex_mask_file, with_header = TRUE, drop_empty_dims = TRUE);
  mni_voxels = freesurferformats::doapply.transform.mtx(coords, freesurferformats::mghheader.ras2vox(cortex_mask)) + 1L;
  if(is.vector(mni_voxels)) {
    mni_voxels = matrix(mni_voxels, ncol = 3, byrow = TRUE);
  }

  return(list("fsaverage_vertices" = verts, "hemi" = hemi, "fsaverage_coords" = fs_coords, "query_mni_coords" = coords, "query_mni_voxels" = mni_voxels));
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


#' @title Map MNI152 coordinates to Colin27 coordinates via fsaverage.
#'
#' @description Convenience function that maps MNI152 RAS coordinates to Colin27 RAS coordinates by chaining through fsaverage: MNI152 coords are first mapped to fsaverage vertices using \code{\link{mni152_coords_to_fsaverage}}, then the fsaverage vertices are mapped to Colin27 coordinates using \code{\link{fsaverage_vertices_to_colin27_coords}}.
#'
#' @inheritParams mni152_coords_to_fsaverage
#'
#' @param simplify logical, whether to return a vector instead of a single-row matrix when only a single query coordinate is given.
#'
#' @return nx3 numeric matrix of Colin27 RAS coordinates. Also see the 'simplify' parameter.
#'
#' @export
mni152_coords_to_colin27_coords <- function(coords, surface='white', fs_home=Sys.getenv("FS_HOME"), silent = TRUE, simplify = FALSE) {
  res = mni152_coords_to_fsaverage(coords, surface = surface, fs_home = fs_home, silent = silent);
  valid_mask = (! is.na(res$fsaverage_vertices)) & (! is.nan(res$fsaverage_vertices));
  c27_coords = matrix(rep(NA_real_, nrow(coords) * 3L), ncol = 3L);
  if(any(valid_mask)) {
    c27_coords[valid_mask, ] = fsaverage_vertices_to_colin27_coords(
      res$fsaverage_vertices[valid_mask], res$hemi[valid_mask],
      fs_home = fs_home, simplify = FALSE);
  }
  if(nrow(coords) == 1L & simplify) {
    return(c27_coords[1L, ]);
  }
  return(c27_coords);
}


#' @title Map Colin27 coordinates to MNI152 coordinates via fsaverage.
#'
#' @description Convenience function that maps Colin27 RAS coordinates to MNI152 RAS coordinates by chaining through fsaverage: Colin27 coords are first mapped to fsaverage vertices using \code{\link{colin27_coords_to_fsaverage}}, then the fsaverage vertices are mapped to MNI152 coordinates using \code{\link{fsaverage_vertices_to_mni152_coords}}.
#'
#' @inheritParams mni152_coords_to_fsaverage
#'
#' @param simplify logical, whether to return a vector instead of a single-row matrix when only a single query coordinate is given.
#'
#' @return nx3 numeric matrix of MNI152 RAS coordinates. Also see the 'simplify' parameter.
#'
#' @export
colin27_coords_to_mni152_coords <- function(coords, surface='white', fs_home=Sys.getenv("FS_HOME"), silent = TRUE, simplify = FALSE) {
  res = colin27_coords_to_fsaverage(coords, surface = surface, fs_home = fs_home, silent = silent);
  valid_mask = (! is.na(res$fsaverage_vertices)) & (! is.nan(res$fsaverage_vertices));
  mni_coords = matrix(rep(NA_real_, nrow(coords) * 3L), ncol = 3L);
  if(any(valid_mask)) {
    mni_coords[valid_mask, ] = fsaverage_vertices_to_mni152_coords(
      res$fsaverage_vertices[valid_mask], res$hemi[valid_mask],
      fs_home = fs_home, simplify = FALSE);
  }
  if(nrow(coords) == 1L & simplify) {
    return(mni_coords[1L, ]);
  }
  return(mni_coords);
}


