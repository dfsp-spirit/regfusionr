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
#' @examples
#' \dontrun{
#'   mni_ras = c(60.0, 0.0, 10.0)
#'   res = mni152_coords_to_fsaverage(mni_ras, surface = "white");
#'   res$fsaverage_vertices;   # 9092
#' }
#'
#' @export
mni152_coords_to_fsaverage <- function(coords, surface='white', fs_home=Sys.getenv("FS_HOME"), silent = TRUE) {
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
  lh_vertex = freesurferformats::read.fs.mgh(lh_map_file, with_header = FALSE, drop_empty_dims = TRUE); # 256x256x256 array
  rh_vertex = freesurferformats::read.fs.mgh(rh_map_file, with_header = FALSE, drop_empty_dims = TRUE); # 256x256x256 array

  # Load cortex mask volume file.
  cortex_mask_file = get_data_file("FSL_MNI152_FS4.5.0_cortex_estimate.nii.gz", subdir = "coordmap");
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
        message(sprintf("Input coord set %d not within MNI cortex mask, returning NaNs.", coord_idx));
      }
      verts[coord_idx] = NaN;
      fs_coords[coord_idx, ] = rep(NaN, 3L);
    }
  }
  return(list("fsaverage_vertices"=verts, "hemi"=hemi, "fsaverage_coords"=fs_coords, "query_mni_coords"=coords, "query_mni_voxels"=mni_voxels));
}


#' @title Map MNI152 voxels of reference file to fsaverage coords and vertices.
#'
#' @description The voxel indices are specific to the reference volume and thus this function is of limited use in general, but it serves as an example on howto achieve voxel mapping.
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


#' @title Map fsaverage vertex indices to MNI152 coordinates.
#'
#' @inheritParams mni152_coords_to_fsaverage
#'
#' @param vertices integer vector of vertex indices (1-based), the \code{n} fsaverage vertices you want to map.
#'
#' @param hemis vector of character strings, the hemispheres of the query vertices. Each entry in the vector has to be either \code{'lh'} or \code{'rh'}. Length must match length of parameter \code{vertices}.
#'
#' @return matrix of dim \code{n x 3}, the MNI152 coordinates for the query vertices, one row per vertex.
#'
#' @export
fsaverage_vertices_to_mni152_coords <- function(vertices, hemis, surface='white', fs_home=Sys.getenv("FS_HOME"), silent = TRUE) {
  if(! is.integer(vertices)) {
    message("Converting vertices to integer values.");
    vertices = as.integer(vertices);
  }
  if(any(!(hemis %in% c("lh", "rh")))) {
    stop("All entries in parameter 'hemis' must be one of 'lh' or 'rh'.");
  }
  num_vertices = length(vertices);
  if(length(hemis) != num_vertices) {
    stop("Lengths of parameters 'vertices' and 'hemis' must match.");
  }
  if(num_vertices < 1L) {
    stop("Must pass at least one query vertex");
  }

  # Load the mapping data.
  mapping_files = list("lh"=NULL, "rh"=NULL); # Filled below.
  ras = list("lh"=NULL, "rh"=NULL);
  if("lh" %in% hemis) { # Only load files we need for speed.
    mapping_files$lh = get_data_file("lh.avgMapping_allSub_RF_ANTs_MNI152_orig_to_fsaverage.txt", subdir = "mappings");
    ras$lh = t(as.matrix(data.table::fread(mapping_files$lh, nrows = 3, header = FALSE)));
  }
  if("rh" %in% hemis) {
    mapping_files$rh = get_data_file("rh.avgMapping_allSub_RF_ANTs_MNI152_orig_to_fsaverage.txt", subdir = "mappings");
    ras$rh = t(as.matrix(data.table::fread(mapping_files$rh, nrows = 3, header = FALSE)));
  }

  mni152_coords = matrix(rep(NA, num_vertices*3L), nrow=num_vertices, ncol=3L);
  for (vertex_local_idx in seq_along(vertices)) {
    vertex_surface_idx = vertices[vertex_local_idx];
    vertex_hemi = hemis[vertex_local_idx];
    mni152_coords[vertex_local_idx, ] = ras[[vertex_hemi]][vertex_surface_idx, ];
  }
  return(mni152_coords);
}

