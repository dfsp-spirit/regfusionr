
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
fsaverage_vertices_to_mni152_coords <- function(vertices, hemis, fs_home=Sys.getenv("FS_HOME"), silent = TRUE) {
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

