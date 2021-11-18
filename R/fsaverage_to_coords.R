
#' @title Map fsaverage vertex indices to MNI152 coordinates.
#'
#' @inheritParams mni152_coords_to_fsaverage
#'
#' @param vertices integer vector of vertex indices (1-based), the \code{n} fsaverage vertices you want to map.
#'
#' @param hemis vector of character strings, the hemispheres of the query vertices. Each entry in the vector has to be either \code{'lh'} or \code{'rh'}. Length must match length of parameter \code{vertices}.
#'
#' @param simplify logical, whether to return a vector instead of a single-row matrix in case only a single query vertex is given.
#'
#' @return matrix of dim \code{n x 3}, the MNI152 coordinates for the query vertices, one row per vertex. Also see the 'simplify' parameter.
#'
#' @export
fsaverage_vertices_to_mni152_coords <- function(vertices, hemis, fs_home=Sys.getenv("FS_HOME"), simplify = FALSE) {
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
  if(length(vertices) == 1L & simplify) {
    return(mni152_coords[1L, ]);
  } else {
    return(mni152_coords);
  }
}


#' @title Find MNI152 coordinate of fsaverage vertex closest to the given MNI305 coordinate.
#'
#' @inheritParams fsaverage_vertices_to_mni152_coords
#'
#' @param coords nx3 numerical matrix, the MNI305 query coordinates.
#'
#' @param surface a character string defining the fsaverage surface to load (like \code{"white"} or \code{"orig"}), or a pre-loaded hemilist of surfaces (i.e., \code{freesurferformats::fs.surface} instances)
#'
#' @param fs_home character string, path of the FreeSurfer directory from which the fsaverage surfaces should be loaded. Ignored if \code{surface} is a hemilist (in that case the surfaces have already been loaded).
#'
#' @return nx3 numerical matrix, the MNI152 coordinates for the vertices closest to the given MNI305 query coordinates. Depending on the distance to the closest vertex, this may be way off. Also see the 'simplify' parameter.
#'
#' @export
mni305_coords_to_mni152_coords <- function(coords, surface = "orig", fs_home=Sys.getenv("FS_HOME"), simplify = FALSE) {
  if(! is.list(surface)) {
    surface_name = surface;
    lh_surf = freesurferformats::read.fs.surface(file.path(fs_home, "subjects", "fsaverage", "surf", sprintf("lh.%s", surface_name)));
    rh_surf = freesurferformats::read.fs.surface(file.path(fs_home, "subjects", "fsaverage", "surf", sprintf("rh.%s", surface_name)));
    surface = list("lh" = lh_surf, "rh" = rh_surf);
  }
  dist_info = coord_closest_vertex(coords, surface);
  return(fsaverage_vertices_to_mni152_coords(dist_info$both_closest_vertex, dist_info$both_hemi, fs_home=fs_home, simplify=simplify));
}


#' @title Find closest surface vertex to a point, using Euclidean distance.
#'
#' @param coordinate \code{nx3} numerical matrix or vector of length 3, the query point coordinates.
#'
#' @param surfaces hemilist of \code{fs.surface} instances
#'
#' @return a data.frame with columns named 'query_x', 'query_y', 'query_z', 'lh_closest_vertex', 'lh_distance', 'rh_closest_vertex', 'rh_distance', 'both_closest_vertex', 'both_distance', 'both_hemi'.
#'
#' @keywords internal
coord_closest_vertex <- function(coordinate, surfaces) {
  if(is.vector(coordinate) & length(coordinate) == 3L) {
    coordinate = matrix(coordinate, ncol = 3, nrow = 1, byrow = TRUE);
  }
  if(! is.matrix(coordinate)) {
    stop("Parameter 'coordinate' must be a vector of length 3 or an nx3 matrix.");
  }

  num_coords = nrow(coordinate);
  lh_closest_vertex = rep(NA, num_coords);
  rh_closest_vertex = rep(NA, num_coords);
  both_closest_vertex = rep(NA, num_coords);
  lh_distance = rep(NA, num_coords);
  rh_distance = rep(NA, num_coords);
  both_distance = rep(NA, num_coords);
  both_hemi = rep(NA, num_coords);

  for (row_idx in seq.int(num_coords)) {
    has_surf = FALSE;
    if(freesurferformats::is.fs.surface(surfaces$lh)) {
      lh_vd = freesurferformats::vertexdists.to.point(surfaces$lh, coordinate[row_idx, ]);
      lh_closest_vertex[row_idx] = which.min(lh_vd);
      lh_distance[row_idx] = lh_vd[lh_closest_vertex[row_idx]];
      both_closest_vertex[row_idx] = lh_closest_vertex[row_idx]; # for now, may change below.
      both_distance[row_idx] = lh_distance[row_idx];             # for now, may change below.
      both_hemi[row_idx] = "lh";
      has_surf = TRUE;
    }
    if(freesurferformats::is.fs.surface(surfaces$rh)) {
      rh_vd = freesurferformats::vertexdists.to.point(surfaces$rh, coordinate[row_idx, ]);
      rh_closest_vertex[row_idx] = which.min(rh_vd);
      rh_distance[row_idx] = rh_vd[rh_closest_vertex[row_idx]];
      if(has_surf) { # whether there was a left surface. in this case we need to compare the values and pick the smaller dist.
        if(lh_distance[row_idx] < rh_distance[row_idx]) {
          both_closest_vertex[row_idx] = lh_closest_vertex[row_idx];
          both_distance[row_idx] = lh_distance[row_idx];
          both_hemi[row_idx] = "lh";
        } else {
          both_closest_vertex[row_idx] = rh_closest_vertex[row_idx];
          both_distance[row_idx] = rh_distance[row_idx];
          both_hemi[row_idx] = "rh";
        }
      } else { # Only the rh surface was given.
        both_closest_vertex[row_idx] = rh_closest_vertex[row_idx];
        both_distance[row_idx] = rh_distance[row_idx];
        both_hemi[row_idx] = "rh";
      }
      has_surf = TRUE;
    }
    if(! has_surf) {
      stop("The hemilist in parameter 'surfaces' must contain at least one fs.surface instance in keys 'lh' or 'rh'.");
    }
  }
  return(data.frame("query_x"=coordinate[,1], "query_y"=coordinate[,2], "query_z"=coordinate[,3], "lh_closest_vertex"=lh_closest_vertex, "lh_distance"=lh_distance, "rh_closest_vertex"=rh_closest_vertex, "rh_distance"=rh_distance, "both_closest_vertex"=both_closest_vertex, "both_distance"=both_distance, "both_hemi"=both_hemi, stringsAsFactors = FALSE));
}

