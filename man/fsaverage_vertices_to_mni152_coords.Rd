% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fsaverage_to_coords.R
\name{fsaverage_vertices_to_mni152_coords}
\alias{fsaverage_vertices_to_mni152_coords}
\title{Map fsaverage vertex indices to MNI152 coordinates.}
\usage{
fsaverage_vertices_to_mni152_coords(
  vertices,
  hemis,
  fs_home = Sys.getenv("FS_HOME"),
  simplify = FALSE
)
}
\arguments{
\item{vertices}{integer vector of vertex indices (1-based), the \code{n} fsaverage vertices you want to map.}

\item{hemis}{vector of character strings, the hemispheres of the query vertices. Each entry in the vector has to be either \code{'lh'} or \code{'rh'}. Length must match length of parameter \code{vertices}.}

\item{fs_home}{character string, path to the FreeSurfer installation. Alternatively, a hemilist of \code{freesurferformats::fs.surface} instances like \code{surface = list("lh"=mysurflh, "rh"=mysurfrh)}. Used to find the surfaces, at \code{<fs_home>/subjects/fsaverage/surf/<hemi>.<surface>}, where hemi is 'lh' and 'rh'. Can be NULL if 'surface' is a hemilist of fs.surface instances.}

\item{simplify}{logical, whether to return a vector instead of a single-row matrix in case only a single query vertex is given.}
}
\value{
matrix of dim \code{n x 3}, the MNI152 coordinates for the query vertices, one row per vertex. Also see the 'simplify' parameter.
}
\description{
Map fsaverage vertex indices to MNI152 coordinates.
}
\seealso{
Use the more general function \code{\link{fsaverage_vertices_to_vol_coords}} for more options.
}
