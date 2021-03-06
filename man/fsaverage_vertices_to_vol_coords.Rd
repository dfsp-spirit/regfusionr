% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fsaverage_to_coords.R
\name{fsaverage_vertices_to_vol_coords}
\alias{fsaverage_vertices_to_vol_coords}
\title{Map fsaverage vertex indices to MNI152 or Colin27 volumne coordinates.}
\usage{
fsaverage_vertices_to_vol_coords(
  vertices,
  hemis,
  fs_home = Sys.getenv("FS_HOME"),
  simplify = FALSE,
  rf_type = "RF_ANTs",
  template_type = "MNI152_orig"
)
}
\arguments{
\item{vertices}{integer vector of vertex indices (1-based), the \code{n} fsaverage vertices you want to map.}

\item{hemis}{vector of character strings, the hemispheres of the query vertices. Each entry in the vector has to be either \code{'lh'} or \code{'rh'}. Length must match length of parameter \code{vertices}.}

\item{fs_home}{character string, path to the FreeSurfer installation. Alternatively, a hemilist of \code{freesurferformats::fs.surface} instances like \code{surface = list("lh"=mysurflh, "rh"=mysurfrh)}. Used to find the surfaces, at \code{<fs_home>/subjects/fsaverage/surf/<hemi>.<surface>}, where hemi is 'lh' and 'rh'. Can be NULL if 'surface' is a hemilist of fs.surface instances.}

\item{simplify}{logical, whether to return a vector instead of a single-row matrix in case only a single query vertex is given.}

\item{rf_type}{the registration fusion type, one of 'RF_ANTs' or 'RF_M3Z'.}

\item{template_type}{the space into which to map. One of 'MNI152_orig', 'MNI152_norm', 'Colin27_orig', 'Colin27_norm'. Note that the 'RF_ANTs' rf_type must be used for \code{_orig} templates, and the 'RF_M3Z' type for \code{_norm} templates.}
}
\value{
matrix of dim \code{n x 3}, the MNI152 or Colin27 coordinates for the query vertices, one row per vertex. Also see the 'simplify' parameter.
}
\description{
Map fsaverage vertex indices to MNI152 or Colin27 volumne coordinates.
}
