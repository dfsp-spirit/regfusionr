# Find Colin27 coordinate of fsaverage vertex closest to the given MNI305 coordinate.

Find Colin27 coordinate of fsaverage vertex closest to the given MNI305
coordinate.

## Usage

``` r
mni305_coords_to_colin27_coords(
  coords,
  surface = "orig",
  fs_home = Sys.getenv("FREESURFER_HOME"),
  simplify = FALSE
)
```

## Arguments

- coords:

  nx3 numerical matrix, the MNI305 query coordinates.

- surface:

  a character string defining the fsaverage surface to load (like
  `"white"` or `"orig"`), or a pre-loaded hemilist of surfaces (i.e.,
  `freesurferformats::fs.surface` instances)

- fs_home:

  character string, path of the FreeSurfer directory from which the
  fsaverage surfaces should be loaded. Ignored if `surface` is a
  hemilist (in that case the surfaces have already been loaded).

- simplify:

  logical, whether to return a vector instead of a single-row matrix in
  case only a single query coordinate is given.

## Value

nx3 numerical matrix, the Colin27 coordinates for the vertices closest
to the given MNI305 query coordinates. Depending on the distance to the
closest vertex, this may be way off. Also see the 'simplify' parameter.

## Examples

``` r
if (FALSE) { # \dontrun{
  # Get MNI305 coords of fsaverage vertex 9092 (right hemisphere, orig surface).
  fsaverage_home = file.path(Sys.getenv("FREESURFER_HOME"), "subjects");
  rh_orig = freesurferformats::read.fs.surface(
    file.path(fsaverage_home, "fsaverage", "surf", "rh.orig"));
  mni305_pt = rh_orig$vertices[9092L, ];
  c27_coords = mni305_coords_to_colin27_coords(mni305_pt);
} # }
```
