# Map Colin27 coordinates to MNI152 coordinates via fsaverage.

Convenience function that maps Colin27 RAS coordinates to MNI152 RAS
coordinates by chaining through fsaverage: Colin27 coords are first
mapped to fsaverage vertices using
[`colin27_coords_to_fsaverage`](https://dfsp-spirit.github.io/regfusionr/reference/colin27_coords_to_fsaverage.md),
then the fsaverage vertices are mapped to MNI152 coordinates using
[`fsaverage_vertices_to_mni152_coords`](https://dfsp-spirit.github.io/regfusionr/reference/fsaverage_vertices_to_mni152_coords.md).

## Usage

``` r
colin27_coords_to_mni152_coords(
  coords,
  surface = "white",
  fs_home = Sys.getenv("FREESURFER_HOME"),
  silent = TRUE,
  simplify = FALSE
)
```

## Arguments

- coords:

  nx3 numeric matrix, the source RAS coordinates in the input image
  which must be in MNI152 space. The coords must be within the cortex,
  otherwise the mapping makes no sense and `NaN` values are returned for
  the respective coords.

- surface:

  character string, the fsaverage surface (brain mesh) to load. Must be
  a valid FreeSurfer surface name like 'white', 'pial', 'orig,
  'inflated'.

- fs_home:

  character string, path to the FreeSurfer installation. Alternatively,
  a hemilist of `freesurferformats::fs.surface` instances like
  `surface = list("lh"=mysurflh, "rh"=mysurfrh)`. Used to find the
  surfaces, at `<fs_home>/subjects/fsaverage/surf/<hemi>.<surface>`,
  where hemi is 'lh' and 'rh'. Can be NULL if 'surface' is a hemilist of
  fs.surface instances.

- silent:

  logical, whether to suppress output messages in case of coords outside
  of cortex.

- simplify:

  logical, whether to return a vector instead of a single-row matrix
  when only a single query coordinate is given.

## Value

nx3 numeric matrix of MNI152 RAS coordinates. Also see the 'simplify'
parameter.

## Examples

``` r
if (FALSE) { # \dontrun{
  c27_ras = c(60.0, 0.0, 10.0);
  mni_coords = colin27_coords_to_mni152_coords(c27_ras, surface = "white");
} # }
```
