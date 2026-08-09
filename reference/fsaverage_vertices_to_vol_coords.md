# Map fsaverage vertex indices to MNI152 or Colin27 volumne coordinates.

Map fsaverage vertex indices to MNI152 or Colin27 volumne coordinates.

## Usage

``` r
fsaverage_vertices_to_vol_coords(
  vertices,
  hemis,
  fs_home = Sys.getenv("FREESURFER_HOME"),
  simplify = FALSE,
  rf_type = "RF_ANTs",
  template_type = "MNI152_orig"
)
```

## Arguments

- vertices:

  integer vector of vertex indices (1-based), the `n` fsaverage vertices
  you want to map.

- hemis:

  vector of character strings, the hemispheres of the query vertices.
  Each entry in the vector has to be either `'lh'` or `'rh'`. Length
  must match length of parameter `vertices`.

- fs_home:

  character string, path to the FreeSurfer installation. Alternatively,
  a hemilist of `freesurferformats::fs.surface` instances like
  `surface = list("lh"=mysurflh, "rh"=mysurfrh)`. Used to find the
  surfaces, at `<fs_home>/subjects/fsaverage/surf/<hemi>.<surface>`,
  where hemi is 'lh' and 'rh'. Can be NULL if 'surface' is a hemilist of
  fs.surface instances.

- simplify:

  logical, whether to return a vector instead of a single-row matrix in
  case only a single query vertex is given.

- rf_type:

  the registration fusion type, one of 'RF_ANTs' or 'RF_M3Z'.

- template_type:

  the space into which to map. One of 'MNI152_orig', 'MNI152_norm',
  'Colin27_orig', 'Colin27_norm'. Note that the 'RF_ANTs' rf_type must
  be used for `_orig` templates, and the 'RF_M3Z' type for `_norm`
  templates.

## Value

matrix of dim `n x 3`, the MNI152 or Colin27 coordinates for the query
vertices, one row per vertex. Also see the 'simplify' parameter.

## Examples

``` r
if (FALSE) { # \dontrun{
  vertices = c(9092L, 100L);
  hemis = c("rh", "lh");
  mni_coords = fsaverage_vertices_to_vol_coords(vertices, hemis,
    rf_type = "RF_ANTs", template_type = "MNI152_orig");
} # }
```
