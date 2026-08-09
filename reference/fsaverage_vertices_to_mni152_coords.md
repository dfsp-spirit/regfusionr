# Map fsaverage vertex indices to MNI152 coordinates.

Map fsaverage vertex indices to MNI152 coordinates.

## Usage

``` r
fsaverage_vertices_to_mni152_coords(
  vertices,
  hemis,
  fs_home = Sys.getenv("FREESURFER_HOME"),
  simplify = FALSE
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

## Value

matrix of dim `n x 3`, the MNI152 coordinates for the query vertices,
one row per vertex. Also see the 'simplify' parameter.

## See also

Use the more general function
[`fsaverage_vertices_to_vol_coords`](https://dfsp-spirit.github.io/regfusionr/reference/fsaverage_vertices_to_vol_coords.md)
for more options.

## Examples

``` r
if (FALSE) { # \dontrun{
  vertices = c(9092L);
  hemis = c("rh");
  mni_coords = fsaverage_vertices_to_mni152_coords(vertices, hemis);
} # }
```
