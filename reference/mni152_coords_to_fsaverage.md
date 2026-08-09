# Map MNI152 coords to fsaverage coords and vertices.

Map MNI152 coords to fsaverage coords and vertices.

## Usage

``` r
mni152_coords_to_fsaverage(
  coords,
  surface = "white",
  fs_home = Sys.getenv("FREESURFER_HOME"),
  silent = TRUE
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

## Value

named list with entries 'fsaverage_vertices': integer vector of
fsaverage surface vertex indices, 'hemi': vector of hemi strings for the
vertices, 'fsaverage_coords': nx3 numeric matrix of target coordinates,
'query_mni_coords': copy of input parameter coords, 'query_mni_voxels':
the voxel indices at the query RAS coords.

## Note

See the more general function
[`vol_coords_to_fsaverage`](https://dfsp-spirit.github.io/regfusionr/reference/vol_coords_to_fsaverage.md)
for more options.

## Author

Tim Schäfer for the R version, Wu Jianxiao and CBIG for the original
Matlab version.

## Examples

``` r
if (FALSE) { # \dontrun{
  mni_ras = c(60.0, 0.0, 10.0)
  res = mni152_coords_to_fsaverage(mni_ras, surface = "white");
  res$fsaverage_vertices;   # 9092
} # }
```
