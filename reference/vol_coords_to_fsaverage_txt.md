# Map volume coords to fsaverage using .txt mapping files (fallback for missing .mgz files).

Uses the per-vertex RAS coordinate .txt mapping files and
nearest-neighbor Euclidean distance to find the closest fsaverage vertex
for each query coordinate. Available for all template_type/rf_type
combinations, but slower than the .mgz-based approach.

## Usage

``` r
vol_coords_to_fsaverage_txt(
  coords,
  lh_surf,
  rh_surf,
  silent,
  rf_type,
  template_type
)
```

## Arguments

- coords:

  nx3 numeric matrix of RAS coordinates in the source volume space.

- lh_surf:

  fs.surface instance for the left hemisphere.

- rh_surf:

  fs.surface instance for the right hemisphere.

- silent:

  logical, whether to suppress messages for coordinates outside cortex.

- rf_type:

  character string, the registration fusion type ('RF_ANTs' or
  'RF_M3Z').

- template_type:

  character string, the template type ('FSL_MNI152' or 'SPM_Colin27').

## Value

named list, same structure as
[`vol_coords_to_fsaverage`](https://dfsp-spirit.github.io/regfusionr/reference/vol_coords_to_fsaverage.md).
