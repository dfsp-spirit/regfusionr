# Map Colin27 voxels of reference file to fsaverage coords and vertices.

The voxel indices are specific to the reference volume and thus this
function is of limited use in general, but it serves as an example on
how to achieve voxel mapping.

## Usage

``` r
colin27_voxels_to_fsaverage(
  voxels,
  surface = "white",
  fs_home = Sys.getenv("FREESURFER_HOME"),
  silent = TRUE
)
```

## Arguments

- voxels:

  integer nx3 matrix, the IJK voxel indices in the 256x256x256 cortex
  mask reference file, 1-based.

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

see `mni152_coords_to_fsaverage`

## Note

The voxel indices used by this function are specific to the cortex mask
reference volume, so it is preferable to use the function. However, the
code of this function illustrates how to get the fsaverage coords based
on the IJK voxel indices for your own image, so take it as a demo.
