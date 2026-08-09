# Project or map per-vertex values from the fsaverage surface to the cortex voxels of an MNI volume.

Applies the Wu et al. regfusion method to obtain MNI volume coordinates,
then interpolates values.

## Usage

``` r
fsaverage_to_vol(
  lh_input,
  rh_input,
  target_space = "FSL_MNI152",
  rf_type = "RF_ANTs",
  interp = "linear",
  out_type = "mgz",
  out_dir = NULL,
  fsaverage_path = NULL,
  silent = TRUE
)
```

## Arguments

- lh_input:

  numerical vector of per-vertex data for the left hemisphere of the
  template subject. Must contain 163842 values for the `fsaverage`
  template. Input for fsaverage6 (40962 values) or fsaverage5 (10242
  values) can also be used and will be upsampled using
  [`nn_interpolate_kdtree`](https://dfsp-spirit.github.io/haze/reference/nn_interpolate_kdtree.html).
  Automatic upsamping is only supported with `interp='linear'`.

- rh_input:

  numerical vector of per-vertex data for the right hemisphere of the
  template subject. Must contain 163842 values for the `fsaverage`
  template. See description for `lh_input` for more details.

- target_space:

  character string, the target template or the space that your output
  volume should be in. One of 'FSL_MNI152' or 'SPM_Colin27'.

- rf_type:

  the `regfusion` type to use, one of 'RF_ANTs' or 'RF_M3Z'.

- interp:

  interpolation method, currently only 'linear' and 'nearest' are
  supported. The performance of the 'linear' method is currently quite
  bad, and it will be rewritten in `C++` when I find the time.

- out_type:

  character string, the format of the output files. One of the
  following: 'mgz' or 'mgh' for FreeSurfer MGZ/MGH format, 'nii' for
  NIFTI v1 format. Ignored unless out_dir is not NULL.

- out_dir:

  optional character string, the path to a writable output directory to
  which the output should be written as volume files. If `NULL`, no data
  is written to files. If out_dir is not `NULL`, the return value
  additionally contains the following keys: 'out_file' and the output
  file format at key 'out_format', (and 'out_file_seg'/'out_format_seg'
  for the respective `seg` versions).

- fsaverage_path:

  character string or NULL, the file system path to the `fsaverage`
  directory (NOT including the 'fsaverage' dir itself). If `NULL`,
  defaults to the return value of
  [`fsbrain::fsaverage.path()`](https://rdrr.io/pkg/fsbrain/man/fsaverage.path.html)
  on the system. This path is used to read the spherical surface (both
  hemisphere meshes) of the template subject.

- silent:

  logical, whether to suppress status messages.

## Value

named list with keys 'projected' and 'projected_seg', each of which
holds an `fs.volume` instance, its 'data' key holds a `256x256x256`
array with the projected data. The data in 'projected_seg' is identical
to the data in 'projected', with the exception that data values
originating from the right hemisphere have been incremented by 1000. See
`out_dir` parameter to easily write results to files. If out_dir is not
`NULL`, the return value additionally contains the following keys:
'out_file' and the output file format at key 'out_format'.

## Note

This function requires the packages 'fsbrain' and 'haze', which are
optional dependencies. The package 'fsbrain' can be installed from CRAN.
For 'haze', see `https://github.com/dfsp-spirit/haze`.

This function is quite expensive computationally, especially when using
`interp = 'linear'`.

## Author

Tim Schäfer for the R version, Wu Jianxiao and CBIG for the original
Matlab version.

## Examples

``` r
if (FALSE) { # \dontrun{
  lh_input = rnorm(163842L, 3.0, 0.2);
  rh_input = rnorm(163842L, 3.0, 0.2);
  res = fsaverage_to_vol(lh_input, rh_input, "FSL_MNI152");
} # }
```
