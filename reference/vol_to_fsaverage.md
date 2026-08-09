# Project or map values from MNI volume to fsaverage surface.

Applies the Wu et al. regfusion method to obtain surface coords, then
interpolates values.

## Usage

``` r
vol_to_fsaverage(
  input_img,
  template_type,
  rf_type = "RF_ANTs",
  interp = "linear",
  out_type = "curv",
  out_dir = "."
)
```

## Arguments

- input_img:

  3D or 4D NIFTI or MGZ image instance of type `fs.volume`, or a
  character string that will be interpreted as a file system path to
  such a volume that should be loaded with
  [`freesurferformats::read.fs.volume`](https://rdrr.io/pkg/freesurferformats/man/read.fs.volume.html).
  If 4D, the 4th dimension is considered the time/subject dimension.

- template_type:

  character string, the source template or the space that your input
  image is in. One of 'MNI152_orig', 'Colin27_orig', 'MNI152_norm',
  'Colin27_norm'.

- rf_type:

  the regfusion type to use, one of 'RF_ANTs' or 'RF_M3Z'.

- interp:

  interpolation method, currently only 'linear' is supported.

- out_type:

  character string, the format of the output files. One of the
  following: 'curv' for FreeSurfer curv format, 'mgz' for FreeSurfer MGZ
  format, 'gii' for GIFTI format.

- out_dir:

  character string, the path to a writable output directory. If `NULL`,
  the returned named list contains the projected data (instead of the
  path of the file it was written to) at the keys 'lh' and 'rh', and the
  parameter 'out_type' is ignored.

## Value

named list of 2 character strings, the output files (for the 2
hemispheres) at keys 'lh' and 'rh'. See the documentation for parameter
'out_dir' if you want the data in R instead.

## Author

Tim Schäfer for the R version, Wu Jianxiao and CBIG for the original
Matlab version.

## Examples

``` r
if (FALSE) { # \dontrun{
  # Project a 3D MNI152 volume to fsaverage surface space, get data in R.
  input_file = system.file("extdata", "testdata",
    "MNI_probMap_ants.central_sulc.nii.gz", package = "regfusionr");
  res = vol_to_fsaverage(input_file, template_type = 'MNI152_orig',
    rf_type = 'RF_ANTs', out_dir = NULL);
  # res$lh contains per-vertex data for the left hemisphere.
} # }
```
