# Perform coord transform and trilinear data interpolation.

Perform coord transform and trilinear data interpolation.

## Usage

``` r
project_data(data, affine, ras, interp = "linear")
```

## Arguments

- data:

  3D or 3D source image values

- affine:

  4x4 double matrix, the ras2vox transformation matrix

- ras:

  nx3 double matrix, the source RAS coordinates

- interp:

  character string, the interpolation mode. Only 'linear' is currently
  supported. We should support 'nearest' as well.

## Value

a numerical vector for 3D input data or a numerical 4D array for 4D
input data, the data values interpolated at the RAS coordinates. For the
4D input case, the 2D output is encoded in a 4D array with two
dimensions of length 1 as follows: The first dimension contains the
projected and interpolated per-vertex data values for a single frame,
dimensions 2 and 3 are both of length 1, and the fourth dimensions is
the timepoint/subject dimension (the frames).
