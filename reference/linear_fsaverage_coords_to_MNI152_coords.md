# Transform MNI305 coords (FreeSurfer fsaverage surface) to MNI152 coordinates using the linear method.

This uses the 4x4 matrix from the FreeSurfer Coordinate Systems
documentation.

## Usage

``` r
linear_fsaverage_coords_to_MNI152_coords(vertex_coords)
```

## Arguments

- vertex_coords:

  nx3 matrix of coordinates, e.g., typically from fsaverage surface
  vertices.

## Value

nx3 numerical matrix if MNI152 coords.

## Note

This is the opposite of using the Wu et al. approach: a linear
transformation matrix is used. This approach is mainly implemented in
this package to allow you to easily check the difference between the
methods.

## Examples

``` r
  fsaverage_verts = matrix(c(10.0, 20.0, 30.0, 0.0, 0.0, 0.0), ncol=3, byrow=TRUE);
  mni152_coords = linear_fsaverage_coords_to_MNI152_coords(fsaverage_verts);
```
