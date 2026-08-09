# Get fsaverage (MNI305) to MNI152 transformation matrix.

This returns the 4x4 matrix from the FreeSurfer Coordinate Systems
documentation.

## Usage

``` r
mni152reg_mtx()
```

## Value

4x4 numeric matrix, the MNI305-to-MNI152 registration matrix.

## Note

This is the opposite of using the Wu et al. approach. It is mainly
implemented in this package to allow you to easily check the difference
between the methods.
