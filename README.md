# regfusionr
R implementation of the registration fusion method for MNI152 and Colin27 to fsaverage/MNI305 mapping.

This package supports easy mapping of neuroimgaging data between the volume and surface templates used by the most common software package for structural neuroimaging in R:

* [FreeSurfer](freesurfer.net/) surface space: the fsaverage template, in MNI305 space
* [FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki) volume space: the MNI152 template
* [SPM](https://www.fil.ion.ucl.ac.uk/spm/software/) volumne space: the Colin27 template

## About

This is an R implementation of [Wu et al. (2018)'s registration fusion methods](https://onlinelibrary.wiley.com/doi/full/10.1002/hbm.24213) to project 3D magnetic resonance imaging (MRI) data from standard space volumetric coordinates, either MNI152 or Colin27, to Freesurfer's fsaverage (MNI305), and the other way around. Using this non-linear approach gives higher accuracy than the linear transformation with a 4x4 matrix. See [the paper](https://onlinelibrary.wiley.com/doi/full/10.1002/hbm.24213) for details.

This R implementation is heavily inspired by [Dan Gale's Python implementation](https://github.com/danjgale/reg-fusion) in the [regfusion pypi package](https://pypi.org/project/regfusion/). A huge thank you to Dan Gale and  Wu *et al.* for making their excellent tools openly available!

## Documentation

### Quickstart

The API of the `regfusionr` package consists of the following functions:

* `mni152_coords_to_fsaverage()` : Map MNI152 RAS coordinates to fsaverage coordinates. For the coordinates, you also get the closest surface vertex and information on which hemisphere it belongs to.
* `vol_to_fsaverage()`: Project the 3D data in an MNI152 or Colin27 volume (in NIFTI or MGH/MGZ format) to fsaverage and obtain per-vertex data (in curv or MGH/MGZ format).
* `fsaverage_vertices_to_mni152_coords()`: Map fsaverage vertex indices to MNI152 coordinates.
* `fsaverage_vertices_to_colin27_coords()`: Map fsaverage vertex indices to Colin27 coordinates.
* `fsaverage_to_vol()`: Project or map per-vertex values from the fsaverage surface to the cortex voxels of an MNI volume. Also supports fsaverage6 and fsaverage5 as sources. This is experimental and work in progress, use with care.

### Usage examples

```
library('regfusionr');

# Get fsaverage coordinates for the MNI152 RAS coordinates 60.0, 0.0, 10.0 and 0.0, 0.0, 0.0:
mni_ras_coords = matrix(c(60, 0, 10, 0, 0, 0), ncol = 3, byrow = TRUE);
res = mni152_coords_to_fsaverage(mni_ras_coords, surface = "white");
```

See the [unit tests](./test/testthat/) for more usage examples, and use the in-built R help (with `?`) to see more details on all the parameters and return values, e.g. `?regfusionr::mni152_coords_to_fsaverage`.


## Limitations

* When projecting volume data to the surface, currently only the 'linear' interpolation method is implemented (which uses trilinear interpolation from the `oce` package). This method is suitable for continuous data. The 'nearest' method, which is required to project labels or atlases (categorical data represented by integers), is not available yet. If you know an R function that does it, please let me know.


## Installation

### Minimal installation

This is a minimal installation with reduced functionality that does not require you to install the optional dependencies, which can be  hard to install.

Note: If you want to project volume data (e.g., a NIFTI or MGZ volume) to fsaverage surface space, follow the Full Installation instructions below.

If all you want to do is to obtain fsaverage coordinates for MNI152 voxels or coordinates, you can get away without installing the system dependencies because you do not need the trilinear interpolation functions, and all you need to do is:

via `remotes`:

```
install.packages('remotes');
remotes::install_github('dfsp-spirit/regfusionr');
```

or using [R universe](https://r-universe.dev/):

```r
options(repos = c(
    dfspspirit = 'https://dfsp-spirit.r-universe.dev',
    CRAN = 'https://cloud.r-project.org'))

install.packages('regfusionr');
```

I prefer R universe.

### Full installation

This is required to use the `vol_to_fsaverage()` function.

#### System Dependencies

You need to install these from your system shell before installing the R package (see below).

For Debian-based Linux distros, run:
```
sudo apt install libudunits2-dev libgdal-dev
```

Required system level packages for other systems:

* rpm-based systems (like Fedora, EPEL, ...): `sudo yum install udunits2-devel gdal-devel`
* MacOS (via [brew](https://brew.sh)): `brew install udunits gdal`

#### R package

From an R session, via `remotes`:

```
install.packages('remotes');
remotes::install_github('dfsp-spirit/regfusionr', dependencies = TRUE);
```

or using [R universe](https://r-universe.dev/):

```r
options(repos = c(
    dfspspirit = 'https://dfsp-spirit.r-universe.dev',
    CRAN = 'https://cloud.r-project.org'))

install.packages('regfusionr', dependencies = TRUE);
```

I prefer R universe.

You can also use `devtools` instead of `remotes` if you already have it installed. For those who haven't, `remotes` is a lot smaller and faster to install though.

It's unlikely that this package will go to CRAN soon, it requires some data files which are about 100 MB in total size, and CRAN only supports 5 MB. I know one can work around that, but my time for this is limited. If you prefer to install without remotes/devtools, you can also get regfusionr from [my R universe repo](https://dfsp-spirit.r-universe.dev).


## Unit tests and continuous integration

See the [development information file](./DEVELOP.md).

![Vis1](./web/output_vol_to_fsaverage.png?raw=true "Projection of a central sulcus probability map from a volume in MNI152 space to the fsaverage surface.")
**Fig. 1** *Visualization of a central sulcus probability map on the fsaverage surface. The data has been obtained by projecting [this probability map in a volume in MNI152 space](./inst/extdata/testdata/MNI_probMap_ants.central_sulc.nii.gz) to the fsaverage surface with the `vol_to_fsaverage` function.*


## Alternatives

If all you need is the coordinate mapping *and* you are fine with a less acurate result, you can use the matrix listed in the FreeSurfer documentation on Coordinate Systems. See [section 8b on this website](https://surfer.nmr.mgh.harvard.edu/fswiki/CoordinateSystems). The method is implemented in this package in the function `linear_fsaverage_coords_to_MNI152_coords()`. The difference between the results is shown below.


![Vis2](./web/regfusionr_vs_linear.png?raw=true "Difference between the regfusionr approach and the linear method. The coordinates of all fsaverage vertices were mapped to MNI152 space using both the regfusionr method and the linear method. The difference between the two methods was computed as the Euclidean distance between the resulting MNI152 coordinates.")
**Fig. 2** *Difference between the regfusionr approach and the linear method. The coordinates of all fsaverage vertices were mapped to MNI152 space using both the regfusionr method and the linear method. The difference between the two methods was computed as the Euclidean distance between the resulting MNI152 coordinates.*

The code used to produce the comparison figure is available [in this unit test](tests/testthat/test-compare_linear_to_regfusion.R).

## Citation

Just cite the original *Wu et al.* paper:

>Wu J, Ngo GH, Greve DN, Li J, He T, Fischl B, Eickhoff SB, Yeo BTT. [**Accurate nonlinear mapping between MNI volumetric and FreeSurfer surface coordinate systems**](http://people.csail.mit.edu/ythomas/publications/2018VolSurfMapping-HBM.pdf), *Human Brain Mapping* 39:3793???3808, 2018.
