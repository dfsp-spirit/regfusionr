# regfusionr
R implementation of the registration fusion method for MNI152/Colin27 to fsaverage MNI305 mapping.

## About

This is an R implementation of [Wu et al. (2018)'s registration fusion methods](https://onlinelibrary.wiley.com/doi/full/10.1002/hbm.24213) to project 3D magnetic resonance imaging (MRI) data from standard space volumetric coordinates, either MNI152 or Colin27, to Freesurfer's fsaverage (MNI305). This R implementation is heavily inspired by [Dan Gale's Python implementation](https://github.com/danjgale/reg-fusion) in the [regfusion pypi package](https://pypi.org/project/regfusion/). A huge thank you to Dan Gale and  *Wu et al* for making their excellent tools openly available!

## Documentation

### Quickstart

Use one the following functions:

* `mni152_coords_to_fsaverage()` : Map MNI152 RAS coordinates to fsaverage coordinates. For the coordinates, you also get the closest surface vertex and information on which hemisphere it belongs to.
* `mni152_voxels_to_fsaverage()` : Map MNI152 voxel indices to fsaverage coordinates.
* `vol_to_fsaverage()`: Project the 3D data in an MNI152 or Colin27 volume (in NIFTI or MGH/MGZ format) to fsaverage and obtain per-vertex data (in curv or MGH/MGZ format).


### Usage examples

```
library('regfusionr');

# Get fsaverage coordinates for the MNI152 RAS coordinates 60.0, 0.0, 10.0 and 0.0, 0.0, 0.0:
mni_ras_coords = matrix(c(60, 0, 10, 0, 0, 0), ncol = 3, byrow = TRUE);
res = mni152_coords_to_fsaverage(mni_ras_coords, surface = "white");

# Get fsaverage coordinates for the MNI152 voxel at IJK index 68, 138, 146: 
mni_voxel_ijk = c(68L, 138L, 146L);
res_in_cortex = mni152_voxels_to_fsaverage(mni_voxel_ijk, surface = "white");
```

See the [unit tests](./test/testthat/) for more usage examples, and use the in-built R help (with `?`) to see more details on all the parameters and return values, e.g. `?regfusionr::mni152_coords_to_fsaverage`.


## Limitations

* When projecting volume data to the surface, currently only the 'linear' interpolation method is implemented (which uses trilinear interpolation from the `oce` package). The 'nearest' method, which is required to project labels or atlases, is not available yet. If you know an R function that does it, please let me know.
* Some parts I don't need myself (like Colin27 coordinates to fsaverage vertex mapping) are not yet implemented. Please open an issue if you need them and I will do it.


## Installation

### Minimal installation

Note: If you want to project volume data (e.g., a NIFTI or MGZ volume) to fsaverage surface space, follow the installation instructions below.

If all you want to do is to obtain fsaverage coordinates for MNI152 voxels or coordinates, you can get away without installing the system dependencies because you do not need the trilinear interpolation functions, and all you need to do is:

```
install.packages('remotes');
remotes::install_github('dfsp-spirit/regfusionr');
```

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

From an R session:

```
install.packages('remotes');
remotes::install_github('dfsp-spirit/regfusionr', dependencies = TRUE);
```

You can also use `devtools` instead of `remotes` if you already have it installed. For those who haven't, `remotes` is a lot smaller and faster to install though.

It's unlikely that this package will go to CRAN soon, it requires some data files which are about 100 MB in total size, and CRAN only supports 5 MB. I know one can work around that, but my time for this is limited.


## Unit tests and continuous integration

See the [development information file](./DEVELOP.md).

## Citation

Just cite the original *Wu et al.* paper:

>Wu J, Ngo GH, Greve DN, Li J, He T, Fischl B, Eickhoff SB, Yeo BTT. [**Accurate nonlinear mapping between MNI volumetric and FreeSurfer surface coordinate systems**](http://people.csail.mit.edu/ythomas/publications/2018VolSurfMapping-HBM.pdf), *Human Brain Mapping* 39:3793–3808, 2018.


   
