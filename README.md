# regfusionr
R implementation of registration fusion method for MNI152/Colin to fsaverage mapping.

## About

This is an R implementation of [Wu et al. (2018)'s registration fusion methods](https://onlinelibrary.wiley.com/doi/full/10.1002/hbm.24213) to project 3D magnetic resonance imaging (MRI) data from standard space volumetric coordinates, either MNI152 or Colin27, to Freesurfer's fsaverage (MNI305). This R implementation is heavily inspired by [Dan Gale's Python implementation](https://github.com/danjgale/reg-fusion) in the [regfusion pypi package](https://pypi.org/project/regfusion/). A huge thank you to Dan Gale and  *Wu et al* for making their excellent tools openly available!

## Usage

This is WIP and not usable yet, come back another day.

## Installation

### System Dependencies

You need to install these from your system shell before installing the R package (see below).

For Debian-based Linux distros, run:
```
sudo apt install libudunits2-dev libgdal-dev
```

Required system level packages for other systems:
 
* rpm-based systems (like Fedora, EPEL, ...): `sudo yum install udunits2-devel gdal-devel`
* MacOS (via [brew](https://brew.sh)): `brew install udunits gdal`

### R package

From an R session:

```
install.packages('remotes');
remotes::install_github('dfsp-spirit/regfusionr');
```

It's unlikely that this package will go to CRAN soon, it requires some data files (about 100 MB), and CRAN only supports 5 MB. I'm sure one could work around that somehow, but my time for this is limited.

## Citation

Just cite the original *Wu et al.* paper:

>Wu J, Ngo GH, Greve DN, Li J, He T, Fischl B, Eickhoff SB, Yeo BTT. [**Accurate nonlinear mapping between MNI volumetric and FreeSurfer surface coordinate systems**](http://people.csail.mit.edu/ythomas/publications/2018VolSurfMapping-HBM.pdf), *Human Brain Mapping* 39:3793–3808, 2018.


   
