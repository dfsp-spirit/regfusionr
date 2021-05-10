# regfusionr
R implementation of registration fusion method for MNI152 - 305 mapping.

## About

This is an R implementation of [Wu et al. (2018)'s registration fusion methods](https://onlinelibrary.wiley.com/doi/full/10.1002/hbm.24213) to project MRI data from standard volumetric coordinates, either MNI152 or Colin27, to Freesurfer's fsaverage. This R implementation is also heavily inspired by Dan Gale's Python implementation in the [regfusion](https://github.com/danjgale/reg-fusion) package. A huge thank you to Dan Gale and  *Wu et al* for making their excellent tools openly available!

## Usage

This is WIP and not usable yet, come back another day.

## Installation

### Sys Deps

For Debian-based Linux distros, run:
```
sudo apt install libudunits2-dev libgdal-dev
```

Required system level packages for other systems:
 
* rpm-based systems (like Fedora, EPEL, ...): `sudo yum install udunits2-devel gdal-devel`
* MacOS (via [brew](https://brew.sh)): `brew install udunits gdal`

## Citation

Just cite the original *Wu et al.* paper:

>Wu J, Ngo GH, Greve DN, Li J, He T, Fischl B, Eickhoff SB, Yeo BTT. [**Accurate nonlinear mapping between MNI volumetric and FreeSurfer surface coordinate systems**](http://people.csail.mit.edu/ythomas/publications/2018VolSurfMapping-HBM.pdf), *Human Brain Mapping* 39:3793–3808, 2018.


   
