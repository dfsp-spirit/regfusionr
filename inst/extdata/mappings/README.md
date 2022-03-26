These files have been converted from Matlab files (.mat format) to MGZ format so one does not need another R package to read them. They have also been split by hemisphere.

E.g.,

```R
cbig_warps_dir = "~/develop/CBIG/stable_projects/registration/Wu2017_RegistrationFusion/bin/final_warps_FS5.3";
wfile = file.path(cbig_warps_dir, "allSub_fsaverage_to_FSL_MNI152_FS4.5.0_RF_ANTs_avgMapping.prop.mat");
dd = R.matlab::readMat(wfile); # requires some RAM. Make sure to install.packages(c("Matrix", "SparseM")) or install R.matlab with dependencies=TRUE. Also use a fresh session on low memory systems.

freesurferformats::write.fs.mgh("inst/extdata/coordmap/allSub_fsaverage_to_FSL_MNI152_FS4.5.0_RF_ANTs_avgMapping.lh.mgz", dd$lh.coord)
freesurferformats::write.fs.mgh("inst/extdata/coordmap/allSub_fsaverage_to_FSL_MNI152_FS4.5.0_RF_ANTs_avgMapping.rh.mgz", dd$rh.coord)
```    
