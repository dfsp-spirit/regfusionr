These files have been converted from Matlab files (.mat format) to CSV tables so one does not need another R package to read them.

E.g.,

```R
cbig_warps_dir = "~/develop/CBIG/stable_projects/registration/Wu2017_RegistrationFusion/bin/final_warps_FS5.3";
wfile = file.path(cbig_warps_dir, "allSub_fsaverage_to_FSL_MNI152_FS4.5.0_RF_ANTs_avgMapping.prop.mat");
dd = R.matlab::readMat(wfile);
```    
