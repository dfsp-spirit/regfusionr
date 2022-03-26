These files have been converted from Matlab files (.mat format) to MGZ format so one does not need another R package to read them. They have also been split by hemisphere.

E.g.,

```R
cbig_warps_dir = "~/develop/CBIG/stable_projects/registration/Wu2017_RegistrationFusion/bin/final_warps_FS5.3";
target_dir = "inst/extdata/coordmap/";
for( rf_type in c('RF_ANTs', 'RF_M3Z')) {
  for( mapping_type in c('FSL_MNI152', 'SPM_Colin27')) {
    # one source file for both hemis.
    source_filename = sprintf("allSub_fsaverage_to_%s_FS4.5.0_%s_avgMapping.prop.mat", mapping_type, rf_type);
    wfile = file.path(cbig_warps_dir, source_filename);
    dd = R.matlab::readMat(wfile); 
    # last command requires some RAM. Make sure to install.packages(c("Matrix", "SparseM")) or 
    # install R.matlab with dependencies=TRUE. Also use a fresh session on low memory systems.
    
    target_file_lh = sprintf("allSub_fsaverage_to_%s_FS4.5.0_%s_avgMapping.lh.mgz", mapping_type, rf_type);
    target_file_rh = sprintf("allSub_fsaverage_to_%s_FS4.5.0_%s_avgMapping.rh.mgz", mapping_type, rf_type);
    
    
    freesurferformats::write.fs.mgh(file.path(target_dir, target_file_lh), dd$lh.coord);
    freesurferformats::write.fs.mgh(file.path(target_dir, target_file_rh), dd$rh.coord);
  }
}
```    
