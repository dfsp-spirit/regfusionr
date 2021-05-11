The files in the coordmap directory are re-encoded versions of Yeo CBIG files:
* `final_warps_FS5.3/allSub_fsaverage_to_FSL_MNI152_FS4.5.0_RF_ANTs_avgMapping.vertex.mat` has been split by hemi and renamed and re-encoded to MGZ format as `FSL_MNI152_FS4.5.0_RF_ANTs_avgMapping.vertex.?h.mgz` (see code below)
* `liberal_cortex_masks_FS5.3/FSL_MNI152_FS4.5.0_cortex_estimate.nii.gz` has been renamed to `FSL_MNI152_FS4.5.0_cortex_estimate.nii.gz` (moved out of subdir only)

These changes were needed because:
1) storing paths with >100 chars in an R pacakge leads to TAR warnings that this is not portable.
2) reading the .mat files requires extra packages in R.

```
dd = R.matlab::readMat("inst/extdata/coordmap/allSub_fsaverage_to_FSL_MNI152_FS4.5.0_RF_ANTs_avgMapping.vertex.mat")
freesurferformats::write.fs.mgh("inst/extdata/coordmap/FSL_MNI152_FS4.5.0_RF_ANTs_avgMapping.vertex.rh.mgz", dd$rh.vertex)
freesurferformats::write.fs.mgh("inst/extdata/coordmap/FSL_MNI152_FS4.5.0_RF_ANTs_avgMapping.vertex.lh.mgz", dd$lh.vertex)
```

