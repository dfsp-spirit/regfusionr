
testthat::test_that("Per-vertex data for fsaverage can be projected, using linear interpolation, to an MNI152 volume.", {

  if(! requireNamespace("fsbrain", quietly = TRUE)) {
    testthat::skip("The optional dependency package 'fsbrain' is required for this unit test.");
  }
  if(! requireNamespace("haze", quietly = TRUE)) {
    testthat::skip("The optional dependency package 'haze' is required for this unit test.");
  }

  # Get some per-vertex data for fsaverage. We use the fsbrain download option.
  fsbrain::download_fsaverage(TRUE);
  lh_pvd_file = fsbrain::get_optional_data_filepath("subjects_dir/fsaverage/surf/lh.curv");
  rh_pvd_file = fsbrain::get_optional_data_filepath("subjects_dir/fsaverage/surf/rh.curv");
  lh_input = freesurferformats::read.fs.curv(lh_pvd_file);
  rh_input = freesurferformats::read.fs.curv(rh_pvd_file);

  # run function
  res = fsaverage_to_vol(lh_input, rh_input, out_dir = NULL);

  # check results
  testthat::expect_true(is.array(res$lh));
  testthat::expect_true(is.array(res$rh));
  testthat::expect_true(is.array(res$both));

  mni152_cortex_vol_dim = c(256L, 256L, 256L);

  testthat::expect_equal(dim(res$lh), mni152_cortex_vol_dim);
  testthat::expect_equal(dim(res$rh), mni152_cortex_vol_dim);
  testthat::expect_equal(dim(res$both), mni152_cortex_vol_dim);
})
