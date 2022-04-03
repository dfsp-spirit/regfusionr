
testthat::test_that("Per-vertex data for fsaverage can be projected, using linear interpolation, to an MNI152 volume.", {

  testthat::skip("too slow");

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
  res = fsaverage_to_vol(lh_input, rh_input, target_space = "FSL_MNI152", out_dir = NULL, interp = "linear");

  # check results
  testthat::expect_true(freesurferformats::is.fs.volume(res$projected));
  testthat::expect_true(is.array(res$projected$data));
  testthat::expect_true(freesurferformats::is.fs.volume(res$projected_seg));
  testthat::expect_true(is.array(res$projected_seg$data));

  mni152_cortex_vol_dim = c(256L, 256L, 256L);

  testthat::expect_equal(dim(res$projected$data), mni152_cortex_vol_dim);
  testthat::expect_equal(dim(res$projected_seg$data), mni152_cortex_vol_dim);
})


testthat::test_that("Per-vertex label data for fsaverage can be projected, using nearest neighbor interpolation, to an MNI152 volume.", {

  if(! requireNamespace("fsbrain", quietly = TRUE)) {
    testthat::skip("The optional dependency package 'fsbrain' is required for this unit test.");
  }
  if(! requireNamespace("haze", quietly = TRUE)) {
    testthat::skip("The optional dependency package 'haze' is required for this unit test.");
  }

  # Get some per-vertex data for fsaverage. We use the fsbrain download option.
  fsbrain::download_fsaverage(TRUE);
  lh_input = sample.int(20, size = 163842L, replace = TRUE); # could load a label instead.
  rh_input = sample.int(20, size = 163842L, replace = TRUE); # could load a label instead.

  # run function
  res = fsaverage_to_vol(lh_input, rh_input, target_space = "FSL_MNI152", out_dir = NULL, interp = "nearest");

  # check results
  testthat::expect_true(freesurferformats::is.fs.volume(res$projected));
  testthat::expect_true(is.array(res$projected$data));
  testthat::expect_true(freesurferformats::is.fs.volume(res$projected_seg));
  testthat::expect_true(is.array(res$projected_seg$data));

  mni152_cortex_vol_dim = c(256L, 256L, 256L);

  testthat::expect_equal(dim(res$projected$data), mni152_cortex_vol_dim);
  testthat::expect_equal(dim(res$projected_seg$data), mni152_cortex_vol_dim);
})
