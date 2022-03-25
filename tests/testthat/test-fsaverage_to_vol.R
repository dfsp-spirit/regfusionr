
testthat::test_that("Per-vertex data for fsaverage can be projected, using linear interpolation, to an MNI152 volume.", {

  if(! requireNamespace("fsbrain", quietly = TRUE)) {
    testthat::skip("The optional dependency package 'fsbrain' is required for this unit test.");
  }
  if(! requireNamespace("haze", quietly = TRUE)) {
    testthat::skip("The optional dependency package 'haze' is required for this unit test.");
  }

  num_fsaverage_verts_per_hemi = 163842L;

  lh_input = rnorm(num_fsaverage_verts_per_hemi, 3.0, 0.2);
  rh_input = rnorm(num_fsaverage_verts_per_hemi, 3.0, 0.2);
  res = fsaverage_to_vol(lh_input, rh_input);

  testthat::expect_true(is.array(res$lh));
  testthat::expect_true(is.array(res$rh));
  testthat::expect_true(is.array(res$brain));

  testthat::expect_equal(dim(res$lh), c(256L, 256L, 256L));
  testthat::expect_equal(dim(res$rh), c(256L, 256L, 256L));
  testthat::expect_equal(dim(res$brain), c(256L, 256L, 256L));
})
