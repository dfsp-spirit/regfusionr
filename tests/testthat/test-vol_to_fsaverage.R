

get_test_file <- function(filename) {
  system.file("extdata", "testdata", filename, package = "regfusionr", mustWork = TRUE);
}

#' Get float comparison tolerance.
tol <- function() {
  1.0
}

test_that("The Colin27 data projection using m3z registration works", {

  if(! requireNamespace("oce", quietly = TRUE)) {
    testthat::skip("The optional dependency package 'oce' is required for this unit test.");
  }

  # Load expected result data.
  expected_res_file_lh = get_test_file('lh.Colin.RF_M3Z_Colin27_norm_to_fsaverage.nii.gz');
  expected_res_file_rh = get_test_file('rh.Colin.RF_M3Z_Colin27_norm_to_fsaverage.nii.gz');
  expected_lh = freesurferformats::read.fs.volume(expected_res_file_lh);
  expected_rh = freesurferformats::read.fs.volume(expected_res_file_rh);
  testthat::expect_true(is.array(expected_lh));
  testthat::expect_true(is.array(expected_rh));
  testthat::expect_equal(length(as.vector(expected_lh)), 163842);
  testthat::expect_equal(length(as.vector(expected_rh)), 163842);

  # Load input image for test.
  input_img = get_test_file('Colin_probMap_ants.central_sulc.nii.gz');
  out_dir = tempdir();
  outfiles = vol_to_fsaverage(input_img, template_type='Colin27_norm', rf_type='RF_M3Z', out_type = 'curv', out_dir = out_dir);

  # Compare results with expected result data.
  testthat::expect_true(file.exists(outfiles$lh));
  testthat::expect_true(file.exists(outfiles$rh));
  actual_lh = freesurferformats::read.fs.morph(outfiles$lh);
  actual_rh = freesurferformats::read.fs.morph(outfiles$rh);
  testthat::expect_true(is.vector(actual_lh));
  testthat::expect_true(is.vector(actual_rh));
  testthat::expect_equal(length(actual_lh), 163842);
  testthat::expect_equal(length(actual_rh), 163842);

  testthat::expect_equal(actual_lh, as.vector(expected_lh), tolerance = tol());
  testthat::expect_equal(actual_rh, as.vector(expected_rh), tolerance = tol());

  #fsbrain::vis.data.on.fsaverage(morph_data_lh = actual_lh, morph_data_rh = actual_rh);
  #
  # fsbrain::vis.symmetric.data.on.subject("~/software/freesurfer/subjects", "fsaverage", morph_data_lh = actual_lh, morph_data_rh = actual_rh, surface = "inflated");
})


test_that("The MNI152 data projection using ANTs registration works using output curv files", {

  if(! requireNamespace("oce", quietly = TRUE)) {
    testthat::skip("The optional dependency package 'oce' is required for this unit test.");
  }

  # Load expected result data.
  expected_res_file_lh = get_test_file('lh.MNI.RF_ANTs_MNI152_orig_to_fsaverage.nii.gz');
  expected_res_file_rh = get_test_file('rh.MNI.RF_ANTs_MNI152_orig_to_fsaverage.nii.gz');
  expected_lh = freesurferformats::read.fs.volume(expected_res_file_lh);
  expected_rh = freesurferformats::read.fs.volume(expected_res_file_rh);
  testthat::expect_true(is.array(expected_lh));
  testthat::expect_true(is.array(expected_rh));
  testthat::expect_equal(length(as.vector(expected_lh)), 163842);
  testthat::expect_equal(length(as.vector(expected_rh)), 163842);

  # Load input image for test.
  input_img = get_test_file('MNI_probMap_ants.central_sulc.nii.gz');
  out_dir = tempdir();
  outfiles = vol_to_fsaverage(input_img, template_type='MNI152_orig', rf_type='RF_ANTs', out_type = 'curv', out_dir = out_dir);

  # Compare results with expected result data.
  testthat::expect_true(file.exists(outfiles$lh));
  testthat::expect_true(file.exists(outfiles$rh));
  actual_lh = freesurferformats::read.fs.morph(outfiles$lh);
  actual_rh = freesurferformats::read.fs.morph(outfiles$rh);
  testthat::expect_true(is.vector(actual_lh));
  testthat::expect_true(is.vector(actual_rh));
  testthat::expect_equal(length(actual_lh), 163842);
  testthat::expect_equal(length(actual_rh), 163842);

  testthat::expect_equal(actual_lh, as.vector(expected_lh), tolerance = tol());
  testthat::expect_equal(actual_rh, as.vector(expected_rh), tolerance = tol());

  #fsbrain::vis.data.on.fsaverage(morph_data_lh = actual_lh, morph_data_rh = actual_rh);
  #
  # fsbrain::vis.symmetric.data.on.subject("~/software/freesurfer/subjects", "fsaverage", morph_data_lh = actual_lh, morph_data_rh = actual_rh, surface = "inflated");
})


test_that("The MNI152 data projection using ANTs registration works returning data", {

  if(! requireNamespace("oce", quietly = TRUE)) {
    testthat::skip("The optional dependency package 'oce' is required for this unit test.");
  }

  # Load expected result data.
  expected_res_file_lh = get_test_file('lh.MNI.RF_ANTs_MNI152_orig_to_fsaverage.nii.gz');
  expected_res_file_rh = get_test_file('rh.MNI.RF_ANTs_MNI152_orig_to_fsaverage.nii.gz');
  expected_lh = freesurferformats::read.fs.volume(expected_res_file_lh);
  expected_rh = freesurferformats::read.fs.volume(expected_res_file_rh);
  testthat::expect_true(is.array(expected_lh));
  testthat::expect_true(is.array(expected_rh));
  testthat::expect_equal(length(as.vector(expected_lh)), 163842);
  testthat::expect_equal(length(as.vector(expected_rh)), 163842);

  # Load input image for test.
  input_img = get_test_file('MNI_probMap_ants.central_sulc.nii.gz');
  out_dir = NULL;
  outdata = vol_to_fsaverage(input_img, template_type='MNI152_orig', rf_type='RF_ANTs', out_type = 'curv', out_dir = out_dir);

  actual_lh = outdata$lh;
  actual_rh = outdata$rh;
  testthat::expect_true(is.vector(actual_lh));
  testthat::expect_true(is.vector(actual_rh));
  testthat::expect_equal(length(actual_lh), 163842);
  testthat::expect_equal(length(actual_rh), 163842);

  testthat::expect_equal(actual_lh, as.vector(expected_lh), tolerance = tol());
  testthat::expect_equal(actual_rh, as.vector(expected_rh), tolerance = tol());

  #fsbrain::vis.data.on.fsaverage(morph_data_lh = actual_lh, morph_data_rh = actual_rh);
  #
  # fsbrain::vis.symmetric.data.on.subject("~/software/freesurfer/subjects", "fsaverage", morph_data_lh = actual_lh, morph_data_rh = actual_rh, surface = "inflated");
  #
  # plot differences expected - actual:
  #fsbrain::vis.symmetric.data.on.subject("~/software/freesurfer/subjects", "fsaverage", morph_data_lh = (expected_lh - actual_lh), morph_data_rh = (expected_rh - actual_rh), surface = "inflated");
})

# TODO: Add tests for 4D input data, should write several output curv files or a single output MGZ file.
