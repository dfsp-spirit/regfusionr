

get_test_file <- function(filename) {
  system.file("extdata", "testdata", filename, package = "regfusionr", mustWork = TRUE);
}

test_that("The colin27 data m3z transformation works", {
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
  outfiles = vol_to_fsaverage(input_img, out_dir, template_type='Colin27_norm', rf_type='RF_M3Z', out_type = 'curv');

  # Compare results with expected result data.
  testthat::expect_true(file.exists(outfiles$lh));
  testthat::expect_true(file.exists(outfiles$rh));
  actual_lh = freesurferformats::read.fs.morph(outfiles$lh);
  actual_rh = freesurferformats::read.fs.morph(outfiles$rh);
  testthat::expect_true(is.vector(actual_lh));
  testthat::expect_true(is.vector(actual_rh));
  testthat::expect_equal(length(actual_lh), 163842);
  testthat::expect_equal(length(actual_rh), 163842);

  testthat::expect_equal(actual_lh, as.vector(expected_lh), tolerance = 1e-5);
  testthat::expect_equal(actual_rh, as.vector(expected_rh), tolerance = 1e-5);

  #fsbrain::vis.data.on.fsaverage(morph_data_lh = actual_lh, morph_data_rh = actual_rh);
  #
  # fsbrain::vis.symmetric.data.on.subject("~/software/freesurfer/subjects", "fsaverage", morph_data_lh = actual_lh, morph_data_rh = actual_rh, surface = "inflated");
})


