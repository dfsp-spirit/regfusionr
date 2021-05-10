

get_test_file <- function(filename) {
  system.file("extdata", "testdata", filename, package = "regusionr", mustWork = TRUE);
}

test_that("The colin27 data m3z transformation works", {
  expected_res_file_lh = get_test_file('lh.Colin_probMap_ants.central_sulc.allSub_RF_M3Z_Colin27_norm_to_fsaverage.nii.gz');
  expected_res_file_rh = get_test_file('rh.Colin_probMap_ants.central_sulc.allSub_RF_M3Z_Colin27_norm_to_fsaverage.nii.gz');

  input_img = get_test_file('Colin_probMap_ants.central_sulc.nii.gz');
  out_dir = tempdir();
  out = vol_to_fsaverage(input_img, out_dir, template_type='Colin27_norm', rf_type='RF_M3Z')

  testthat::expect_true(file.exists(out$lh));
  testthat::expect_true(file.exists(out$rh));
})


