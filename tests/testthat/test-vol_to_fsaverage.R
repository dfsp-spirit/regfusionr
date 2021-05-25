

get_test_file <- function(filename) {
  system.file("extdata", "testdata", filename, package = "regfusionr", mustWork = TRUE);
}

#' Get float comparison tolerance.
tol <- function() {
  1.0
}

# Create a 4D volume from a 3D volume by copying the 3D volume several times (treating the copies as new frames).
vol3dto4d <- function(input_volume_3d, num_frames = 3L) {
  input_volume_4d = input_volume_3d;
  input_volume_4d$header$voldim[4] = num_frames;
  input_volume_4d$header$voldim_orig[4] = num_frames;
  data_4d = array(rep(NA, prod(dim(input_volume_3d$data))*num_frames), dim=c(dim(input_volume_3d$data),num_frames));
  for(frame_idx in seq_len(num_frames)) {
    data_4d[,,,frame_idx] = input_volume_3d$data;
  }
  input_volume_4d$data = data_4d;
  return(input_volume_4d);
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
})



test_that("The MNI152 data projection using ANTs registration works using output gii files", {

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
  outfiles = vol_to_fsaverage(input_img, template_type='MNI152_orig', rf_type='RF_ANTs', out_type = 'gii', out_dir = out_dir);

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
})


test_that("The MNI152 data projection using ANTs registration works returning data (from 3D input)", {

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


  # Show the source probability map as a volume in 3D:
  # vol = freesurferformats::read.fs.volume(input_img, drop = TRUE);
  # fsbrain::volvis.contour(vol, level = 0.2, color = "gray");
  #
  #fsbrain::vis.data.on.fsaverage(morph_data_lh = actual_lh, morph_data_rh = actual_rh);
  #
  # fsbrain::vis.symmetric.data.on.subject("~/software/freesurfer/subjects", "fsaverage", morph_data_lh = actual_lh, morph_data_rh = actual_rh, surface = "inflated");
  #
  # plot differences expected - actual:
  #fsbrain::vis.symmetric.data.on.subject("~/software/freesurfer/subjects", "fsaverage", morph_data_lh = (expected_lh - actual_lh), morph_data_rh = (expected_rh - actual_rh), surface = "inflated");
})


test_that("The MNI152 data projection using ANTs registration works returning data (from 4D input)", {

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

  # Load input image for test and construct 4D image from it.
  input_img = get_test_file('MNI_probMap_ants.central_sulc.nii.gz');
  input_volume_3d = freesurferformats::read.fs.volume(input_img, with_header = TRUE, drop_empty_dims = TRUE);
  testthat::expect_equal(length(dim(input_volume_3d$data)), 3L);
  # Construct 4D image by using several copies of the 3D image as frames.
  num_frames = 5L;
  input_volume_4d = vol3dto4d(input_volume_3d, num_frames = num_frames);
  testthat::expect_equal(length(dim(input_volume_4d$data)), 4L);


  out_dir = NULL;
  outdata = vol_to_fsaverage(input_volume_4d, template_type='MNI152_orig', rf_type='RF_ANTs', out_type = 'mgz', out_dir = out_dir);

  actual_lh = outdata$lh;
  actual_rh = outdata$rh;
  testthat::expect_true(is.array(actual_lh));
  testthat::expect_true(is.array(actual_rh));
  testthat::expect_equal(dim(actual_lh), c(163842,1,1,num_frames));
  testthat::expect_equal(dim(actual_rh), c(163842,1,1,num_frames));
})





test_that("Projecting input 3D and 4D data leads to the expected return types.", {

  if(! requireNamespace("oce", quietly = TRUE)) {
    testthat::skip("The optional dependency package 'oce' is required for this unit test.");
  }


  # Load input image for test.
  input_img = get_test_file('MNI_probMap_ants.central_sulc.nii.gz');
  input_volume_3d = freesurferformats::read.fs.volume(input_img, with_header = TRUE, drop_empty_dims = TRUE);
  testthat::expect_equal(length(dim(input_volume_3d$data)), 3L);

  # Construct 4D image by using several copies of the 3D image as frames.
  num_frames = 5L;
  input_volume_4d = vol3dto4d(input_volume_3d, num_frames = num_frames);
  testthat::expect_equal(length(dim(input_volume_4d$data)), 4L);

  out_dir = NULL;

  mapping = ".avgMapping_allSub_RF_ANTs_MNI152_orig_to_fsaverage.txt";
  mapping_file = system.file("extdata", "mappings", sprintf("%s%s", "lh", mapping), package = "regfusionr", mustWork = TRUE);
  ras = t(as.matrix(data.table::fread(mapping_file, nrows = 3, header = FALSE)));
  affine = freesurferformats::mghheader.ras2vox(input_volume_3d$header);

  projected_from_3d = regfusionr:::project_data(input_volume_3d$data, affine, ras);
  projected_from_4d = regfusionr:::project_data(input_volume_4d$data, affine, ras);

  testthat::expect_true(is.vector(projected_from_3d));
  testthat::expect_true(is.array(projected_from_4d));
  testthat::expect_equal(dim(projected_from_4d), c(163842L, 1, 1, num_frames));
})
