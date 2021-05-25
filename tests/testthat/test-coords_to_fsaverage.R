

test_that("MNI152 coords can be mapped to fsaverage space", {

  # The next line is a setup for Tim's test system only and should be removed once this package is official.
  Sys.setenv("FS_HOME"=file.path(Sys.getenv("HOME"), "software/freesurfer/"));

  if(nchar(Sys.getenv("FS_HOME")) == 0L) {
    testthat::skip("No FreeSurfer installation found or FS_HOME environment variable not set correctly.");
  }
  if(! dir.exists(Sys.getenv("FS_HOME"))) {
    testthat::skip("No FreeSurfer installation found at path given in FS_HOME environment variable.");
  }

  mni_coord_in_cortex = c(60, 0, 10);
  mni_coord_outside_cortex = c(0, 0, 0);

  res_in_cortex = mni152_coords_to_fsaverage(mni_coord_in_cortex, surface = "white");
  res_outside_cortex = mni152_coords_to_fsaverage(mni_coord_outside_cortex, surface = "white");

  testthat::expect_equal(res_in_cortex$fsaverage_vertices, c(9092));
  testthat::expect_equal(res_in_cortex$query_mni_voxels, matrix(c(68L, 138L, 146L), ncol=3, byrow = TRUE));

  testthat::expect_equal(res_outside_cortex$fsaverage_vertices, c(NaN));
})


test_that("MNI152 coord mapping works with matrix iunput of several coords at once", {

  # The next line is a setup for Tim's test system only and should be removed once this package is official.
  Sys.setenv("FS_HOME"=file.path(Sys.getenv("HOME"), "software/freesurfer/"));

  if(nchar(Sys.getenv("FS_HOME")) == 0L) {
    testthat::skip("No FreeSurfer installation found or FS_HOME environment variable not set correctly.");
  }
  if(! dir.exists(Sys.getenv("FS_HOME"))) {
    testthat::skip("No FreeSurfer installation found at path given in FS_HOME environment variable.");
  }

  mni_coords = matrix(c(60, 0, 10, 0, 0, 0), ncol = 3, byrow = TRUE);

  res = mni152_coords_to_fsaverage(mni_coords, surface = "white");
  testthat::expect_equal(res$fsaverage_vertices, c(9092, NaN));
  testthat::expect_equal(length(res$fsaverage_vertices), 2L);
  testthat::expect_equal(length(res$hemi), 2L);
  testthat::expect_equal(nrow(res$query_mni_coords), 2L);
  testthat::expect_equal(nrow(res$query_mni_voxels), 2L);
  testthat::expect_equal(nrow(res$fsaverage_coords), 2L);
})


# Keep in mind that the voxel indices are specific for the template file, and thus of very
# limited use in general.
test_that("MNI152 voxels based on the demo input file can be mapped to fsaverage space", {

  # The next line is a setup for Tim's test system only and should be removed once this package is official.
  Sys.setenv("FS_HOME"=file.path(Sys.getenv("HOME"), "software/freesurfer/"));

  if(nchar(Sys.getenv("FS_HOME")) == 0L) {
    testthat::skip("No FreeSurfer installation found or FS_HOME environment variable not set correctly.");
  }
  if(! dir.exists(Sys.getenv("FS_HOME"))) {
    testthat::skip("No FreeSurfer installation found at path given in FS_HOME environment variable.");
  }

  mni_voxel_ijk = c(68L, 138L, 146L);

  res_in_cortex = regfusionr:::mni152_voxels_to_fsaverage(mni_voxel_ijk, surface = "white");

  testthat::expect_equal(res_in_cortex$fsaverage_vertices, c(9092));
  testthat::expect_equal(res_in_cortex$query_mni_voxels, matrix(mni_voxel_ijk, ncol=3, byrow = TRUE));
})


