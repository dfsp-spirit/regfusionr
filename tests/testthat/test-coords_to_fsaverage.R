

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
  res_outside_cortex = mni152_coords_to_fsaverage(mni_coord_in_cortex, surface = "white");

  testthat::expect_equal(res_in_cortex$fsaverage_vertices, c(9092));
  testthat::expect_equal(res_outside_cortex$fsaverage_vertices, c(NaN));
})


