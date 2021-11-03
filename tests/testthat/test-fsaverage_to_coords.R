


test_that("MNI152 coords can be mapped to fsaverage space", {

  # The next line is a setup for Tim's test system only and should be removed once this package is official.
  Sys.setenv("FS_HOME"=file.path(Sys.getenv("HOME"), "software/freesurfer/"));

  if(nchar(Sys.getenv("FS_HOME")) == 0L) {
    testthat::skip("No FreeSurfer installation found or FS_HOME environment variable not set correctly.");
  }
  if(! dir.exists(Sys.getenv("FS_HOME"))) {
    testthat::skip("No FreeSurfer installation found at path given in FS_HOME environment variable.");
  }

  query_fsaverage_vertex = c(9092L);
  query_hemis = c("rh");

  mni_coord = fsaverage_vertices_to_mni152_coords(vertices = query_fsaverage_vertex, hemis = query_hemis);


  testthat::expect_equal(length(mni_coord), 3L);
})
