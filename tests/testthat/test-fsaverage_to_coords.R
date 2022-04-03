


test_that("fsaverage vertices can be mapped to MNI152 space", {

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

  testthat::expect_true(is.matrix(mni_coord));
  testthat::expect_equal(length(mni_coord), 3L);
})


test_that("MNI305 coordinates can be mapped to MNI152 space", {

  # The next line is a setup for Tim's test system only and should be removed once this package is official.
  Sys.setenv("FS_HOME"=file.path(Sys.getenv("HOME"), "software/freesurfer/"));

  if(nchar(Sys.getenv("FS_HOME")) == 0L) {
    testthat::skip("No FreeSurfer installation found or FS_HOME environment variable not set correctly.");
  }
  if(! dir.exists(Sys.getenv("FS_HOME"))) {
    testthat::skip("No FreeSurfer installation found at path given in FS_HOME environment variable.");
  }

  fsaverage_home = file.path(Sys.getenv("FS_HOME"), "subjects");
  query_fsaverage_vertex = c(9092L);
  query_hemis = c("rh");
  rh_orig_surf = freesurferformats::read.fs.surface(file.path(fsaverage_home, "fsaverage", "surf", "rh.orig"));
  query_vertex_coords = rh_orig_surf$vertices[query_fsaverage_vertex, ];

  mni_coord = mni305_coords_to_mni152_coords(query_vertex_coords);

  testthat::expect_true(is.matrix(mni_coord));
  testthat::expect_equal(length(mni_coord), 3L);
})


test_that("MNI305 coordinates can be mapped to Colin27 space", {

  # The next line is a setup for Tim's test system only and should be removed once this package is official.
  Sys.setenv("FS_HOME"=file.path(Sys.getenv("HOME"), "software/freesurfer/"));

  if(nchar(Sys.getenv("FS_HOME")) == 0L) {
    testthat::skip("No FreeSurfer installation found or FS_HOME environment variable not set correctly.");
  }
  if(! dir.exists(Sys.getenv("FS_HOME"))) {
    testthat::skip("No FreeSurfer installation found at path given in FS_HOME environment variable.");
  }

  fsaverage_home = file.path(Sys.getenv("FS_HOME"), "subjects");
  query_fsaverage_vertex = c(9092L);
  query_hemis = c("rh");
  rh_orig_surf = freesurferformats::read.fs.surface(file.path(fsaverage_home, "fsaverage", "surf", "rh.orig"));
  query_vertex_coords = rh_orig_surf$vertices[query_fsaverage_vertex, ];

  colin_coord = mni305_coords_to_colin27_coords(query_vertex_coords);

  testthat::expect_true(is.matrix(colin_coord));
  testthat::expect_equal(length(colin_coord), 3L);
})



