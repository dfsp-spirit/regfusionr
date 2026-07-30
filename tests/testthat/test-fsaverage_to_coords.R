


test_that("fsaverage vertices can be mapped to MNI152 space", {

  query_fsaverage_vertex = c(9092L);
  query_hemis = c("rh");

  mni_coord = fsaverage_vertices_to_mni152_coords(vertices = query_fsaverage_vertex, hemis = query_hemis);

  testthat::expect_true(is.matrix(mni_coord));
  testthat::expect_equal(length(mni_coord), 3L);
})


test_that("fsaverage vertices can be mapped to Colin27 space", {

  query_fsaverage_vertex = c(9092L);
  query_hemis = c("rh");

  c27_coord = fsaverage_vertices_to_colin27_coords(vertices = query_fsaverage_vertex, hemis = query_hemis);

  testthat::expect_true(is.matrix(c27_coord));
  testthat::expect_equal(length(c27_coord), 3L);
})


test_that("MNI305 coordinates can be mapped to MNI152 space", {

  fs_home = Sys.getenv("FREESURFER_HOME");
  if(nchar(fs_home) == 0L || ! dir.exists(fs_home)) {
    testthat::skip("FreeSurfer not available: FREESURFER_HOME not set or invalid.");
  }

  fsaverage_surf_file = file.path(fs_home, "subjects", "fsaverage", "surf", "rh.orig");
  if(! file.exists(fsaverage_surf_file)) {
    testthat::skip(sprintf("fsaverage surface file not found at '%s'.", fsaverage_surf_file));
  }

  query_fsaverage_vertex = c(9092L);
  query_hemis = c("rh");
  rh_orig_surf = freesurferformats::read.fs.surface(fsaverage_surf_file);
  query_vertex_coords = rh_orig_surf$vertices[query_fsaverage_vertex, ];

  mni_coord = mni305_coords_to_mni152_coords(query_vertex_coords, fs_home = fs_home);

  testthat::expect_true(is.matrix(mni_coord));
  testthat::expect_equal(length(mni_coord), 3L);
})


test_that("MNI305 coordinates can be mapped to Colin27 space", {

  fs_home = Sys.getenv("FREESURFER_HOME");
  if(nchar(fs_home) == 0L || ! dir.exists(fs_home)) {
    testthat::skip("FreeSurfer not available: FREESURFER_HOME not set or invalid.");
  }

  fsaverage_surf_file = file.path(fs_home, "subjects", "fsaverage", "surf", "rh.orig");
  if(! file.exists(fsaverage_surf_file)) {
    testthat::skip(sprintf("fsaverage surface file not found at '%s'.", fsaverage_surf_file));
  }

  query_fsaverage_vertex = c(9092L);
  query_hemis = c("rh");
  rh_orig_surf = freesurferformats::read.fs.surface(fsaverage_surf_file);
  query_vertex_coords = rh_orig_surf$vertices[query_fsaverage_vertex, ];

  colin_coord = mni305_coords_to_colin27_coords(query_vertex_coords, fs_home = fs_home);

  testthat::expect_true(is.matrix(colin_coord));
  testthat::expect_equal(length(colin_coord), 3L);
})



