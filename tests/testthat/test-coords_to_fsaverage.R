

test_that("MNI152 coords can be mapped to fsaverage space", {
  fs_home = Sys.getenv("FREESURFER_HOME");
  if(nchar(fs_home) == 0L || ! dir.exists(fs_home)) {
    testthat::skip("FreeSurfer not available: FREESURFER_HOME not set or invalid.");
  }

  mni_coord_in_cortex = c(60.0, 0.0, 10.0);
  mni_coord_outside_cortex = c(0.0, 0.0, 0.0);

  res_in_cortex = mni152_coords_to_fsaverage(mni_coord_in_cortex, surface = "white");
  res_outside_cortex = mni152_coords_to_fsaverage(mni_coord_outside_cortex, surface = "white");

  testthat::expect_equal(res_in_cortex$fsaverage_vertices, c(9092));
  testthat::expect_equal(res_in_cortex$query_mni_voxels, matrix(c(68L, 138L, 146L), ncol=3, byrow = TRUE));

  testthat::expect_equal(res_outside_cortex$fsaverage_vertices, c(NaN));
})


# test_that("Colin27 coords can be mapped to fsaverage space", {
#
#   # The next line is a setup for Tim's test system only and should be removed once this package is official.
#   Sys.setenv("FS_HOME"=file.path(Sys.getenv("HOME"), "software/freesurfer/"));
#
#   if(nchar(Sys.getenv("FS_HOME")) == 0L) {
#     testthat::skip("No FreeSurfer installation found or FS_HOME environment variable not set correctly.");
#   }
#   if(! dir.exists(Sys.getenv("FS_HOME"))) {
#     testthat::skip("No FreeSurfer installation found at path given in FS_HOME environment variable.");
#   }
#
#   mni_coord_in_cortex = c(60.0, 0.0, 10.0);
#   mni_coord_outside_cortex = c(0.0, 0.0, 0.0);
#
#   res_in_cortex = colin27_coords_to_fsaverage(mni_coord_in_cortex, surface = "white");
#   res_outside_cortex = colin27_coords_to_fsaverage(mni_coord_outside_cortex, surface = "white");
#
#   testthat::expect_equal(res_in_cortex$fsaverage_vertices, c(9092));
#   testthat::expect_equal(res_in_cortex$query_mni_voxels, matrix(c(68L, 138L, 146L), ncol=3, byrow = TRUE));
#
#   testthat::expect_equal(res_outside_cortex$fsaverage_vertices, c(NaN));
# })


test_that("MNI152 coord mapping works with matrix input of several coords at once", {
  fs_home = Sys.getenv("FREESURFER_HOME");
  if(nchar(fs_home) == 0L || ! dir.exists(fs_home)) {
    testthat::skip("FreeSurfer not available: FREESURFER_HOME not set or invalid.");
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
  fs_home = Sys.getenv("FREESURFER_HOME");
  if(nchar(fs_home) == 0L || ! dir.exists(fs_home)) {
    testthat::skip("FreeSurfer not available: FREESURFER_HOME not set or invalid.");
  }

  mni_voxel_ijk = c(68L, 138L, 146L);

  res_in_cortex = regfusionr:::mni152_voxels_to_fsaverage(mni_voxel_ijk, surface = "white");

  testthat::expect_equal(res_in_cortex$fsaverage_vertices, c(9092));
  testthat::expect_equal(res_in_cortex$query_mni_voxels, matrix(mni_voxel_ijk, ncol=3, byrow = TRUE));
})


test_that("Vertex-indexing bug is fixed: fsaverage_coords correspond to fsaverage_vertices", {
  fs_home = Sys.getenv("FREESURFER_HOME");
  if(nchar(fs_home) == 0L || ! dir.exists(fs_home)) {
    testthat::skip("FreeSurfer not available: FREESURFER_HOME not set or invalid.");
  }

  mni_coord_in_cortex = c(60.0, 0.0, 10.0);
  res = mni152_coords_to_fsaverage(mni_coord_in_cortex, surface = "white");

  # The returned fsaverage_coords should match the surface vertex coordinates at the returned vertex index.
  # If the bug were present, fsaverage_coords would be for vertex 1 instead of vertex 9092.
  fsaverage_surf_file = file.path(fs_home, "subjects", "fsaverage", "surf", "rh.white");
  if(! file.exists(fsaverage_surf_file)) {
    testthat::skip(sprintf("fsaverage surface file not found at '%s'.", fsaverage_surf_file));
  }
  rh_surf = freesurferformats::read.fs.surface(fsaverage_surf_file);
  expected_coords = rh_surf$vertices[res$fsaverage_vertices, ];

  testthat::expect_equal(res$hemi, c("rh"));
  testthat::expect_equal(as.vector(res$fsaverage_coords), as.vector(expected_coords), tolerance = 1e-4);
})


test_that("Colin27 coords can be mapped to fsaverage space via txt fallback", {
  fs_home = Sys.getenv("FREESURFER_HOME");
  if(nchar(fs_home) == 0L || ! dir.exists(fs_home)) {
    testthat::skip("FreeSurfer not available: FREESURFER_HOME not set or invalid.");
  }

  # Use a coordinate known to be in the Colin27 cortex. The coordinate 60,0,10 in MNI152
  # maps to fsaverage vertex 9092. In Colin27 space the mapping will differ.
  colin_coord_in_cortex = c(60.0, 0.0, 10.0);
  colin_coord_outside_cortex = c(200.0, 200.0, 200.0);  # far outside head

  res_in_cortex = colin27_coords_to_fsaverage(colin_coord_in_cortex, surface = "white");
  res_outside_cortex = colin27_coords_to_fsaverage(colin_coord_outside_cortex, surface = "white");

  # The in-cortex coord should return a valid vertex index (not NaN)
  testthat::expect_false(is.na(res_in_cortex$fsaverage_vertices));
  testthat::expect_false(is.nan(res_in_cortex$fsaverage_vertices));
  testthat::expect_true(res_in_cortex$fsaverage_vertices > 0L);
  testthat::expect_true(res_in_cortex$hemi %in% c("lh", "rh"));

  # The outside-cortex coord should return NaN
  testthat::expect_true(is.nan(res_outside_cortex$fsaverage_vertices));

  # Check return list structure has the renamed fields
  testthat::expect_true("query_coords" %in% names(res_in_cortex));
  testthat::expect_true("query_voxels" %in% names(res_in_cortex));
  testthat::expect_false("query_mni_coords" %in% names(res_in_cortex));
})


test_that("MNI152 coords can be mapped to Colin27 coords via convenience function", {
  fs_home = Sys.getenv("FREESURFER_HOME");
  if(nchar(fs_home) == 0L || ! dir.exists(fs_home)) {
    testthat::skip("FreeSurfer not available: FREESURFER_HOME not set or invalid.");
  }

  mni_coord_in_cortex = c(60.0, 0.0, 10.0);

  c27_coords = mni152_coords_to_colin27_coords(mni_coord_in_cortex, surface = "white");

  testthat::expect_true(is.matrix(c27_coords));
  testthat::expect_equal(ncol(c27_coords), 3L);
  testthat::expect_equal(nrow(c27_coords), 1L);
  testthat::expect_false(any(is.na(c27_coords)));

  # Test with simplify=TRUE
  c27_vec = mni152_coords_to_colin27_coords(mni_coord_in_cortex, surface = "white", simplify = TRUE);
  testthat::expect_true(is.vector(c27_vec));
  testthat::expect_equal(length(c27_vec), 3L);

  # Test that coordinates outside cortex return NA
  mni_outside = c(0.0, 0.0, 0.0);
  c27_outside = mni152_coords_to_colin27_coords(mni_outside, surface = "white");
  testthat::expect_true(all(is.na(c27_outside)));
})


test_that("Colin27 coords can be mapped to MNI152 coords via convenience function", {
  fs_home = Sys.getenv("FREESURFER_HOME");
  if(nchar(fs_home) == 0L || ! dir.exists(fs_home)) {
    testthat::skip("FreeSurfer not available: FREESURFER_HOME not set or invalid.");
  }

  colin_coord_in_cortex = c(60.0, 0.0, 10.0);

  mni_coords = colin27_coords_to_mni152_coords(colin_coord_in_cortex, surface = "white");

  testthat::expect_true(is.matrix(mni_coords));
  testthat::expect_equal(ncol(mni_coords), 3L);
  testthat::expect_equal(nrow(mni_coords), 1L);
  testthat::expect_false(any(is.na(mni_coords)));

  # Test with simplify=TRUE
  mni_vec = colin27_coords_to_mni152_coords(colin_coord_in_cortex, surface = "white", simplify = TRUE);
  testthat::expect_true(is.vector(mni_vec));
  testthat::expect_equal(length(mni_vec), 3L);

  # Test that coordinates outside cortex return NA
  c27_outside = c(200.0, 200.0, 200.0);  # far outside head
  mni_outside = colin27_coords_to_mni152_coords(c27_outside, surface = "white");
  testthat::expect_true(all(is.na(mni_outside)));
})


test_that("vol_coords_to_fsaverage works with RF_M3Z via txt fallback", {
  fs_home = Sys.getenv("FREESURFER_HOME");
  if(nchar(fs_home) == 0L || ! dir.exists(fs_home)) {
    testthat::skip("FreeSurfer not available: FREESURFER_HOME not set or invalid.");
  }

  mni_coord_in_cortex = c(60.0, 0.0, 10.0);

  # RF_M3Z has no .mgz file, so this exercises the .txt fallback path.
  res = vol_coords_to_fsaverage(mni_coord_in_cortex, surface = "white", rf_type = "RF_M3Z", template_type = "FSL_MNI152", silent = FALSE);

  testthat::expect_false(is.na(res$fsaverage_vertices));
  testthat::expect_false(is.nan(res$fsaverage_vertices));
  testthat::expect_true(res$fsaverage_vertices > 0L);
  testthat::expect_true(res$hemi %in% c("lh", "rh"));
  testthat::expect_equal(ncol(res$fsaverage_coords), 3L);

  # The RF_M3Z result should be close to the RF_ANTs result (same coordinate should map nearby).
  res_ants = mni152_coords_to_fsaverage(mni_coord_in_cortex, surface = "white");
  # Vertex indices may differ between methods, but both should return valid vertices.
  testthat::expect_true(res_ants$fsaverage_vertices > 0L);
})


