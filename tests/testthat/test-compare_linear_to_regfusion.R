
test_that("We can compare the linear FreeSurfer method versus the regfusion approach", {

  testthat::skip_on_cran();

  fs_home = Sys.getenv("FREESURFER_HOME");
  if(nchar(fs_home) == 0L || ! dir.exists(fs_home)) {
    testthat::skip("FreeSurfer not available: FREESURFER_HOME not set or invalid.");
  }

  lh_orig_file = file.path(fs_home, "subjects", "fsaverage", "surf", "lh.orig");
  rh_orig_file = file.path(fs_home, "subjects", "fsaverage", "surf", "rh.orig");
  if(! file.exists(lh_orig_file)) {
    testthat::skip(sprintf("Could not find fsaverage lh.orig surface file at '%s'.\n", lh_orig_file));
  }
  if(! file.exists(rh_orig_file)) {
    testthat::skip(sprintf("Could not find fsaverage rh.orig surface file at '%s'.\n", rh_orig_file));
  }
  lh_orig = freesurferformats::read.fs.surface(lh_orig_file);
  rh_orig = freesurferformats::read.fs.surface(rh_orig_file);

  num_vertices = nrow(lh_orig$vertices);          # 163842, same for the right hemi.
  query_fsaverage_vertices = seq.int(num_vertices);

  mni_coords_linear_lh = linear_fsaverage_coords_to_MNI152_coords(lh_orig$vertices);
  mni_coords_regfusionr_lh = fsaverage_vertices_to_mni152_coords(vertices = query_fsaverage_vertices, hemis = rep("lh", num_vertices));
  mni_coords_linear_rh = linear_fsaverage_coords_to_MNI152_coords(rh_orig$vertices);
  mni_coords_regfusionr_rh = fsaverage_vertices_to_mni152_coords(vertices = query_fsaverage_vertices, hemis = rep("rh", num_vertices));

  dist_diff_lh = rep(NA, num_vertices);
  dist_diff_rh = rep(NA, num_vertices);
  for(row_idx in seq.int(num_vertices)) {
    dist_diff_lh[row_idx] = sqrt(sum((mni_coords_linear_lh[row_idx, ] - mni_coords_regfusionr_lh[row_idx, ]) ^ 2));
    dist_diff_rh[row_idx] = sqrt(sum((mni_coords_linear_rh[row_idx, ] - mni_coords_regfusionr_rh[row_idx, ]) ^ 2));
  }

  do_plot = FALSE;
  if(do_plot) {
    #if(requireNamespace("fsbrain", quietly = TRUE)) {
      #cm = fsbrain::vis.data.on.fsaverage(morph_data_lh = dist_diff_lh, morph_data_rh = dist_diff_rh, draw_colorbar = TRUE);
      #fsbrain::export(cm, colorbar_legend = "Difference between regfusion and linear method by vertex [mm].", output_img = "~/regfusionr_vs_linear.png");
    #}
  }

  testthat::expect_equal(1L ,1L); # Without an expect, the test will be skipped.
})
