# functions for projecting fsaverage(7/6/5) per-vertex data to MNI152/Colin27 volumes.


# lh_input = rh_input = rnorm(163842L, 3.0, 0.2); template_type="MNI152_orig"; rf_type='RF_ANTs'; interp='linear'; out_type='mgz'; out_dir=NULL; fsaverage_path=NULL;


#' @title Project or map per-vertex values from the fsaverage surface to the cortex voxels of an MNI volume.
#'
#' @description Applies the Wu et al. regfusion method to obtain MNI volume coordinates, then interpolates values.
#'
#' @param lh_input numerical vector of per-vertex data for the left hemisphere of the template subject. Must 163842 values for the \code{fsaverage} template.
#'
#' @param rh_input numerical vector of per-vertex data for the right hemisphere of the template subject. Must 163842 values for the \code{fsaverage} template.
#'
#' @param template_type character string, the target template or the space that your output volume should be in. One of 'MNI152_orig', 'Colin27_orig', 'MNI152_norm', 'Colin27_norm'.
#'
#' @param rf_type the regfusion type to use, one of 'RF_ANTs' or 'RF_M3Z'.
#'
#' @param interp interpolation method, currently only 'linear' is supported.
#'
#' @param out_type character string, the format of the output files. One of the following: 'mgz' or 'mgh' for FreeSurfer MGZ/MGH format, 'nii' for NIFTI v1 format.
#'
#' @param out_dir character string, the path to a writable output directory. If \code{NULL}, the returned named list contains the projected data (instead of the path of the file it was written to) at key 'out_data', and the parameter 'out_type' is ignored. The 'out_data' is a named list with keys 'lh', 'rh' and 'both', each of which holds a \code{256x256x256} array with the data.
#'
#' @param fsaverage_path character string or NULL, the file system path to the fsaverage directory (NOT including the 'fsaverage' dir itself). If \code{NULL}, defaults to the return value of \code{fsbrain::fsaverage.path()} on the system. This path is used to read the spherical surface (both hemisphere meshes) of the template subject.
#'
#' @return see \code{out_dir} parameter. If out_dir is not \code{NULL}, the return value is instead a named list of character strings, the output file is at keys 'out_file' and the output file format at key 'out_format'. See the documentation for parameter 'out_dir' if you want the data in R instead.
#'
#' @author Tim Schäfer for the R version, Wu Jianxiao and CBIG for the original Matlab version.
#'
#' @note THIS FUNCTION IS CURRENTLY WORK-IN-PROGRESS AND NOT PART OF THE OFFICIAL API, DO NOT USE IT OR BE PREPARED FOR UN-ANNOUNCED BREAKING CHANGES ANY TIME.
#' @note This function requires the packages 'fsbrain' and 'haze', which are optional dependencies. The package 'fsbrain' can be installed from CRAN. For 'haze', see \code{https://github.com/dfsp-spirit/haze}.
#'
#' @importFrom freesurferformats write.fs.morph read.fs.volume read.fs.mgh
#'
#' @examples
#' \dontrun{
#'   lh_input = rnorm(163842L, 3.0, 0.2);
#'   rh_input = rnorm(163842L, 3.0, 0.2);
#'   res = fsaverage_to_vol(lh_input, rh_input);
#' }
#'
#' @export
fsaverage_to_vol <- function(lh_input, rh_input, template_type="MNI152_orig", rf_type='RF_ANTs', interp='linear', out_type='mgz', out_dir=NULL, fsaverage_path=NULL) {

  if(requireNamespace("fsbrain", quietly = TRUE)) {
    if(requireNamespace("haze", quietly = TRUE)) {

      template_subject = "fsaverage";
      if(is.null(fsaverage_path)) {
        fsaverage_path = fsbrain::fsaverage.path(allow_fetch = TRUE);
      }

      num_template_vertices = 327684L; # for fsaverage
      num_template_vertices_per_hemi = num_template_vertices / 2L;

      # We support automatic up-sampling of fsaverage6 and fsaverage5 data to fsaverage.
      if((! is.vector(lh_input)) | (! is.numeric(lh_input))) {
        stop(sprintf("Parameter 'lh_input' must be a numerical vector containing the per-vertex data for the fsaverage/fsaverage6/fsaverage5 left hemisphere vertices.\n"));
      }
      if((! is.vector(rh_input)) | (! is.numeric(rh_input))) {
        stop(sprintf("Parameter 'rh_input' must be a numerical vector containing the per-vertex data for the fsaverage/fsaverage6/fsaverage5 right hemisphere vertices.\n"));
      }

      if(length(lh_input) != length(rh_input)) {
        stop(sprintf("The lengths of the lh_input and rh_input vectors must be identical for all supported templates (FreeSurfer fsaverage templates), but lengths $d and %d differ.\n", length(lh_input), length(rh_input)));
      }

      template_meshes_surface = fsbrain::subject.surface(fsaverage_path, template_subject, surface = "sphere");

      if(length(lh_input) != num_template_vertices_per_hemi) {
        if(length(lh_input) == 40962L) {
          # Automatic up-sampling of input data from fsaverage6 mesh.
          template_orig_meshes = fsbrain::subject.surface(fsaverage_path, "fsaverage6", surface = "sphere");
          lh_input = haze::nn_interpolate_kdtree(template_meshes_surface$lh$vertices, template_orig_meshes$lh, lh_input);
          rh_input = haze::nn_interpolate_kdtree(template_meshes_surface$lh$vertices, template_orig_meshes$lh, rh_input);
        } else if(length(lh_input) == 10242L) {
          # Automatic up-sampling of input data from fsaverage5 mesh.
          template_orig_meshes = fsbrain::subject.surface(fsaverage_path, "fsaverage5", surface = "sphere");
          lh_input = haze::nn_interpolate_kdtree(template_meshes_surface$lh$vertices, template_orig_meshes$lh, lh_input);
          rh_input = haze::nn_interpolate_kdtree(template_meshes_surface$lh$vertices, template_orig_meshes$lh, rh_input);
        } else {
          stop(sprintf("Unsupported number of input vertices: %d. Value must match vertex count for one of the templates fsaverage (163842), fsaverage6 (40962), or fsaverage5 (10242). Please map data to an MNI305 template using recon-all/qcache.\n", length(lh_input)));
        }
      }

      valid_out_types = c('mgh', 'mgz', 'nii');
      if(! (out_type %in% valid_out_types)) {
        stop(sprintf("Parameter 'out_type' must be one of: %s.", paste(valid_out_types(), collapse=", ")));
      }

      check_rf_and_template(template_type = template_type, rf_type = rf_type);
      # Currently only MNI152 is supported.
      if(template_type != "MNI152_orig") {
        stop("Currently the only supported 'template_type' is 'MNI152_orig'.");
      }

      mapping = "FSL_MNI152";    # One of "FSL_MNI152" or "SPM_Colin27".
      # Load cortex mask volume file.
      cortex_mask_file_volume = get_data_file(sprintf("%s_FS4.5.0_cortex_estimate.nii.gz", mapping), subdir = "coordmap");
      cortex_mask_volume = freesurferformats::read.fs.volume(cortex_mask_file_volume, with_header = TRUE, drop_empty_dims = TRUE);


      cortex_label_surface = fsbrain::subject.label(fsaverage_path, template_subject, label = "cortex", hemi = "both");
      # We do not need a mask, the label is fine. So these lines are currently commented out.
      #cortex_mask_surface = list();
      #cortex_mask_surface$lh = fsbrain::mask.from.labeldata.for.hemi(cortex_label_surface$lh,  num_template_vertices_per_hemi);
      #cortex_mask_surface$rh = fsbrain::mask.from.labeldata.for.hemi(cortex_label_surface$rh,  num_template_vertices_per_hemi);

      if(! is.null(out_dir)) {
        if(! dir.exists(out_dir)) {
          dir.create(out_dir, recursive = FALSE);
        }
      }

      # The 3D arrays in the following files assign to each voxel a vertex index (integer).
      ## ALL OF THIS IS RUBBISH: we need the reverse mapping.
      #lh_map_file = get_data_file(sprintf("FSL_MNI152_FS4.5.0_%s_avgMapping.vertex.lh.mgz", rf_type), subdir = "coordmap");
      #rh_map_file = get_data_file(sprintf("FSL_MNI152_FS4.5.0_%s_avgMapping.vertex.rh.mgz", rf_type), subdir = "coordmap");
      #lh_coord = freesurferformats::read.fs.mgh(lh_map_file, with_header = FALSE, drop_empty_dims = TRUE); # 256x256x256 array
      #rh_coord = freesurferformats::read.fs.mgh(rh_map_file, with_header = FALSE, drop_empty_dims = TRUE); # 256x256x256 array

      # Binary volume masks for projection. Voxels from lh_coord and rh_coord with value (0, 0, 0) mapping will be masked out and can be ignored as they are not part of the cortex.
      ## ALL OF THIS IS RUBBISH, see the comment on the mapping above.
      #lh_mask = array(data = rep(FALSE, prod(dim(lh_coord))), dim = dim(lh_coord));
      #lh_mask[which(lh_coord != 0L, arr.ind = TRUE)] = TRUE;
      #rh_mask = array(data = rep(FALSE, prod(dim(rh_coord))), dim = dim(rh_coord));
      #rh_mask[which(rh_coord != 0L, arr.ind = TRUE)] = TRUE;


      lh_map_file = regfusionr:::get_data_file(sprintf("allSub_fsaverage_to_%s_FS4.5.0_%s_avgMapping.lh.mgz", mapping, rf_type), subdir = "coordmap");
      rh_map_file = regfusionr:::get_data_file(sprintf("allSub_fsaverage_to_%s_FS4.5.0_%s_avgMapping.rh.mgz", mapping, rf_type), subdir = "coordmap");
      lh_coord = freesurferformats::read.fs.mgh(lh_map_file, with_header = FALSE, drop_empty_dims = TRUE); # 3x16777216 matrix
      rh_coord = freesurferformats::read.fs.mgh(rh_map_file, with_header = FALSE, drop_empty_dims = TRUE); # 3x16777216 matrix

      # TODO: create lh_mask and rh_mask
      lh_mask = ; # 1x16777216 matrix (or vector if dropped), logical
      rh_mask = ; # 1x16777216 matrix (or vector if dropped), logical

      out = list();

      projected_vol_data = list();
      if(interp == 'linear') {
        lh_interp_res = haze::linear_interpolate_kdtree(template_meshes_surface$lh$vertices[cortex_label_surface$lh, ], template_meshes_surface$lh, lh_input);
        projected_vol_data$lh = lh_interp_res$interp_values;
        rh_interp_res = haze::linear_interpolate_kdtree(template_meshes_surface$rh$vertices[cortex_label_surface$rh, ], template_meshes_surface$rh, lh_input);
        projected_vol_data$rh = rh_interp_res$interp_values;
      } else if(interp == 'nearest') {
        # we could use haze::nn_interpolate_kdtree/haze::find_nv_kdtree and pracma::interp1 to implement this.
        stop("The 'nearest' method is not implemented yet.");
      } else {
        stop("Currently the only supported interpolation method is 'linear'.");
      }

      # TODO: Convert the 2D surface vector into 3D volume data.
      # TODO: Apply the volume mask to the result, and combine the results of the hemispheres.
      # see https://github.com/ThomasYeoLab/CBIG/blob/master/stable_projects/registration/Wu2017_RegistrationFusion/bin/scripts_final_proj/CBIG_RF_projectfsaverage2Vol_single.m

      # ...

      projected_vol_data$both = projected_vol_data$lh; # TODO: merge lh and rh.

      if(is.null(out_dir)) {
        out$out_data = projected_vol_data;
      } else {
        out_file = file.path(out_dir, sprintf("projected_%s_to_%s_both.%s", template_subject, mapping, out_type));
        freesurferformats::write.fs.morph(out_file, projected_vol_data$both);
        out$out_file = out_file;
        out$out_format = out_type;
      }
      return(out);
    } else {
      stop("This functionality requires the 'haze' package. Please install it, see 'https://github.com/dfsp-spirit/haze' for details.");
    }
  } else {
    stop("This functionality requires the 'fsbrain' package. Please install it.");
  }
}

