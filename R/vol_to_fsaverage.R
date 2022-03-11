# Main functions for vol to surf mapping for an input 3D/4D image.


#' @title Project or map values from MNI volume to fsaverage surface.
#'
#' @description Applies the Wu et al. regfusion method to obtain surface coords, then interpolates values.
#'
#' @param input_img 3D or 4D NIFTI or MGZ image instance of type \code{fs.volume}, or a character string that will be interpreted as a file system path to such a volume that should be loaded with \code{freesurferformats::read.fs.volume}. If 4D, the 4th dimension is considered the time/subject dimension.
#'
#' @param template_type character string, the source template or the space that your input image is in. One of 'MNI152_orig', 'Colin27_orig', 'MNI152_norm', 'Colin27_norm'.
#'
#' @param rf_type the regfusion type to use, one of 'RF_ANTs' or 'RF_M3Z'.
#'
#' @param interp interpolation method, currently only 'linear' is supported.
#'
#' @param out_type character string, the format of the output files. One of the following: 'curv' for FreeSurfer curv format, 'mgz' for FreeSurfer MGZ format, 'gii' for GIFTI format.
#'
#' @param out_dir character string, the path to a writable output directory. If \code{NULL}, the returned named list contains the projected data (instead of the path of the file it was written to) at the keys 'lh' and 'rh', and the parameter 'out_type' is ignored.
#'
#' @return named list of 2 character strings, the output files (for the 2 hemispheres) at keys 'lh' and 'rh'. See the documentation for parameter 'out_dir' if you want the data in R instead.
#'
#' @author Tim Schäfer for the R version, Wu Jianxiao and CBIG for the original Matlab version.
#'
#' @importFrom data.table fread
#'
#' @export
vol_to_fsaverage <- function(input_img, template_type, rf_type='RF_ANTs', interp='linear', out_type='curv', out_dir=".") {

  if(! freesurferformats::is.fs.volume(input_img)) {
    if(is.character(input_img)) {
      input_img = freesurferformats::read.fs.volume(input_img, with_header = TRUE);
    } else {
      stop("Parameter 'input_img' must be fs.volume instance or a valid path to a volume file.");
    }
  }

  check_rf_and_template(template_type = template_type, rf_type = rf_type);

  valid_out_types = c('curv', 'mgz', 'gii');
  if(! (out_type %in% valid_out_types)) {
    stop(sprintf("Parameter 'out_type' must be one of: %s.", paste(valid_out_types(), collapse=", ")));
  }

  mapping = sprintf(".avgMapping_allSub_%s_%s_to_fsaverage.txt", rf_type, template_type);

  if(! is.null(out_dir)) {
    if(! dir.exists(out_dir)) {
      dir.create(out_dir, recursive = FALSE);
    }
  }

  out = list('lh' = NULL, 'rh' = NULL);
  for (hemi in c('lh', 'rh')) {
    mapping_file = system.file("extdata", "mappings", sprintf("%s%s", hemi, mapping), package = "regfusionr", mustWork = TRUE);
    ras = t(as.matrix(data.table::fread(mapping_file, nrows = 3, header = FALSE)));
    affine = freesurferformats::mghheader.ras2vox(input_img$header);
    projected = project_data(input_img$data, affine, ras, interp);

    if(is.null(out_dir)) {
      out[[hemi]] = projected;
    } else {
      if(is.array(projected) && out_type %in% c("curv", "gii")) {
        stop("The 'curv' and 'gii' output formats are not supported for 4D input data, try 'mgz' instead.");
        # We could write one output curv/gifti file per frame, but I guess simply using MGZ is better.
      }
      out_file = file.path(out_dir, sprintf("%s%s.%s", hemi, mapping, out_type));
      freesurferformats::write.fs.morph(out_file, projected);
      out[[hemi]] = out_file;
    }
  }

  return(out);
}


#' @title Perform coord transform and trilinear data interpolation.
#'
#' @param data 3D or 3D source image values
#'
#' @param affine 4x4 double matrix, the ras2vox transformation matrix
#'
#' @param ras nx3 double matrix, the source RAS coordinates
#'
#' @param interp character string, the interpolation mode. Only 'linear' is currently supported. We should support 'nearest' as well.
#'
#' @return a numerical vector for 3D input data or a numerical 4D array for 4D input data, the data values interpolated at the RAS coordinates. For the 4D input case, the 2D output is encoded in a 4D array with two dimensions of length 1 as follows: The first dimension contains the projected and interpolated per-vertex data values for a single frame, dimensions 2 and 3 are both of length 1, and the fourth dimensions is the timepoint/subject dimension (the frames).
#'
#' @keywords internal
project_data <- function(data, affine, ras, interp='linear') {

  if(requireNamespace("oce", quietly = TRUE)) {

    supported_interp = c("linear");
    if(!(interp %in% supported_interp)) {
      stop("Parameter 'interp' must be 'linear', others not supported yet.");
    }

    if(is.null(data)) {
      stop("Parameter 'data' must not be NULL (expected numerical array).");
    }

    check_affine(affine);

    coords = freesurferformats::doapply.transform.mtx(ras, affine);

    data = drop(data); # remove empty dimensions (typically the last one)

    if(length(dim(data)) == 3L) {
      nvols = 1L;
      x = seq_len(dim(data)[1]);
      y = seq_len(dim(data)[2]);
      z = seq_len(dim(data)[3]);
      proj_data = oce::approx3d(x, y, z, data, coords[,1], coords[,2], coords[,3]);
    } else if (length(dim(data)) == 4L) {
      nvols = dim(data)[4];
      proj_data = array(data = rep(NA, (nrow(ras)*nvols)), dim = c(nrow(ras), 1L, 1L, nvols));
      for(vol_idx in seq_len(nvols)) {
        x = seq_len(dim(data)[1]);
        y = seq_len(dim(data)[2]);
        z = seq_len(dim(data)[3]);
        approx_dta = oce::approx3d(x, y, z, data[,,,vol_idx], coords[,1], coords[,2], coords[,3]);
        proj_data[,1L,1L,vol_idx] = approx_dta;
      }
    } else {
      stop("Only 3D or 4D input data supported.");
    }
    return(proj_data);
  } else {
    stop("The 'oce' package must be installed to use this functionality. See https://github.com/dfsp-spirit/regfusionr for installation instructions.");
  }
}


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
#' @param out_dir character string, the path to a writable output directory. If \code{NULL}, the returned named list contains the projected data (instead of the path of the file it was written to) at key 'out_data', and the parameter 'out_type' is ignored.
#'
#' @param fsaverage_path character string or NULL, the file system path to the fsaverage directory (NOT including the 'fsaverage' dir itself). If \code{NULL}, defaults to the return value of \code{fsbrain::fsaverage.path()} on the system. This path is used to read the spherical surface (both hemisphere meshes) of the template subject.
#'
#' @return named list of character strings, the output file is at keys 'out_file' and the output file format at key 'out_format'. See the documentation for parameter 'out_dir' if you want the data in R instead.
#'
#' @author Tim Schäfer for the R version, Wu Jianxiao and CBIG for the original Matlab version.
#'
#' @note THIS FUNCTION IS CURRENTLY WORK-IN-PROGRESS AND NOT PART OF THE OFFICIAL API, DO NOT USE IT OR BE PREPARED FOR UN-ANNOUNCED BREAKING CHANGES ANY TIME.
#'
#' @keywords internal
fsaverage_to_vol <- function(lh_input, rh_input, template_type, rf_type='RF_ANTs', interp='linear', out_type='mgz', out_dir=".", fsaverage_path=NULL) {

  template_subject = "fsaverage";
  if(is.null(fsaverage_path)) {
    fsaverage_path = fsbrain::fsaverage.path(allow_fetch = TRUE);
  }

  num_template_vertices = 327684L; # for fsaverage
  num_template_vertices_per_hemi = num_template_vertices / 2L;

  # We currently do not support up-sampling from fsaverage6/5/..., the user has to do that manually using FreeSurfer tools at this time if the data were computed on a down-sampled surface.
  if((! is.vector(lh_input)) | (! is.numeric(lh_input))) {
    stop(sprintf("Parameter 'lh_input' must be a numerical vector containing %d values for the fsaverage left hemisphere vertices.\n", num_template_vertices_per_hemi));
  }
  if((! is.vector(rh_input)) | (! is.numeric(rh_input))) {
    stop(sprintf("Parameter 'rh_input' must be a numerical vector containing %d values for the fsaverage right hemisphere vertices.\n", num_template_vertices_per_hemi));
  }

  if(length(lh_input) != length(rh_input)) {
    stop(sprintf("The lengths of the lh_input and rh_input vectors must be identical for all supported templates (FreeSurfer fsaverage templates), but lengths $d and %d differ.\n", length(lh_input), length(rh_input)));
  }

  template_meshes_surface = fsbrain::subject.surface(fsaverage_path, template_subject, surface = "sphere");

  if(length(lh_input) != num_template_vertices_per_hemi) {
    if(length(lh_input) == 40962L) {
      message("Automatic up-sampling of input data from fsaverage6 mesh not supported yet.");
      template_orig_meshes = fsbrain::subject.surface(fsaverage_path, "fsaverage6", surface = "sphere");
      lh_input = haze::nn_interpolate_kdtree(template_meshes_surface$lh$vertices, template_orig_meshes$lh, lh_input);
      rh_input = haze::nn_interpolate_kdtree(template_meshes_surface$lh$vertices, template_orig_meshes$lh, rh_input);
    } else if(length(lh_input) == 10242L) {
      message("Automatic up-sampling of input data from fsaverage5 mesh not supported yet.");
      template_orig_meshes = fsbrain::subject.surface(fsaverage_path, "fsaverage5", surface = "sphere");
      lh_input = haze::nn_interpolate_kdtree(template_meshes_surface$lh$vertices, template_orig_meshes$lh, lh_input);
      rh_input = haze::nn_interpolate_kdtree(template_meshes_surface$lh$vertices, template_orig_meshes$lh, rh_input);
    } else {
      stop(sprintf("Unsupported number of input vertices: %d. Value must match vertex count for one of the templates fsaverage (163842), fsaverage6 (40962), or fsaverage5 (10242).\n", length(lh_input)));
    }

    stop(sprintf("The input vectors lh_input and rh_input must contain exactly %d values each.\n", num_template_vertices_per_hemi));
  }


  check_rf_and_template(template_type = template_type, rf_type = rf_type);

  valid_out_types = c('mgh', 'mgz', 'nii');
  if(! (out_type %in% valid_out_types)) {
    stop(sprintf("Parameter 'out_type' must be one of: %s.", paste(valid_out_types(), collapse=", ")));
  }

  mapping = "FSL_MNI152";    # One of "FSL_MNI152" or "SPM_Colin27".
  # Load cortex mask volume file.
  cortex_mask_file_volume = get_data_file(sprintf("%s_FS4.5.0_cortex_estimate.nii.gz", mapping), subdir = "coordmap");
  cortex_mask_volume = freesurferformats::read.fs.volume(cortex_mask_file_volume, with_header = TRUE, drop_empty_dims = TRUE);


  cortex_label_surface = fsbrain::subject.label(fsaverage_path, template_subject, label = "cortex", hemi = "both");
  #cortex_mask_surface = list();
  #cortex_mask_surface$lh = fsbrain::mask.from.labeldata.for.hemi(cortex_label_surface$lh,  num_template_vertices_per_hemi);
  #cortex_mask_surface$rh = fsbrain::mask.from.labeldata.for.hemi(cortex_label_surface$rh,  num_template_vertices_per_hemi);

  lh_map_file = get_data_file("FSL_MNI152_FS4.5.0_RF_ANTs_avgMapping.vertex.lh.mgz", subdir = "coordmap");
  rh_map_file = get_data_file("FSL_MNI152_FS4.5.0_RF_ANTs_avgMapping.vertex.rh.mgz", subdir = "coordmap");
  lh_vertex = freesurferformats::read.fs.mgh(lh_map_file, with_header = FALSE, drop_empty_dims = TRUE); # 256x256x256 array
  rh_vertex = freesurferformats::read.fs.mgh(rh_map_file, with_header = FALSE, drop_empty_dims = TRUE); # 256x256x256 array


  if(! is.null(out_dir)) {
    if(! dir.exists(out_dir)) {
      dir.create(out_dir, recursive = FALSE);
    }
  }

  out = list();

  # Find the vertices on

  projected_vol_data = list();
  projected_vol_data$lh = rep(0.0, num_template_vertices_per_hemi);
  projected_vol_data$rh = rep(0.0, num_template_vertices_per_hemi);
  if(interp == 'linear') {
    projected_vol_data$lh = haze::linear_interpolate_kdtree(template_meshes_surface$lh$vertices[cortex_label_surface$lh, ], template_meshes_surface$lh, lh_input);
    projected_vol_data$rh = haze::linear_interpolate_kdtree(template_meshes_surface$rh$vertices[cortex_label_surface$rh, ], template_meshes_surface$rh, lh_input);
  } else {
    stop("Currently the only supported interpolation method is 'linear'.");
  }

  # TODO: Apply the volume mask to the result, and combine the results of the hemispheres.
  # see https://github.com/ThomasYeoLab/CBIG/blob/master/stable_projects/registration/Wu2017_RegistrationFusion/bin/scripts_final_proj/CBIG_RF_projectfsaverage2Vol_single.m


  if(is.null(out_dir)) {
    out$out_data = projected_vol_data;
  } else {
    out_file = file.path(out_dir, sprintf("projected_%s_to_%s.%s", template_subject, mapping, out_type));
    freesurferformats::write.fs.morph(out_file);
    out$out_file = out_file;
    out$out_format = out_type;
  }
  return(out);
}

