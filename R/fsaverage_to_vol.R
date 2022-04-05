# functions for projecting fsaverage(7/6/5) per-vertex data to MNI152/Colin27 volumes.


# lh_input = rh_input = rnorm(163842L, 3.0, 0.2); target_space="FSL_MNI152"; rf_type='RF_ANTs'; interp='linear'; out_type='mgz'; out_dir=NULL; fsaverage_path=NULL;


#' @title Project or map per-vertex values from the fsaverage surface to the cortex voxels of an MNI volume.
#'
#' @description Applies the Wu et al. regfusion method to obtain MNI volume coordinates, then interpolates values.
#'
#' @param lh_input numerical vector of per-vertex data for the left hemisphere of the template subject. Must contain 163842 values for the \code{fsaverage} template. Input for fsaverage6 (40962 values) or fsaverage5 (10242 values) can also be used and will be upsampled using \code{\link[haze]{nn_interpolate_kdtree}}. Automatic upsamping is only supported with \code{interp='linear'}.
#'
#' @param rh_input numerical vector of per-vertex data for the right hemisphere of the template subject. Must contain 163842 values for the \code{fsaverage} template. See description for \code{lh_input} for more details.
#'
#' @param target_space character string, the target template or the space that your output volume should be in. One of 'FSL_MNI152' or  'SPM_Colin27'.
#'
#' @param rf_type the \code{regfusion} type to use, one of 'RF_ANTs' or 'RF_M3Z'.
#'
#' @param interp interpolation method, currently only 'linear' and 'nearest' are supported. The performance of the 'linear' method is currently quite bad, and it will be rewritten in \code{C++} when I find the time.
#'
#' @param out_type character string, the format of the output files. One of the following: 'mgz' or 'mgh' for FreeSurfer MGZ/MGH format, 'nii' for NIFTI v1 format. Ignored unless out_dir is not NULL.
#'
#' @param out_dir optional character string, the path to a writable output directory to which the output should be written as volume files. If \code{NULL}, no data is written to files. If out_dir is not \code{NULL}, the return value additionally contains the following keys: 'out_file' and the output file format at key 'out_format', (and 'out_file_seg'/'out_format_seg' for the respective \code{seg} versions).
#'
#' @param fsaverage_path character string or NULL, the file system path to the \code{fsaverage} directory (NOT including the 'fsaverage' dir itself). If \code{NULL}, defaults to the return value of \code{fsbrain::fsaverage.path()} on the system. This path is used to read the spherical surface (both hemisphere meshes) of the template subject.
#'
#' @param silent logical, whether to suppress status messages
#'
#' @return named list with keys 'projected' and 'projected_seg', each of which holds an \code{fs.volume} instance, its 'data' key holds a \code{256x256x256} array with the projected data. The data in 'projected_seg' is identical to the data in 'projected', with the exception that data values originating from the right hemisphere have been incremented by 1000. See \code{out_dir} parameter to easily write results to files. If out_dir is not \code{NULL}, the return value additionally contains the following keys: 'out_file' and the output file format at key 'out_format'.
#'
#' @author Tim Schäfer for the R version, Wu Jianxiao and CBIG for the original Matlab version.
#'
#' @note This function requires the packages 'fsbrain' and 'haze', which are optional dependencies. The package 'fsbrain' can be installed from CRAN. For 'haze', see \code{https://github.com/dfsp-spirit/haze}.
#' @note This function is quite expensive computationally, especially when using \code{interp = 'linear'}.
#'
#' @importFrom freesurferformats write.fs.morph read.fs.volume read.fs.mgh
#' @importFrom pracma interp1
#'
#' @examples
#' \dontrun{
#'   lh_input = rnorm(163842L, 3.0, 0.2);
#'   rh_input = rnorm(163842L, 3.0, 0.2);
#'   res = fsaverage_to_vol(lh_input, rh_input, "FSL_MNI152");
#' }
#'
#' @export
fsaverage_to_vol <- function(lh_input, rh_input, target_space="FSL_MNI152", rf_type='RF_ANTs', interp='linear', out_type='mgz', out_dir=NULL, fsaverage_path=NULL, silent = TRUE) {

  rh_seg_start = 1000; # offset to identify right hemisphere projected values in whole brain volume.

  if(!(target_space %in% c("FSL_MNI152", "SPM_Colin27"))) {
    stop("Parameter 'target_space' must be one of 'FSL_MNI152' or 'SPM_Colin27'.");
  }

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
        stop(sprintf("The lengths of the lh_input and rh_input vectors must be identical for all supported templates (FreeSurfer fsaverage templates), but lengths %d and %d differ.\n", length(lh_input), length(rh_input)));
      }

      template_meshes_surface = fsbrain::subject.surface(fsaverage_path, template_subject, surface = "sphere");

      # Check input data length. May require up-samling if from fsaverage6, or fsaverage5.
      if(length(lh_input) != num_template_vertices_per_hemi) {
        if(interp != "linear") {
          stop("Automatic upsampling is only supported for 'interp'='linear'.");
        }
        if(length(lh_input) == 40962L) {
          # Automatic up-sampling of input data from fsaverage6 mesh.
          if(! silent) {
            cat(sprintf("Upsampling fsaverage6 per-vertex data to fsaverage mesh.\n"));
          }
          template_orig_meshes = fsbrain::subject.surface(fsaverage_path, "fsaverage6", surface = "sphere");
          lh_input = haze::nn_interpolate_kdtree(template_meshes_surface$lh$vertices, template_orig_meshes$lh, lh_input);
          rh_input = haze::nn_interpolate_kdtree(template_meshes_surface$lh$vertices, template_orig_meshes$lh, rh_input);
        } else if(length(lh_input) == 10242L) {
          if(! silent) {
            cat(sprintf("Upsampling fsaverage5 per-vertex data to fsaverage mesh.\n"));
          }
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

      if(! is.null(out_dir)) {
        if(! dir.exists(out_dir)) {
          dir.create(out_dir, recursive = FALSE);
        }
      }

      lh_map_file = regfusionr:::get_data_file(sprintf("allSub_fsaverage_to_%s_FS4.5.0_%s_avgMapping.lh.mgz", target_space, rf_type), subdir = "coordmap");
      rh_map_file = regfusionr:::get_data_file(sprintf("allSub_fsaverage_to_%s_FS4.5.0_%s_avgMapping.rh.mgz", target_space, rf_type), subdir = "coordmap");
      lh_coord = freesurferformats::read.fs.mgh(lh_map_file, with_header = FALSE, drop_empty_dims = TRUE); # 3x16777216 matrix
      rh_coord = freesurferformats::read.fs.mgh(rh_map_file, with_header = FALSE, drop_empty_dims = TRUE); # 3x16777216 matrix

      # Create a mask that allows us to ignore background voxels which do not map to the surface.
      lh_mask = which(colSums(lh_coord) == 0L); # 1x16777216 matrix (or vector if dropped), logical
      rh_mask = which(colSums(rh_coord) == 0L); # 1x16777216 matrix (or vector if dropped), logical
      num_voxels = length(lh_mask);

      out = list();

      projected_vol_data = list();
      if(interp == 'linear') { # for continuous data like thickness
        if(! silent) {
          cat(sprintf("Using 'linear' interpolation, suitable for continuous per-vertex data, for %d lh values.\n", length(lh_mask)));
        }
        projected_vol_data$lh = rep(0.0, num_voxels);
        projected_vol_data$lh[lh_mask] = haze::linear_interpolate_kdtree(t(lh_coord[, lh_mask]), template_meshes_surface$lh, lh_input)$interp_values;
        if(! silent) {
          cat(sprintf("Using 'linear' interpolation, suitable for continuous per-vertex data, for %d rh values.\n", length(rh_mask)));
        }
        projected_vol_data$rh = rep(0.0, num_voxels);
        projected_vol_data$rh[rh_mask] = haze::linear_interpolate_kdtree(t(rh_coord[, rh_mask]), template_meshes_surface$rh, rh_input)$interp_values;
      } else if(interp == 'nearest') { # for labels
        if(! silent) {
          cat(sprintf("Using 'nearest' interpolation, suitable for label data (integers), for %d lh values.\n", length(lh_mask)));
        }
        lh_vertex = rep(0.0, num_voxels);
        lh_vertex[lh_mask] = haze::find_nv_kdtree(t(lh_coord[, lh_mask]), template_meshes_surface$lh)$index;
        projected_vol_data$lh = rep(0.0, num_voxels);
        projected_vol_data$lh[lh_mask] = pracma::interp1(seq.int(length(lh_input)), lh_input, lh_vertex[lh_mask], method = "nearest");
        if(! silent) {
          cat(sprintf("Using 'nearest' interpolation, suitable for label data (integers), for %d rh values.\n", length(rh_mask)));
        }
        rh_vertex = rep(0.0, num_voxels);
        rh_vertex[rh_mask] = haze::find_nv_kdtree(t(rh_coord[, rh_mask]), template_meshes_surface$rh)$index;
        projected_vol_data$rh = rep(0.0, num_voxels);
        projected_vol_data$rh[rh_mask] = pracma::interp1(seq.int(length(rh_input)), rh_input, rh_vertex[rh_mask], method = "nearest");
      } else {
        stop("Currently the only supported interpolation methods are 'linear' and 'nearest'.");
      }

      if(! silent) {
        cat(sprintf("Applying volume mask and merging projected hemisphere data into a single volume.\n"));
      }

      # Load cortex mask volume file and apply it. Data projected to non-cortical brain regions is useless/wrong, so we
      # remove it with a cortical volume mask.
      cortex_mask_file_volume = regfusionr:::get_data_file(sprintf("%s_FS4.5.0_cortex_estimate.nii.gz", target_space), subdir = "coordmap");
      cortex_mask_fs_volume = freesurferformats::read.fs.volume(cortex_mask_file_volume, with_header = TRUE, drop_empty_dims = TRUE);
      cortex_mask_volume = cortex_mask_fs_volume$data;
      #cortex_label_surface = fsbrain::subject.label(fsaverage_path, template_subject, label = "cortex", hemi = "both");
      # Apply loaded vol mask.
      projected_vol_data$lh[cortex_mask_volume==0L] = 0L;
      projected_vol_data$rh[cortex_mask_volume==0L] = 0L;

      # Combine results of the two hemispheres
      projected = cortex_mask_fs_volume;
      projected$data = array(data = (projected_vol_data$lh + projected_vol_data$rh), dim = dim(cortex_mask_volume)); # reshape 1x16777216 vector to 256x256x256 array.

      # Create 2nd version with RH values incremented by 1000 offset.
      projected_vol_data_rh_shifted = projected_vol_data$rh + rh_seg_start; # TODO: this is wrong, it increments all values. we only want to increment the rh values, not the background.
      projected_seg = cortex_mask_fs_volume;
      projected_seg$data = array(data = projected_vol_data$lh + projected_vol_data_rh_shifted, dim = dim(cortex_mask_volume)); # same as above, but RH values have been incremented.

      out$projected = projected;
      out$projected_seg = projected_seg;
      if(! is.null(out_dir)) { # also write the data to disk in requested format.
        out_file = file.path(out_dir, sprintf("projected_%s_to_%s.%s", template_subject, target_space, out_type));
        freesurferformats::write.fs.morph(out_file, projected);
        out$out_file = out_file;
        out$out_format = out_type;
        out_file_seg = file.path(out_dir, sprintf("projectedseg_%s_to_%s.%s", template_subject, target_space, out_type));
        freesurferformats::write.fs.morph(out_file_seg, projected_seg);
        out$out_file_seg = out_file_seg;
        out$out_format_seg = out_type;
      }
      return(out);
    } else {
      stop("This functionality requires the 'haze' package. Please install it, see 'https://github.com/dfsp-spirit/haze' for details.");
    }
  } else {
    stop("This functionality requires the 'fsbrain' package. Please install it.");
  }
}


