library(spant)

# random number generator seed used for noise samples in lb_renoise function
set.seed(1)

# 3T data
input_path_3t      <- file.path("simulated_data", "fmrs_bold_3T")
mrs_data_in_paths  <- find_bids_mrs(input_path_3t)$path
mrs_data_out_paths <- gsub("fmrs_bold_3T", "fmrs_bold_3T-preproc-lb",
                           mrs_data_in_paths)

# generate predicted BOLD lb
onsets      <- c(320, 960)
durations   <- rep(320, 2)

bold_rf_3t  <- gen_bold_reg(onsets, durations, tr = 2.5, Ndyns = 640)
bold_lb_dyn <- bold_rf_3t$stim_bold * 0.2
lb_vec_3t   <- bold_lb_dyn + Mod(min(bold_lb_dyn))

for (n in 1:length(mrs_data_out_paths)) {
  dir.create(dirname(mrs_data_out_paths[n]), recursive = TRUE)
  mrs_data <- read_mrs(mrs_data_in_paths[n])
  mrs_data_proc <- lb_renoise(mrs_data, lb_vec_3t)
  write_mrs(mrs_data_proc, mrs_data_out_paths[n], force = TRUE)
}

# 7T data
input_path_7t      <- file.path("simulated_data", "fmrs_bold_7T")
mrs_data_in_paths  <- find_bids_mrs(input_path_7t)$path
mrs_data_out_paths <- gsub("fmrs_bold_7T", "fmrs_bold_7T-preproc-lb",
                           mrs_data_in_paths)

bold_rf_7t  <- gen_bold_reg(onsets, durations, tr = 5, Ndyns = 320)
bold_lb_dyn <- bold_rf_7t$stim_bold * 0.46
lb_vec_7t   <- bold_lb_dyn + Mod(min(bold_lb_dyn))

for (n in 1:length(mrs_data_out_paths)) {
  dir.create(dirname(mrs_data_out_paths[n]), recursive = TRUE)
  mrs_data <- read_mrs(mrs_data_in_paths[n])
  mrs_data_proc <- lb_renoise(mrs_data, lb_vec_7t)
  write_mrs(mrs_data_proc, mrs_data_out_paths[n], force = TRUE)
}
