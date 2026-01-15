library(spant)
library(foreach)
library(doParallel)

n_cores <- 32

fit_run <- function(field_strength, av_scheme, dyn_basis_lb, av_name,
                    do_dyn_basis_analysis) {
  
  mrs_data_paths <- find_bids_mrs(file.path("simulated_data",
                                            paste0("fmrs_bold_",
                                                   field_strength)))$path
  
  cluster <- makeCluster(n_cores)
  registerDoParallel(cluster)
  
  foreach(n = 1:length(mrs_data_paths), .packages = "spant") %dopar% {
    
    system(paste0("echo 'Processing ", n, " of ", length(mrs_data_paths),".'"))
    
    n_lab <- sprintf("%03d", n)
    
    out_dir <- file.path("fit_results", paste0("abfit-reg_static-basis_",
                                               field_strength, "_", av_name), 
                                               n_lab)
    fit_svs(mrs_data_paths[n], dyn_av_scheme = av_scheme, TE1 = 0.008,
            TE2 = 0.011, TE3 = 0.009, TE = 0.028, pul_seq = "slaser",
            fit_method = "abfit-reg", append_basis = c("peth", "asc", "gly"),
            output_dir = out_dir)
    
    out_dir <- file.path("fit_results", paste0("lcmodel_static-basis_",
                                               field_strength, "_", av_name),
                                               n_lab)
    fit_svs(mrs_data_paths[n], dyn_av_scheme = av_scheme, TE1 = 0.008,
            TE2 = 0.011, TE3 = 0.009, TE = 0.028, pul_seq = "slaser",
            fit_method = "lcmodel", append_basis = c("peth", "asc", "gly"),
            output_dir = out_dir)
    
    if (do_dyn_basis_analysis) {
      
      out_dir <- file.path("fit_results", paste0("abfit-reg_dyn-basis_",
                                                 field_strength, "_", av_name),
                                                 n_lab)
      fit_svs(mrs_data_paths[n], dyn_av_scheme = av_scheme, TE1 = 0.008,
              TE2 = 0.011, TE3 = 0.009, TE = 0.028, pul_seq = "slaser",
              fit_method = "abfit-reg", append_basis = c("peth", "asc", "gly"),
              output_dir = out_dir, dyn_basis_lb = dyn_basis_lb)
      
      out_dir <- file.path("fit_results", paste0("lcmodel_dyn-basis_",
                                                 field_strength, "_", av_name),
                                                 n_lab)
      fit_svs(mrs_data_paths[n], dyn_av_scheme = av_scheme, TE1 = 0.008,
              TE2 = 0.011, TE3 = 0.009, TE = 0.028, pul_seq = "slaser",
              fit_method = "lcmodel", append_basis = c("peth", "asc", "gly"),
              output_dir = out_dir, dyn_basis_lb = dyn_basis_lb)
    }
  }
  
  stopCluster(cl = cluster)

}

fit_run("3T", rep(c(1, 2, 1, 2, 1),  each = 128), c(0.2, 0), "02-blocks", TRUE)
fit_run("3T", rep(1:5,  each = 128), c(0.2, 0, 0.2, 0, 0.2), "05-blocks", TRUE)
fit_run("3T", rep(1:10, each = 64),  rep(c(0.2, 0, 0.2, 0, 0.2), each = 2),
        "10-blocks", TRUE)

fit_run("7T", rep(c(1, 2, 1, 2, 1),  each = 64), c(0.46, 0), "02-blocks", TRUE)
fit_run("7T", rep(1:5,  each = 64),  c(0.46, 0, 0.46, 0, 0.46), "05-blocks",
        TRUE)
fit_run("7T", rep(1:10, each = 32),  rep(c(0.46, 0, 0.46, 0, 0.46), each = 2),
        "10-blocks", TRUE)

fit_run("3T-preproc-lb", rep(c(1, 2, 1, 2, 1), each = 128), c(0.2, 0),
        "02-blocks", FALSE)
fit_run("3T-preproc-lb", rep(1:5,  each = 128), c(0.2, 0, 0.2, 0, 0.2),
        "05-blocks", FALSE)
fit_run("3T-preproc-lb", rep(1:10, each = 64), rep(c(0.2, 0, 0.2, 0, 0.2),
                                               each = 2), "10-blocks", FALSE)

fit_run("7T-preproc-lb", rep(c(1, 2, 1, 2, 1), each = 64), c(0.46, 0),
        "02-blocks", FALSE)
fit_run("7T-preproc-lb", rep(1:5,  each = 64), c(0.46, 0, 0.46, 0, 0.46),
        "05-blocks", FALSE)
fit_run("7T-preproc-lb", rep(1:10, each = 32), rep(c(0.46, 0, 0.46, 0, 0.46),
                                               each = 2), "10-blocks", FALSE)

# combine results
fit_runs <- basename(Sys.glob("fit_results/*"))
for (n in 1:length(fit_runs)) {
  fit_svs_group_results(search_path = file.path("fit_results", fit_runs[n]), 
                    output_dir = file.path("fit_results_combined", fit_runs[n]))
}
