library(spant)
library(cowplot)
library(ggplot2)
theme_set(theme_bw())

sim_dataset <- function(seq_tr, N_scans, bz_inhom_lb, bold_lb_hz, ss_spec_snr,
                        subjects, ft, out_dir)
{
  # Make a data frame containing a single row of basis signal amplitudes.
  # Metabolite values are for visual cortex listed in Bednarik et al 2015 Table 1.
  # Note Alanine and Glycine are not listed in the table and therefore set to 0.
  basis_amps <- data.frame("ala"    = 0.00, "asc"    = 0.96, "asp"   = 3.58,
                           "cr"     = 4.22, "gaba"   = 1.03, "glc"   = 0.62,
                           "gln"    = 2.79, "gly"    = 0.00, "gsh"   = 1.09,
                           "glu"    = 8.59, "gpc"    = 0.54, "ins"   = 6.08,
                           "lac"    = 1.01, "naa"    = 11.9, "naag"  = 1.32,
                           "pch"    = 0.40, "pcr"    = 3.34, "peth"  = 0.93,
                           "sins"   = 0.27, "tau"    = 1.27, "lip09" = 0.00,
                           "lip13a" = 0.00, "lip13b" = 0.00, "lip20" = 0.00,
                           "mm09"   = 4.00, "mm12"   = 4.00, "mm14"  = 4.00,
                           "mm17"   = 4.00, "mm20"   = 4.00)
  
  # Duplicate the row N_scans times to make a table of values
  basis_amps <- basis_amps[rep(1, N_scans),]
  
  # simulate two 120 second blocks of stimulation starting at 100 and 500
  # seconds
  onsets    <- c(320, 960)
  durations <- rep(320, 2)
  
  # generate a dummy mrs_data object for generating the regressors
  mrs_data_dummy <- sim_zero(dyns = N_scans) |> set_tr(seq_tr) |> 
    set_Ntrans(N_scans)
  
  bold_rf <- gen_bold_reg(onsets, durations, mrs_data = mrs_data_dummy)
  
  # simulate a typical basis for TE=28ms semi-LASER acquisition at 3T / 7T
  acq_paras <- def_acq_paras(ft = ft, N = 2048, fs = 6000)
  basis     <- sim_basis(names(basis_amps), pul_seq = seq_slaser_ideal,
                         TE1 = 0.008, TE2 = 0.011, TE3 = 0.009,
                         acq_paras = acq_paras)
  
  dir.create("basis_sets", showWarnings = FALSE)
  lcm_basis_file <- file.path("basis_sets", paste0("slaser_",
                                                   round(ft * 1e-6),
                                                   "_mhz.basis"))
  nii_basis_dir  <- file.path("basis_sets", paste0("slaser_",
                                                   round(ft * 1e-6), "_mhz"))
  write_basis(basis, lcm_basis_file)
  write_basis_niidir(basis, nii_basis_dir)
  
  # apply basis amplitudes to the basis set to generate a simulated fMRS dataset
  mrs_dyn_orig <- basis2dyn_mrs_data(basis, basis_amps, seq_tr)
  
  # broaden basis to simulate B0 inhomogeneity, apply any addition BOLD related 
  # narrowing
  bold_lb_dyn <- (1 - bold_rf$stim_bold) * bold_lb_hz
  bold_lb_dyn <- bold_lb_dyn + Mod(min(bold_lb_dyn))
  
  dir.create("regressors", showWarnings = FALSE)
  write.table(bold_lb_dyn, file.path("regressors", paste0("lorentz_lb_",
                                                 round(ft * 1e-6), "_mhz.csv")),
            col.names = FALSE, sep = ",")
  
  mrs_dyn <- mrs_dyn_orig |> lb(bz_inhom_lb) |> lb(bold_lb_dyn, 0) 
  
  # duplicate the data to generate multiple subjects with different noise
  # samples
  mrs_dyn_list <- rep(list(mrs_dyn), subjects)
  
  # add noise
  mrs_dyn_list <- mrs_dyn_list |> add_noise_spec_snr(ss_spec_snr)
  
  # export to BIDS structure
  mr_data2bids(mrs_dyn_list, "svs", out_dir, skip_existing = FALSE)
}

# random number generator seed used for noise samples
set.seed(1)

n_subjects <- 128

# 3T simulation
sim_dataset(seq_tr = 2.5, N_scans = 640, bz_inhom_lb = 6, bold_lb_hz = 0.2,
            ss_spec_snr = 45, subjects = n_subjects, ft = 123e6,
            out_dir = file.path("simulated_data", "fmrs_bold_3T"))

# 7T simulation
sim_dataset(seq_tr = 5.0, N_scans = 320, bz_inhom_lb = 11, bold_lb_hz = 0.46,
            ss_spec_snr = 80, subjects = n_subjects, ft = 297e6,
            out_dir = file.path("simulated_data", "fmrs_bold_7T"))

# example spectra plot
path_3t <- file.path("simulated_data", "fmrs_bold_3T", "sub-001", "mrs",
                     "sub-001_svs.nii.gz")

eg_mrs_3t <- read_mrs(path_3t)

path_7t <- file.path("simulated_data", "fmrs_bold_7T", "sub-001", "mrs",
                     "sub-001_svs.nii.gz")

eg_mrs_7t <- read_mrs(path_7t)

part_a <- \() eg_mrs_3t |> get_dyns(1) |> zf() |> plot(xlim = c(4, 0.5))
part_b <- \() eg_mrs_3t |> mean_dyns() |> zf() |> plot(xlim = c(4, 0.5))
part_c <- \() eg_mrs_7t |> get_dyns(1) |> zf() |> plot(xlim = c(4, 0.5))
part_d <- \() eg_mrs_7t |> mean_dyns() |> zf() |> plot(xlim = c(4, 0.5))

png("eg_spec.png", width = 1000, height = 800, res = 120)
plot_grid(part_a, part_b, part_c, part_d, labels = c('A', 'B', 'C', 'D'),
          label_size = 12)
dev.off()

# BOLD LW plot
bold_lw_3t <- read.csv(file.path("regressors", "lorentz_lb_123_mhz.csv"),
                       header = FALSE)
bold_lw_7t <- read.csv(file.path("regressors", "lorentz_lb_297_mhz.csv"),
                       header = FALSE)

time_3t <- seq(from = 0, by = 2.5, length.out = 640)
time_7t <- seq(from = 0, by = 5, length.out = 320)

bold_df <- data.frame(time = c(time_3t, time_7t),
                      lw = c(bold_lw_3t$V2, bold_lw_7t$V2),
                      field_strength = c(rep("3T", 640), rep("7T", 320)))

bold_df |> ggplot(aes(x = time, y = lw, col = field_strength)) + geom_line() +
           ylab("Applied Lorentzian linebroadening (Hz)") + xlab("Time (s)") + 
           labs(color = "Field strength")

ggsave("BOLD_LB.png", width = 6, height = 3)