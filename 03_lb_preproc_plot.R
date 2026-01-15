library(spant)
library(tidyverse)
library(cowplot)

theme_set(theme_bw())

eg_data <- file.path("simulated_data", "fmrs_bold_7T", "sub-001", "mrs", 
                     "sub-001_svs.nii.gz")

mrs <- read_mrs(eg_data)

orig <- mrs

# generate predicted BOLD lb
onsets      <- c(320, 960)
durations   <- rep(320, 2)
bold_rf_7t  <- gen_bold_reg(onsets, durations, tr = 5, Ndyns = 320)
bold_lb_dyn <- bold_rf_7t$stim_bold * 0.46
lb_vec_7t   <- bold_lb_dyn + Mod(min(bold_lb_dyn))

lb_std <- mrs |> lb(lb_vec_7t, 0)

set.seed(1)
lb_renoise <- mrs |> lb_renoise(lb_vec_7t)

task_change <- c(64, 128, 192, 256)

a_plot <- \(x) orig |> image(xlim = c(4, 0.5), hline = task_change)
b_plot <- \(x) orig |> sub_mean_dyns() |>
  image(xlim = c(4, 0.5), hline = task_change)
c_plot <- \(x) lb_std |> sub_mean_dyns() |>
  image(xlim = c(4, 0.5), hline = task_change)
d_plot <- \(x) lb_renoise |> sub_mean_dyns() |>
  image(xlim = c(4, 0.5), hline = task_change)

png("lb_preproc.png", width = 1000, height = 1100, res = 120)
plot_grid(a_plot, b_plot, c_plot, d_plot, labels = c('A', 'B', 'C', 'D'),
          label_size = 12)
dev.off()


# orig_noise_snr <- calc_spec_snr(orig, full_output = TRUE)$snr |> as.numeric()
# lb_only_snr    <- calc_spec_snr(lb_std, full_output = TRUE)$snr |> 
#     as.numeric()
# lb_renoise_snr <- calc_spec_snr(lb_renoise, full_output = TRUE)$snr |> 
#     as.numeric()
# noise_snr_df <- data.frame(dynamic = 1:320, orig = orig_noise_snr,
#                            lb_only = lb_only_snr, lb_renoise = lb_renoise_snr)
# 
# noise_snr_long <- noise_snr_df |> pivot_longer(cols = orig:lb_renoise)
# 
# noise_snr_long$name <- factor(noise_snr_long$name,
#                               levels = c("orig", "lb_only", "lb_renoise"))
# 
# noise_snr_long |> ggplot(aes(x = dynamic, y = value)) + geom_line() + 
#   facet_grid(cols = vars(name)) + 
#   geom_hline(yintercept = mean(noise_snr_df$orig), col = "red") + ylab("SNR")
