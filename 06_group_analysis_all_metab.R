library(spant)
library(tidyverse)
library(ggsignif)
library(cowplot)
library(tableHTML)

theme_set(theme_bw())

common_cols <- c("path", "Dynamic", "Glu", "Lac", "NAA", "NAAG", "Cr", "PCr", 
                 "GABA", "GPC", "PCh", "Asp", "Glc", "Asc", "Gln", "GSH",
                 "Ins", "PEth", "sIns", "Tau")

results_glob <- file.path("fit_results_combined", "*",
                          "fit_res_group_unscaled_conc.csv")

all_csvs <- sort(Sys.glob(results_glob))

# read from file and extract subset of useful columns before binding runs
proc_csv <- function(x) {
  res <- read.csv(x)
  res <- res[, common_cols]
  return(res)
}

csv_list <- lapply(all_csvs, proc_csv)
all_res  <- do.call("rbind", csv_list)

# add some combination columns
all_res$tNAA <- all_res$NAA + all_res$NAAG
all_res$tCr  <- all_res$Cr  + all_res$PCr
all_res$tCho <- all_res$GPC + all_res$PCh

# extract run info from file path
split_sep     <- str_split(all_res$path, pattern = .Platform$file.sep)
all_res$run   <- sapply(split_sep, \(x) x[3])
protocol      <- sapply(split_sep, \(x) x[2])
prot_split    <- str_split(protocol, pattern = "_")
protocol_info <- do.call("rbind", prot_split)

colnames(protocol_info) <- c("fit_method", "basis_type", "field_strength",
                             "n_blocks")

all_res <- cbind(all_res, protocol_info)

preproc_lb_3t <- which(all_res$field_strength == "3T-preproc-lb")
all_res[preproc_lb_3t,]$field_strength <- "3T"
all_res[preproc_lb_3t,]$basis_type     <- "preproc-lb"
preproc_lb_7t <- which(all_res$field_strength == "7T-preproc-lb")
all_res[preproc_lb_7t,]$field_strength <- "7T"
all_res[preproc_lb_7t,]$basis_type     <- "preproc-lb"

all_res$Dynamic    <- factor(all_res$Dynamic)

all_res$basis_type <- factor(all_res$basis_type,
                             levels = c("static-basis", "preproc-lb",
                                        "dyn-basis"))

all_res$run_id <- paste(all_res$fit_method, all_res$basis_type,
                        all_res$field_strength, all_res$n_blocks, sep = "_")

levels(all_res$basis_type) <- c("default\nanalysis", "preproc lw\nmatching",
                                "basis lw\nmatching")

block_states    <- c("REST", "TASK", "REST", "TASK", "REST")
all_res_2 <- all_res |> filter(n_blocks == "02-blocks") |> 
  cbind(state = c("REST", "TASK"))

all_res_5 <- all_res |> filter(n_blocks == "05-blocks") |>
  cbind(state = block_states)

all_res_10 <- all_res |> filter(n_blocks == "10-blocks") |>
  cbind(state = rep(block_states, each = 2))

all_res <- rbind(all_res_2, all_res_5, all_res_10)

all_res$n_blocks <- factor(all_res$n_blocks,
                           levels = c("02-blocks", "05-blocks", "10-blocks"))
levels(all_res$n_blocks) <- c("2 blocks", "5 blocks", "10 blocks")

all_res$state <- factor(all_res$state)

all_res$ft_basis <- paste(all_res$field_strength, all_res$basis_type)

# calculate percentage change
run_ids  <- unique(all_res$run_id)
id_split <- str_split(run_ids, "_")
table_1  <- as.data.frame(do.call("rbind", id_split))
colnames(table_1) <- c("Fit method", "Analysis type", "Field strength",
                       "Blocks")

all_res$Glu_perc  <- NA
all_res$Lac_perc  <- NA
all_res$Asp_perc  <- NA
all_res$Glc_perc  <- NA
all_res$GABA_perc <- NA
all_res$tNAA_perc <- NA
all_res$tCr_perc  <- NA
all_res$tCho_perc <- NA
all_res$Asc_perc  <- NA
all_res$Gln_perc  <- NA
all_res$GSH_perc  <- NA

all_res$Ins_perc  <- NA
all_res$PEth_perc <- NA
all_res$sIns_perc <- NA
all_res$Tau_perc  <- NA

glu_perc_change    <- rep(NA, length = length(run_ids))
glu_perc_change_sd <- rep(NA, length = length(run_ids))
glu_perc_change_cd <- rep(NA, length = length(run_ids))

for (n in 1:length(run_ids)) {
  id <- run_ids[n]
  
  # Glu
  av_rest <- mean(all_res[all_res$run_id == id & all_res$state == "REST",]$Glu)
  Glu_perc <- all_res[all_res$run_id == id,]$Glu / av_rest * 100 - 100
  all_res[all_res$run_id == id,]$Glu_perc <- Glu_perc
  
  # Asc
  av_rest <- mean(all_res[all_res$run_id == id & all_res$state == "REST",]$Asc)
  Asc_perc <- all_res[all_res$run_id == id,]$Asc / av_rest * 100 - 100
  all_res[all_res$run_id == id,]$Asc_perc <- Asc_perc
  
  # Lac
  av_rest <- mean(all_res[all_res$run_id == id & all_res$state == "REST",]$Lac)
  Lac_perc <- all_res[all_res$run_id == id,]$Lac / av_rest * 100 - 100
  all_res[all_res$run_id == id,]$Lac_perc <- Lac_perc
  
  # GABA
  av_rest <- mean(all_res[all_res$run_id == id & all_res$state == "REST",]$GABA)
  GABA_perc <- all_res[all_res$run_id == id,]$GABA / av_rest * 100 - 100
  all_res[all_res$run_id == id,]$GABA_perc <- GABA_perc
  
  # Asp
  av_rest <- mean(all_res[all_res$run_id == id & all_res$state == "REST",]$Asp)
  Asp_perc <- all_res[all_res$run_id == id,]$Asp / av_rest * 100 - 100
  all_res[all_res$run_id == id,]$Asp_perc <- Asp_perc
  
  # Glc
  av_rest <- mean(all_res[all_res$run_id == id & all_res$state == "REST",]$Glc)
  Glc_perc <- all_res[all_res$run_id == id,]$Glc / av_rest * 100 - 100
  all_res[all_res$run_id == id,]$Glc_perc <- Glc_perc
  
  # tNAA
  av_rest <- mean(all_res[all_res$run_id == id & all_res$state == "REST",]$tNAA)
  tNAA_perc <- all_res[all_res$run_id == id,]$tNAA / av_rest * 100 - 100
  all_res[all_res$run_id == id,]$tNAA_perc <- tNAA_perc
  
  # tCr
  av_rest <- mean(all_res[all_res$run_id == id & all_res$state == "REST",]$tCr)
  tCr_perc <- all_res[all_res$run_id == id,]$tCr / av_rest * 100 - 100
  all_res[all_res$run_id == id,]$tCr_perc <- tCr_perc
  
  # tCho
  av_rest <- mean(all_res[all_res$run_id == id & all_res$state == "REST",]$tCho)
  tCho_perc <- all_res[all_res$run_id == id,]$tCho / av_rest * 100 - 100
  all_res[all_res$run_id == id,]$tCho_perc <- tCho_perc
  
  # Gln
  av_rest <- mean(all_res[all_res$run_id == id & all_res$state == "REST",]$Gln)
  Gln_perc <- all_res[all_res$run_id == id,]$Gln / av_rest * 100 - 100
  all_res[all_res$run_id == id,]$Gln_perc <- Gln_perc
  
  # GSH
  av_rest <- mean(all_res[all_res$run_id == id & all_res$state == "REST",]$GSH)
  GSH_perc <- all_res[all_res$run_id == id,]$GSH / av_rest * 100 - 100
  all_res[all_res$run_id == id,]$GSH_perc <- GSH_perc
  
  # Ins
  av_rest <- mean(all_res[all_res$run_id == id & all_res$state == "REST",]$Ins)
  Ins_perc <- all_res[all_res$run_id == id,]$Ins / av_rest * 100 - 100
  all_res[all_res$run_id == id,]$Ins_perc <- Ins_perc
  
  # PEth
  av_rest <- mean(all_res[all_res$run_id == id & all_res$state == "REST",]$PEth)
  PEth_perc <- all_res[all_res$run_id == id,]$PEth / av_rest * 100 - 100
  all_res[all_res$run_id == id,]$PEth_perc <- PEth_perc
  
  # sIns
  av_rest <- mean(all_res[all_res$run_id == id & all_res$state == "REST",]$sIns)
  sIns_perc <- all_res[all_res$run_id == id,]$sIns / av_rest * 100 - 100
  all_res[all_res$run_id == id,]$sIns_perc <- sIns_perc
  
  # Tau
  av_rest <- mean(all_res[all_res$run_id == id & all_res$state == "REST",]$Tau)
  Tau_perc <- all_res[all_res$run_id == id,]$Tau / av_rest * 100 - 100
  all_res[all_res$run_id == id,]$Tau_perc <- Tau_perc
  
  task_set <- all_res$run_id == id & all_res$state == "TASK"
  glu_perc_change[n]    <- round(mean(all_res[task_set,]$Glu_perc), 2)
  glu_perc_change_sd[n] <- round(sd(all_res[task_set,]$Glu_perc), 2)
  glu_perc_change_cd[n] <- round(mean(all_res[task_set,]$Glu_perc) / 
                                 sd(all_res[task_set,]$Glu_perc), 2)
}

metab_cols <- c("Glu_perc", "Lac_perc",  "Asp_perc",  "Glc_perc",  "GABA_perc",
                "tCr_perc", "tCho_perc", "tNAA_perc", "Asc_perc",  "Gln_perc",
                "GSH_perc", "Ins_perc",  "PEth_perc", "sIns_perc", "Tau_perc")

sel_cols <- c("field_strength", "basis_type", metab_cols)

# all_res_perc_change_abf <- all_res |>
#                        filter(state == "TASK", fit_method == "abfit-reg",
#                               n_blocks == "5 blocks") |>
#                               pivot_longer(cols = all_of(metab_cols))

all_res_perc_change_abf <- all_res |>
  filter(state == "TASK", fit_method == "abfit-reg",
         n_blocks == "5 blocks") |>
  select(all_of(sel_cols)) |>
  pivot_longer(cols = all_of(metab_cols))
                               

all_res_perc_change_abf$name <- gsub('.{5}$', '', all_res_perc_change_abf$name)

# mean
mean_p <- all_res_perc_change_abf |> ggplot(aes(x = name, y = value)) + 
  stat_summary(fun = mean, aes(col = basis_type), geom = "point") + 
  facet_wrap(~ field_strength, ncol = 2) + geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Mean bias (%)") + xlab("") + theme(legend.title=element_blank())

# sd
sd_p <- all_res_perc_change_abf |> ggplot(aes(x = name, y = value)) + 
  stat_summary(fun = sd, aes(col = basis_type), geom = "point") + 
  facet_wrap(~ field_strength, ncol = 2) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("s.d. bias (%)") + xlab("") + theme(legend.title=element_blank())

cohens_d <- function(x) Mod(mean(x)) / sd(x)

# cd
cd_p <- all_res_perc_change_abf |> ggplot(aes(x = name, y = value)) + 
  stat_summary(fun = cohens_d, aes(col = basis_type), geom = "point") + 
  facet_wrap(~ field_strength, ncol = 2) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("|Cohen's d|") + xlab("") + theme(legend.title=element_blank())

plot_grid(mean_p, sd_p, cd_p, nrow = 3, labels = c("A", "B", "C"), align = "v")
ggsave("abfit_metab_bias_cd.png", width = 8, height = 8)

# plot without Cohen's d
# plot_grid(mean_p, sd_p, nrow = 2, labels = c("A", "B"), align = "v")
# ggsave("abfit_metab_bias.png", width = 8, height = 6)

# all_res_perc_change_abf |>
#   filter(field_strength == "3T", basis_type == "default\nanalysis") |> 
#   ggplot(aes(x = name, y = value)) + 
#   geom_hline(yintercept = 0) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#   ylab("Mean bias (%)") + xlab("") + theme(legend.title=element_blank()) +
#   geom_jitter(alpha = 0.2) +
#   stat_summary(fun = mean, aes(col = basis_type), geom = "point") + 
#   ylim(-10, 10)

all_res_perc_change_lcm <- all_res |> 
                       filter(state == "TASK", fit_method == "lcmodel", 
                              n_blocks == "5 blocks") |>
                       pivot_longer(cols = all_of(metab_cols))

all_res_perc_change_lcm$name <- gsub('.{5}$', '', all_res_perc_change_lcm$name)

# mean
mean_p <- all_res_perc_change_lcm |> ggplot(aes(x = name, y = value)) + 
  stat_summary(fun = mean, aes(col = basis_type), geom = "point") + 
  facet_wrap(~ field_strength, ncol = 2) + geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Mean bias (%)") + xlab("") + theme(legend.title=element_blank())

# sd
sd_p <- all_res_perc_change_lcm |> ggplot(aes(x = name, y = value)) + 
  stat_summary(fun = sd, aes(col = basis_type), geom = "point") + 
  facet_wrap(~ field_strength, ncol = 2) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("s.d. bias (%)") + xlab("") + theme(legend.title=element_blank())

# cd
cd_p <- all_res_perc_change_lcm |> ggplot(aes(x = name, y = value)) + 
  stat_summary(fun = cohens_d, aes(col = basis_type), geom = "point") + 
  facet_wrap(~ field_strength, ncol = 2) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("|Cohen's d|") + xlab("") + theme(legend.title=element_blank())

plot_grid(mean_p, sd_p, cd_p, nrow = 3, labels = c("A", "B", "C"), align = "v")
ggsave("lcmodel_metab_bias_cd.png", width = 8, height = 8)

# plot without Cohen's d
# plot_grid(mean_p, sd_p, nrow = 2, labels = c("A", "B"), align = "v")
# ggsave("lcmodel_metab_bias.png", width = 8, height = 6)
