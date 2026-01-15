library(spant)
library(tidyverse)
library(ggsignif)
library(cowplot)

theme_set(theme_bw())

common_cols <- c("path", "Dynamic", "Glu", "Lac", "NAA", "NAAG", "Cr", "PCr", 
                 "GABA")

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
all_res <- do.call("rbind", csv_list)

# extract run info from file path
split_sep <- str_split(all_res$path, pattern = .Platform$file.sep)
all_res$run <- sapply(split_sep, \(x) x[3])
protocol   <- sapply(split_sep, \(x) x[2])
prot_split <- str_split(protocol, pattern = "_")
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
colnames(table_1) <- c("Fit method", "Analysis type", "Field strength", "Blocks")
all_res$Glu_perc <- NA

glu_perc_change    <- rep(NA, length = length(run_ids))
glu_perc_change_sd <- rep(NA, length = length(run_ids))
for (n in 1:length(run_ids)) {
  id <- run_ids[n]
  av_rest <- mean(all_res[all_res$run_id == id & all_res$state == "REST",]$Glu)
  Glu_perc <- all_res[all_res$run_id == id,]$Glu / av_rest * 100 - 100
  all_res[all_res$run_id == id,]$Glu_perc <- Glu_perc
  task_set <- all_res$run_id == id & all_res$state == "TASK"
  glu_perc_change[n] <- round(mean(all_res[task_set,]$Glu_perc), 2)
  glu_perc_change_sd[n] <- round(sd(all_res[task_set,]$Glu_perc), 2)
}

table_1 <- cbind(table_1, "Mean Glu perc change" = glu_perc_change, 
                 "Mean Glu perc change s.d." = glu_perc_change_sd)

table_1 |> arrange(`Field strength`) |> arrange((`Analysis type`))

mean_dyn_bias   <- table_1 |> filter(`Analysis type` == "dyn-basis") |>
  pull(`Mean Glu perc change`) |> Mod() |> mean() 
mean_static_bias <- table_1 |> filter(`Analysis type` == "static-basis") |>
  pull(`Mean Glu perc change`) |> Mod() |> mean()
mean_preproc_lb <- table_1 |> filter(`Analysis type` == "preproc-lb") |>
  pull(`Mean Glu perc change`) |> Mod() |> mean()

mean_static_bias / mean_dyn_bias    # ~9.0 fold bias reduction
mean_static_bias / mean_preproc_lb  # ~9.5 fold bias reduction

lcm_glu <- all_res |> filter(fit_method == "lcmodel") |>
  ggplot(aes(x = basis_type, y = Glu_perc)) +
  facet_grid(n_blocks ~ field_strength) +
  geom_boxplot(aes(col = state)) + 
  ylab("Glu % change from REST") + xlab("") +
  ggtitle("LCModel") + ylim(c(-6, 6.4)) + theme(legend.position="none")

abfit_glu <- all_res |> filter(fit_method == "abfit-reg") |>
  ggplot(aes(x = basis_type, y = Glu_perc)) +
  facet_grid(n_blocks ~ field_strength) +
  geom_boxplot(aes(col = state)) + 
  ylab("") + xlab("") +
  ggtitle("ABfit") + ylim(c(-6, 6.4))

plot_grid(lcm_glu, abfit_glu, rel_widths = c(1, 1.2))

ggsave("lcm_abfit_glu_bias.png", width = 11, height = 6)

# 5 block analysis

table_1 |> filter(Blocks == "05-blocks") |> arrange(`Fit method`,
                                                    `Analysis type`,
                                                    `Field strength`)

table_2 <- table_1 |> filter(Blocks == "05-blocks") |> arrange(`Fit method`,
                                                    `Analysis type`,
                                                    `Field strength`)

library(tableHTML)
write_tableHTML(tableHTML(table_2), file = 'table.html')
