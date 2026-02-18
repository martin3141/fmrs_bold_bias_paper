library(tidyverse)
library(cowplot)

theme_set(theme_bw())

dyn_res <- read_csv("dynamic_fitting_results.csv", show_col_types = FALSE) |> 
           select(-1)

dyn_res$Model <- factor(dyn_res$Model, levels = c("fixed", "correct", "exact"))

mean_p <- dyn_res |> filter(value_type == "Mean") |> 
          ggplot(aes(x = Metabolite, y = value, col = Model)) + geom_point() +
          facet_wrap(~ Field, ncol = 2) +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
          ylab("Mean bias (%)") + xlab("") +
          theme(legend.title=element_blank()) +
          geom_hline(yintercept = 0)

sd_p <- dyn_res |> filter(value_type == "SD") |> 
        ggplot(aes(x = Metabolite, y = value, col = Model)) + geom_point() +
        facet_wrap(~ Field, ncol = 2) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        ylab("s.d. bias (%)") + xlab("") +
        theme(legend.title=element_blank())

cd_p <- dyn_res |> filter(value_type == "Cohens_d") |> 
        ggplot(aes(x = Metabolite, y = value, col = Model)) + geom_point() +
        facet_wrap(~ Field, ncol = 2) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        ylab("Cohen's d") + xlab("") +
        theme(legend.title=element_blank())

plot_grid(mean_p, sd_p, cd_p, nrow = 3, labels = c("A", "B", "C"), align = "v")
ggsave("dyn_fitting_metab_bias_cd.png", width = 8, height = 8)
