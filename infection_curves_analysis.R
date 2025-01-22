library(tidyverse)
library(ggrepel)
library(scales)
library(colorBlindness)
library(knitr)
library(survival)
library(tibble)
library(Rcpp)
#install.packages("survminer")
library(survminer)

raw_bite_iv <- read_tsv("C:/Users/Computer/Desktop/Results/bite_back_iv_data_clean.txt", col_types = cols(.default="c")) %>% #import data 
 # select(!inf) %>% # keep inf out
  mutate_at(.vars =c(paste0("field",seq(1,5)), "mean"), #
            .funs = str_replace,
            pattern = ",", 
            replacement = ".") %>% 
  filter(!(mouse == "482" & field3 == "#DIV/0!")) %>% 
  filter(!(mouse == "472" & field3 == "#DIV/0!")) %>% 
  filter(!(mouse == "492" & field3 == "#DIV/0!"))



i_v_mice <- c(481,482, 483, 471, 1, 2, 3, 4)  # best practice: no hard coding, better: additional input table (import)
bb_mice <- c(491, 492, 493, 472, 2891, 2892, 2901, 2902)

bite_iv <- raw_bite_iv %>% # add column for i.v. and bite back 
  mutate(inf = case_when(
    mouse %in% i_v_mice ~ "i_v",
    mouse %in% bb_mice ~ "biteback",
    TRUE ~ NA_character_
  ))%>%
  mutate(mean = case_when( # change Excels NA in Rs NA 
    mean == "#DIV/0!" ~ NA,
    #is.na(as.numeric(mean)) ~ 3, # this was just to check if it works 
    TRUE ~ as.numeric(mean) # convert all mean numbers into integers
  )) %>% 
  mutate(date_1 = as.Date(date, format = "%d.%m.%Y")) %>% # change the date format into the standard to order it 
  arrange(date_1) %>% # order ??
  group_by(construct, inf) %>% # group to sort the dates for each construct (WT and TSRpoint)
  arrange(date_1) %>% # order within one group to get day 1 for each construct
  mutate(d_p_i = dense_rank(date_1)) %>% # add column to say which day post infection 
  ungroup() %>% 
  group_by(construct, inf) %>% 
  mutate(group = cur_group_id()) %>%  # guve id for each group characterized by construct and inf (i.v. and bite back)
  mutate(plotnames = case_when(
    construct == "PbANKA" ~ "wt",
    construct == "TSRswap" ~ "trp1-tsr-swap",
    construct == "TRP1-TSRpoint" ~ "trp1-tsr-point",
    construct == "TRP1-delTSR" ~ paste0("trp1-","\u0394","tsr"),
    TRUE ~ "nothing"
        
  ))


pl_df_iv_bb <- bite_iv %>% # select only parts of the dataframe needed for the plot
  filter(counted == "rate") %>% # only rows containing the parasitemia 
  select(date_1, construct, mouse, mean, inf) %>%
  mutate(mean = as.numeric(mean)) %>% 
  mutate(mean_perc =  percent(mean, accuracy = 0.01)) # add a column for the percentage (nicer numbers in plot)


line_plot_data <- bite_iv %>% # define the data again for the line plot (to overlay)?????
  filter(counted == "rate",
         !is.na(mean)) %>%  # remove all NAs to make the lines continue, data points are missing because stainings were too bad
  mutate(d_p_i = as.factor(d_p_i)) %>% 
  mutate(group = as.factor(group))

# data for mean plots with errorbars 
mp_data <- bite_iv %>% 
  filter(counted == "rate") %>% 
  group_by(group, d_p_i, construct, inf) %>% 
  summarise(mean_p = mean(mean, na.rm = TRUE),
            sd_p = ifelse(n() > 1, sd(mean, na.rm = TRUE), 0)) %>% # calculate standard deviation as new column 
  mutate(group = as.factor(group)) %>% # groups are not numeric (fix x axis to whole days)
  mutate(d_p_i = as.factor(d_p_i)) 
  
blues <- c("#004166", "#00517f" ,"#006299", "#0072b2", "#0082cc", "#0093e5", "#00a3ff")
reds <- c("#893c00", "#a24700", "#bc5300", "#d55e00", "#ef6900", "#ff7609", "#ff8423")

gene_names_italics <- lapply(bite_iv$plotnames, function(x) bquote(italic(.(x)))) %>% 
  unique() %>% 
  as.character()


###############STATISTICS#######################

km_data <- bite_iv %>% 
  filter(counted == "rate") %>% 
  mutate(status = case_when(
    mean == 0 ~ 1,
    mean > 0 ~ 2,
    TRUE ~ NA_real_
  )) #%>%
  # select(group, d_p_i, construct, inf, mean, status) %>% 
  #drop_na()#remove all rows containing any NAs

km_data_bb <- km_data %>% 
  filter(inf == "biteback")

km_data_iv <- km_data %>% 
  filter(inf == "i_v") 


surv_diff_bb <- survdiff(Surv(d_p_i, status) ~ construct, data = km_data_bb)

surv_diff_iv <- survdiff(Surv(d_p_i, status) ~ construct, data = km_data_iv)

# Fit the Kaplan-Meier survival models for each subset
surv_fit_bb <- survfit(Surv(d_p_i, status) ~ construct, data = km_data_bb)
surv_fit_iv <- survfit(Surv(d_p_i, status) ~ construct, data = km_data_iv)

# Generate Kaplan-Meier plot with solid lines and no linetype legend
kaplan_meier_biteback1 <- ggsurvplot(
  surv_fit_bb,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = TRUE, # Add risk table
  risk.table.col = "strata", # Change risk table color by groups
  surv.median.line = "hv", # Specify median survival
  ggtheme = theme_bw(), # Change ggplot2 theme
  palette = c("#0082cc", "#ff8423"), # Specify color palette
  legend.title = "Parasite line",
  legend.labs = c("WT" = "wt", 
                  "TRP1_TSRswap" = "trp1-tsr-swap"), # Use italic labels
  xlab = "Days post infection",
  ylab = "Non-infected proportion (blood stage)"
)

# Modify the ggplot object for solid lines and italic legend labels
kaplan_meier_biteback1$plot <- kaplan_meier_biteback1$plot +
  scale_color_manual(
    values = c("#0082cc", "#ff8423"), # Set colors for groups
    labels = c("wt", "trp1-tsr-swap") # Italic labels
  ) +
  guides(linetype = "none") + # Remove the linetype legend
  theme(legend.text = element_text(face = "italic")) # Italicize legend text if needed


# Generate Kaplan-Meier plot with solid lines and italicized labels
kaplan_meier_iv1 <- ggsurvplot(
  surv_fit_iv,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = TRUE, # Add risk table
  risk.table.col = "strata", # Change risk table color by groups
  surv.median.line = "hv", # Specify median survival
  ggtheme = theme_bw(), # Change ggplot2 theme
  palette = c("#0082cc", "#ff8423"), # Specify color palette
  legend.title = "Parasite line",
  legend.labs = c("wt", 
                  "trp1-tsr-swap"), # Italic labels
  xlab = "Days post infection",
  ylab = "Non-infected proportion (blood stage)"
)

# Modify the ggplot object for solid lines and italic legend labels
kaplan_meier_iv1$plot <- kaplan_meier_iv1$plot +
  scale_color_manual(
    values = c("#0082cc", "#ff8423"), # Set colors for groups
    labels = c("wt", "trp1-tsr-swap") # Italic labels
  ) +
  guides(linetype = "none") + # Remove the linetype legend
  theme(legend.text = element_text(face = "italic")) # Italicize legend text








p_val_transassay <- data.frame(
  comparison = c("TRP1-TSRswap vs WT (biteback)", "TRP1-TSRswap vs WT (iv)"),
  p_value = c(surv_diff_bb$pvalue, surv_diff_iv$pvalue)
)

##################PLOTS###################

#mean plots 

only_mean_plot <- ggplot(mp_data, aes(x = d_p_i, y = mean_p, color = construct))+
                  geom_point(size = 6, alpha = 1, shape = 18) + 
                   geom_line(aes(group = group), size = 1) +
                  geom_errorbar(aes(ymin = mean_p - (1/2)*sd_p, ymax = mean_p + (1/2)*sd_p), width = 0.2, size = 1)+
                  facet_wrap(ncol = 2,
                             facets = "inf",
                             scales = "free_x",
                             labeller =as_labeller(c("biteback" = "natural transmission",
                                                     "i_v" = "1000 sporozoites i.v.")))+
                  labs(x = "days post infection",
                         y = "parasitemia",
                         color = " ")+
                  scale_color_manual(values = c(
                         # "#00517f", 
                          "#0082cc", 
                       # "#d55e00", 
                          "#e53e00" # TSRswap
                    ), 
                 labels = c(expression(italic("wt")), expression(italic("trp1-tsr-swap"))))+
  theme_bw(
    base_size = 16)+
  theme(axis.text = element_text(size = 12))

only_mean_plot
cvdPlot(only_mean_plot)


# test: jeden tag alle werte miteinander vergleichen oder wann positiv logrank 
mean_plus_mice_plot <- ggplot(mp_data, aes(x = d_p_i, y = mean_p, color = construct, shape = "mean")) +
  geom_point(size = 6, alpha = 1, shape = 18) +  # Filled diamond points for the mean
  geom_line(aes(group = group), size = 1) +
  
  facet_grid(~ inf,
             scales = "free_x",
             space = "free",
             labeller = as_labeller(c("biteback" = "Natural transmission",
                                       "i_v" = "1000 sporozoites i.v."))) +
  labs(x = "Days post infection", 
       y = "Parasitemia", 
       color = "Parasite Line", 
       shape = " ") +  # Adding shape legend
  scale_color_manual(values = c(
   # "#00517f", 
    "#0082cc", 
   # "#d55e00", 
    "#e53e00" # TSRswap
  ), 
    labels = c(expression(italic("wt")), 
              expression(italic("trp1-tsr-swap")))) +
  scale_shape_manual(values = c("mean" = 18, "mice" = 16),  # Set different shapes
                     labels = c("Mean", "Single mice")) +  # Label the shapes
  geom_point(data = line_plot_data, aes(x = d_p_i, y = mean, color = construct, shape = "mice"), 
             size = 4, alpha = 0.3) +   # Set shape for the single mice data points
  geom_line(data = line_plot_data, aes(x = d_p_i, y = mean, color  = construct, group = mouse), alpha = 0.3,
            size = 1)+
  theme_bw(
    base_size = 16)+
  theme(axis.text = element_text(size = 12))
  

mean_plus_mice_plot
cvdPlot(mean_plus_mice_plot)

######## TSRdel and TSRpoint data #######
colors_TSR <- c(
  "#0082cc", # WT
  "#711b38",# delTSR
  "#e53e00", # TSRswap
  "#d8c240" # TSRpoint
)



raw_bite_iv_delpoint <- read_tsv("C:/Users/Computer/Desktop/Results/bite_back_iv_data_point_del.txt", col_types = cols(.default="c")) %>%
  select(
    "construct",
    "date", 
    mouse = "mouse Nr",
    d_p_i = "dpi",
    "mean",
    "conditions") %>%
  mutate_at(.vars ="mean", #
            .funs = str_replace,
            pattern = ",",
            replacement = ".")


wt_data <- bite_iv %>% 
  filter(construct == "PbANKA", counted == "rate") %>% 
  select(
    "construct",
    "date",
    "mouse",
    "d_p_i",
    "inf",
    "date_1",
    "group",
    "plotnames",
    "mean"
  ) %>% 
  mutate(replicate = 1, conditions =  case_when(
    mouse %in% bb_mice ~ "bb",
    mouse %in% i_v_mice ~ "iv",
    TRUE ~ NA_character_
  ))

bite_iv_dp1 <- raw_bite_iv_delpoint %>% # add column for i.v. and bite back 
  mutate(mean = as.numeric(mean)) %>%  # convert all mean numbers into integers 
  mutate(date_1 = as.Date(date, format = "%d.%m.%Y")) %>% # change the date format into the standard to order it 
  arrange(date_1) %>% # order ??
  group_by(construct, conditions) %>% # group to sort the dates for each construct (WT and TSRpoint)
  arrange(date_1) %>% # order within one group to get day 1 for each construct
  mutate(d_p_i = dense_rank(date_1)) %>% # add column to say which day post infection 
  ungroup() %>% 
  group_by(construct, conditions) %>% 
  mutate(group = cur_group_id()) %>%  # guve id for each group characterized by construct and inf (i.v. and bite back)
  mutate(plotnames = case_when(
    construct == "PbANKA" ~ "wt",
    construct == "TSRswap" ~ "trp1-tsr-swap",
    construct == "TRP1-TSRpoint" ~ "trp1-tsr-point",
    construct == "TRP1-delTSR" ~ paste0("trp1-","\u0394","tsr"),
    TRUE ~ "nothing"
    
  )) %>% 
  mutate(inf = case_when(
    conditions %in% c("bb1", "bb2", "bb") ~ "biteback",
    conditions %in% c("iv1", "iv2", "iv") ~ "i_v",
    TRUE ~ NA_character_
  ) ) %>% 
  mutate(replicate = case_when(
    conditions %in% c("bb1","iv1", "iv", "bb") ~ 1,
    conditions %in% c("bb2", "iv2") ~ 2,
    TRUE ~ NA_real_
  ))


bite_iv_dp <- bind_rows(bite_iv_dp1, wt_data) %>% 
  unique() %>% 
  group_by(construct, conditions) %>% 
  mutate(group = cur_group_id())

pl_df_iv_bb_dp <- bite_iv_dp %>% # select only parts of the dataframe needed for the plot
  mutate(mean_perc =  percent(mean, accuracy = 0.01)) # add a column for the percentage (nicer numbers in plot)


line_plot_data_dp <- bite_iv_dp %>% # define the data again for the line plot (to overlay)?????
 # filter(counted == "rate",
  #       !is.na(mean)) %>%  # remove all NAs to make the lines continue, data points are missing because stainings were too bad
  mutate(d_p_i = as.factor(d_p_i)) %>% 
  mutate(group = as.factor(group))

# data for mean plots with errorbars 
mp_data_dp_point <- bite_iv_dp %>% 
  filter(construct %in% c("PbANKA","TRP1-TSRpoint")) %>% 
  group_by(group, d_p_i, inf, conditions, construct) %>% 
  summarise(mean_p = mean(mean, na.rm = TRUE),
            sd_p = ifelse(n() > 1, sd(mean, na.rm = TRUE), 0)) %>% # calculate standard deviation as new column 
  mutate(group = as.factor(group)) %>% # groups are not numeric (fix x axis to whole days)
  mutate(d_p_i = as.factor(d_p_i)) 


mp_data_dp_point_simp <- bite_iv_dp %>% 
  filter(construct %in% c("PbANKA","TRP1-TSRpoint")) %>% 
  mutate(group = case_when(
    conditions %in% c("bb1", "bb2") ~ 3,
    conditions %in% c("iv1", "iv2") ~ 4,
    conditions == "bb" ~ 1,
    conditions == "iv" ~ 2,
    TRUE ~ NA_real_
  )) %>% 
  group_by(group, d_p_i, inf, construct) %>% 
  summarise(mean_p = mean(mean, na.rm = TRUE),
            sd_p = ifelse(n() > 1, sd(mean, na.rm = TRUE), 0)) %>% # calculate standard deviation as new column 
  mutate(group = as.factor(group)) %>% # groups are not numeric (fix x axis to whole days)
  mutate(d_p_i = as.factor(d_p_i)) 



mp_data_dp_del <- bite_iv_dp %>% 
  filter(construct  %in% c("PbANKA","TRP1-delTSR")) %>% 
  group_by(group, d_p_i, inf, conditions, construct) %>% 
  summarise(mean_p = mean(mean, na.rm = TRUE),
            sd_p = ifelse(n() > 1, sd(mean, na.rm = TRUE), 0)) %>% # calculate standard deviation as new column 
  mutate(group = as.factor(group)) %>% # groups are not numeric (fix x axis to whole days)
  mutate(d_p_i = as.factor(d_p_i)) 

gene_names_italics_dp <- lapply(bite_iv_dp$plotnames, function(x) bquote(italic(.(x)))) %>% 
  unique() %>% 
  as.character()


###############STATISTICS#######################

km_data_dp <- bite_iv_dp %>% 
  #filter(counted == "rate") %>% 
  mutate(status = case_when(
    mean == 0 ~ 1,
    mean > 0 ~ 2,
    TRUE ~ NA_real_
  )) #%>%
# select(group, d_p_i, construct, inf, mean, status) %>% 
#drop_na()#remove all rows containing any NAs

km_data_bb_dp <- km_data_dp %>% 
  filter(inf == "biteback")

km_data_iv_dp <- km_data_dp %>% 
  filter(inf == "i_v") 

km_data_bb_point <- km_data_bb_dp %>% 
  filter(construct %in% c("PbANKA","TRP1-TSRpoint"))

km_data_iv_point <- km_data_iv_dp %>% 
  filter(construct %in% c("PbANKA","TRP1-TSRpoint")) 

km_data_bb_del <- km_data_bb_dp %>% 
  filter(construct  %in% c("PbANKA","TRP1-delTSR"))

km_data_iv_del <- km_data_iv_dp %>% 
  filter(construct  %in% c("PbANKA","TRP1-delTSR")) 



surv_diff_bb_point <- survdiff(Surv(d_p_i, status) ~ construct, data = km_data_bb_point)

surv_diff_iv_point <- survdiff(Surv(d_p_i, status) ~ construct, data = km_data_iv_point)


surv_diff_bb_del <- survdiff(Surv(d_p_i, status) ~ construct, data = km_data_bb_del)

surv_diff_iv_del <- survdiff(Surv(d_p_i, status) ~ construct, data = km_data_iv_del)

# Fit the Kaplan-Meier survival models for each subset
surv_fit_bb_point <- survfit(Surv(d_p_i, status) ~ construct, data = km_data_bb_point)
surv_fit_iv_point <- survfit(Surv(d_p_i, status) ~ construct, data = km_data_iv_point)

surv_fit_bb_del <- survfit(Surv(d_p_i, status) ~ construct, data = km_data_bb_del)
surv_fit_iv_del <- survfit(Surv(d_p_i, status) ~ construct, data = km_data_iv_del)

# Generate Kaplan-Meier plot with solid lines and no linetype legend
kaplan_meier_biteback1_point <- ggsurvplot(
  surv_fit_bb_point,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = TRUE, # Add risk table
  risk.table.col = "strata", # Change risk table color by groups
  surv.median.line = "hv", # Specify median survival
  ggtheme = theme_bw(), # Change ggplot2 theme
  palette = c("#0082cc", "#d8c240"), # Specify color palette
  legend.title = "Parasite line",
  legend.labs = c("WT" = "wt", 
                  "TRP1_TSRpoint" = "trp1-tsr-point"), # Use italic labels
  xlab = "Days post infection",
  ylab = "Non-infected proportion (blood stage)"
)
# Modify the ggplot object for solid lines and italic legend labels
kaplan_meier_biteback1_point$plot <- kaplan_meier_biteback1_point$plot +
  scale_color_manual(
    values = c("#0082cc", "#d8c240"), # Set colors for groups
    labels = c("wt", "trp1-tsr-point") # Italic labels
  ) +
  guides(linetype = "none") + # Remove the linetype legend
  theme(legend.text = element_text(face = "italic")) # Italicize legend text if needed


kaplan_meier_biteback1_del <- ggsurvplot(
  surv_fit_bb_del,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = TRUE, # Add risk table
  risk.table.col = "strata", # Change risk table color by groups
  surv.median.line = "hv", # Specify median survival
  ggtheme = theme_bw(), # Change ggplot2 theme
  palette = c("#0082cc", "#711b38"), # Specify color palette
  legend.title = "Parasite line",
  legend.labs = c("WT" = "wt", 
                  "TRP1_TSRdel" = paste0("trp1-","\u0394","tsr")), # Use italic labels
  xlab = "Days post infection",
  ylab = "Non-infected proportion (blood stage)"
)
# Modify the ggplot object for solid lines and italic legend labels
kaplan_meier_biteback1_del$plot <- kaplan_meier_biteback1_del$plot +
  scale_color_manual(
    values = c("#0082cc", "#711b38"), # Set colors for groups
    labels = c("wt", paste0("trp1-","\u0394","tsr")) # Italic labels
  ) +
  guides(linetype = "none") + # Remove the linetype legend
  theme(legend.text = element_text(face = "italic")) # Italicize legend text if needed


# Generate Kaplan-Meier plot with solid lines and italicized labels
kaplan_meier_iv1_point <- ggsurvplot(
  surv_fit_iv_point,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = TRUE, # Add risk table
  risk.table.col = "strata", # Change risk table color by groups
  surv.median.line = "hv", # Specify median survival
  ggtheme = theme_bw(), # Change ggplot2 theme
  palette = c("#0082cc", "#d8c240"), # Specify color palette
  legend.title = "Parasite line",
  legend.labs = c("wt", 
                  "trp1-tsr-point"), # Italic labels
  xlab = "Days post infection",
  ylab = "Non-infected proportion (blood stage)"
)

# Modify the ggplot object for solid lines and italic legend labels
kaplan_meier_iv1_point$plot <- kaplan_meier_iv1_point$plot +
  scale_color_manual(
    values = c("#0082cc", "#d8c240"), # Set colors for groups
    labels = c("wt", "trp1-tsr-point") # Italic labels
  ) +
  guides(linetype = "none") + # Remove the linetype legend
  theme(legend.text = element_text(face = "italic")) # Italicize legend text

cvdPlot(kaplan_meier_iv1_del$plot)

 # Generate Kaplan-Meier plot with solid lines and italicized labels
kaplan_meier_iv1_del <- ggsurvplot(
  surv_fit_iv_del,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = TRUE, # Add risk table
  risk.table.col = "strata", # Change risk table color by groups
  surv.median.line = "hv", # Specify median survival
  ggtheme = theme_bw(), # Change ggplot2 theme
  palette = c("#0082cc", "#711b38"), # Specify color palette
  legend.title = "Parasite line",
  legend.labs = c("wt", 
                  paste0("trp1-","\u0394","tsr")), # Italic labels
  xlab = "Days post infection",
  ylab = "Non-infected proportion (blood stage)"
)

# Modify the ggplot object for solid lines and italic legend labels
kaplan_meier_iv1_del$plot <- kaplan_meier_iv1_del$plot +
  scale_color_manual(
    values = c("#0082cc", "#711b38"), # Set colors for groups
    labels = c("wt", paste0("trp1-","\u0394","tsr")) # Italic labels
  ) +
  guides(linetype = "none") + # Remove the linetype legend
  theme(legend.text = element_text(face = "italic")) # Italicize legend text






paste0("trp1-","\u0394","tsr")


p_val_transassay <- data.frame(
  comparison = c("TRP1-TSRswap vs WT (biteback)", "TRP1-TSRswap vs WT (iv)"),
  p_value = c(surv_diff_bb$pvalue, surv_diff_iv$pvalue)
)

##################PLOTS###################


#mean plots 

only_mean_plot_point <- ggplot(mp_data_dp_point_simp, aes(x = d_p_i, y = mean_p, color = construct))+
  geom_point(size = 6, alpha = 1, shape = 18) + 
  geom_line(aes(group = group), size = 1) +
  geom_errorbar(aes(ymin = mean_p - (1/2)*sd_p, ymax = mean_p + (1/2)*sd_p), width = 0.2, size = 1)+
  facet_wrap(ncol = 2,
             facets = "inf",
             scales = "free_x",
             labeller =as_labeller(c("biteback" = "natural transmission",
                                     "i_v" = "1000 sporozoites i.v.")))+
  labs(x = "days post infection",
       y = "parasitemia",
       color = " ")+
  scale_color_manual(values = c(
    # "#00517f", 
    "#0082cc", 
    # "#d55e00", 
    "#d8c240" # TSRpoint
  ), 
  labels = c(expression(italic("wt")), expression(italic("trp1-tsr-point"))))+
  theme_bw(
    base_size = 16)+
  theme(axis.text = element_text(size = 12))

only_mean_plot_point
cvdPlot(only_mean_plot)


mean_plus_mice_plot_point <- ggplot(mp_data_dp_point_simp, aes(x = d_p_i, y = mean_p, color = construct, shape = "mean")) +
  geom_point(size = 6, alpha = 1, shape = 18) +  # Filled diamond points for the mean
  geom_line(aes(group = group), size = 1) +
  
  facet_grid(~ inf,
             scales = "free_x",
             space = "free",
             labeller = as_labeller(c("biteback" = "Natural transmission",
                                      "i_v" = "1000 sporozoites i.v."))) +
  labs(x = "Days post infection", 
       y = "Parasitemia", 
       color = "Parasite Line", 
       shape = " ") +  # Adding shape legend
  scale_color_manual(values = c(
    # "#00517f", 
    "#0082cc", 
    # "#d55e00", 
    "#d8c240" # TSRpoint
  ), 
  labels = c(expression(italic("wt")), 
             expression(italic("trp1-tsr-point")))) +
  scale_shape_manual(values = c("mean" = 18, "mice" = 16),  # Set different shapes
                     labels = c("Mean", "Single mice")) +  # Label the shapes
  geom_point(data = line_plot_data_dp %>% filter(construct %in% c("PbANKA","TRP1-TSRpoint")), aes(x = d_p_i, y = mean, color = construct, shape = "mice"), 
             size = 4, alpha = 0.3) +   # Set shape for the single mice data points
  geom_line(data = line_plot_data_dp %>% filter(construct %in% c("PbANKA","TRP1-TSRpoint")), aes(x = d_p_i, y = mean, color  = construct, group = mouse), alpha = 0.3,
            size = 1)+
  theme_bw(
    base_size = 16)+
  theme(axis.text = element_text(size = 12))




mean_plus_mice_plot_point
cvdPlot(mean_plus_mice_plot_point)


only_mean_plot_del <- ggplot(mp_data_dp_del, aes(x = d_p_i, y = mean_p, color = construct))+
  geom_point(size = 6, alpha = 1, shape = 18) + 
  geom_line(aes(group = group), size = 1) +
  geom_errorbar(aes(ymin = mean_p - (1/2)*sd_p, ymax = mean_p + (1/2)*sd_p), width = 0.2, size = 1)+
  facet_wrap(ncol = 2,
             facets = "inf",
             scales = "free_x",
             labeller =as_labeller(c("biteback" = "natural transmission",
                                     "i_v" = "1000 sporozoites i.v.")))+
  labs(x = "days post infection",
       y = "parasitemia",
       color = " ")+
  scale_color_manual(values = c(
    # "#00517f", 
    "#0082cc", 
    # "#d55e00", 
    "#711b38"), 
  labels = c(expression(italic("wt")), bquote(italic("trp1-" * Delta * "tsr"))))+
  theme_bw(
    base_size = 16)+
  theme(axis.text = element_text(size = 12))

only_mean_plot_del
cvdPlot(only_mean_plot_del)


mean_plus_mice_plot_del <- ggplot(mp_data_dp_del, aes(x = d_p_i, y = mean_p, color = construct, shape = "mean")) +
  geom_point(size = 6, alpha = 1, shape = 18) +  # Filled diamond points for the mean
  geom_line(aes(group = group), size = 1) +
  
  facet_grid(~ inf,
             scales = "free_x",
             space = "free",
             labeller = as_labeller(c("biteback" = "Natural transmission",
                                      "i_v" = "1000 sporozoites i.v."))) +
  labs(x = "Days post infection", 
       y = "Parasitemia", 
       color = " ",
       #color = "Parasite Line", 
       shape = " ") +  # Adding shape legend
  scale_color_manual(values = c(
    # "#00517f", 
    "#0082cc", 
    # "#d55e00", 
    "#711b38"), 
  labels = c(expression(italic("wt")), 
             bquote(italic("trp1-" * Delta * "tsr")))) +
  scale_shape_manual(values = c("mean" = 18, "mice" = 16),  # Set different shapes
                     labels = c("Mean", "Single mice")) +  # Label the shapes
  geom_point(data = line_plot_data_dp %>% filter(construct %in% c("PbANKA","TRP1-delTSR")), aes(x = d_p_i, y = mean, color = construct, shape = "mice"), 
             size = 4, alpha = 0.3) +   # Set shape for the single mice data points
  geom_line(data = line_plot_data_dp %>% filter(construct %in% c("PbANKA","TRP1-delTSR")), aes(x = d_p_i, y = mean, color  = construct, group = mouse), alpha = 0.3,
            size = 1)+
  theme_bw(
    base_size = 16)+
  theme(axis.text = element_text(size = 12))


# Save plots as PDF files 


pdf(file = "D:/Documents/UNI/MASTER/THESIS/Figures Thesis/Exported14012025/transas_del.pdf",
    title = "transas_del",
    height = 4,
    width = 9)
print(mean_plus_mice_plot_del)
dev.off()

pdf(file = "D:/Documents/UNI/MASTER/THESIS/Figures Thesis/Exported14012025/transas_point.pdf",
    title = "transas_point",
    height = 4,
    width = 9)
print(mean_plus_mice_plot_point)
dev.off()

pdf(file = "D:/Documents/UNI/MASTER/THESIS/Figures Thesis/Exported14012025/transas_swap.pdf",
    title = "transas_swap",
    height = 4,
    width = 9)
print(mean_plus_mice_plot)
dev.off()
