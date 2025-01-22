library(tidyverse)
library(ggrepel)
library(scales)
library(colorBlindness)
library(ggsignif)
library(dunn.test)
#source("")

################# midgut infection rate ##########

# load dataframe
raw_data_inf <- read_tsv("C:/Users/Computer/Desktop/Results/infection_rate.txt") %>%
  select(
    construct = "clone", # rename columns to match the ones from the other dataframe
    midgut = "MG",
    oocyst_no = "Oocyst",
    "rep",
    "date") %>% 
  mutate(construct = ifelse(construct == "784", "TRP1-flag2", construct)) %>% # rename the flag2 clone 784 
  filter(!construct == "741-X") # remove data of 741-X since its wt-like 

# count the n for each construct (mosquitos that were dissected)
count_inf <- raw_data_inf %>% 
  group_by(construct) %>% 
  tally()


# vectors for sorting and filtering with all constructs 
TSRorder <- c("WT", "TRP1-delTSR", "TRP1-TSRswap", "TRP1-TSRpoint")
restorder <- c("WT", "TRP1-TMD-GFP", "TRP1-flag2")
other <- c("TRP1-TMD-GFP", "TRP1-flag2")

# split dataframe into two (topic)
# define df containing TSR data
data_inf_tsr <- raw_data_inf %>% 
  filter(!construct %in% other) %>% 
  mutate(construct = factor(construct, levels = TSRorder))
# define df containing tagged data (flag2 and GFP)
data_inf_tags <- raw_data_inf %>% 
  filter(construct %in% restorder) %>% 
  mutate(construct = factor(construct, levels = restorder))

kruskal_tags <- kruskal.test(data_inf_tags$oocyst_no, data_inf_tags$construct)
tags_p_val <- kruskal_tags$p.value

# adapt colours for each construct 
#############colors#############
colors_TSR <- c(
  "#0082cc", # WT
  "#711b38",# delTSR
  "#e53e00", # TSRswap
  "#d8c240" # TSRpoint
)

"#ef6900" # TSRswap

colors_tags <- c(
  "#0082cc", # WT
  "#6eb931", # GFP
  "#115e4d" # flag2
  )

########### Plots infection rate ##############

# Calculate the fraction of infected mosquitoes
fraction_inf_tsr <- data_inf_tsr %>%
  group_by(construct) %>%
  summarize(
    infected_frac = mean(oocyst_no > 0, na.rm = TRUE), # Fraction where oocyst_no > 0
    n = n() # Number of samples
  )

# Perform Kruskal-Wallis test
kruskal_tsr <- kruskal.test(oocyst_no ~ construct, data = data_inf_tsr)
tsr_p_val <-round(kruskal_tsr$p.value, digits = 5)

# Perform Dunn's test
dunn_tsr <- dunn.test(data_inf_tsr$oocyst_no, data_inf_tsr$construct, 
                          kw = TRUE, label = TRUE, wrap = TRUE, 
                          method = "bh") # Adjust p-values for multiple comparisons

################


significance_data1 <- tibble(
  comparisons = dunn_tsr$comparisons,  # Get the comparisons from Dunn's test
  p_value = dunn_tsr$P.adjusted  # Get the adjusted p-values
) %>%
  mutate(
    stars = case_when(
      p_value <= 0.001 ~ "***",
      p_value <= 0.01 ~ "**",
      p_value <= 0.05 ~ "*",
      TRUE ~ "ns"  # If the p-value is greater than 0.05, no significance stars
    ),
    comparisons = str_split(comparisons, " - ")  # Split the string into a 2-element vector
    ) %>% 
  filter(stars != "ns") %>% 
  mutate(y = 900 + row_number() * 50) # Position for the stars (adjust as needed)


# Create the first plot (with boxplot)
ooc_pl_TSR <- ggplot(data_inf_tsr, aes(x = construct, y = oocyst_no, color = construct)) + 
  geom_jitter(width = 0.1, height = 0, size = 5, alpha = 0.5) + # Avoid overplotting by jittering dots
  geom_boxplot(alpha = 0.3) + # Add boxplot
  labs(y = "Number of oocysts", x = "") + # Y-axis label
  scale_color_manual(values = colors_TSR) + # Set colors
  guides(color = "none") +
  geom_text(data = fraction_inf_tsr, 
            aes(x = construct, y = -50, label = paste0("n = ", n)),
            inherit.aes = FALSE) +
  geom_text(data = fraction_inf_tsr %>% 
              filter(construct == "WT"), 
            aes(x = construct, y = 850, label = paste0( "infected: ", round(infected_frac, 2))),
            inherit.aes = FALSE) +# Add n and infected fraction above jitter points
  geom_text(data = fraction_inf_tsr %>% 
              filter(construct != "WT"), 
            aes(x = construct, y = 850, label =  round(infected_frac, 2)),
            inherit.aes = FALSE)+
  scale_x_discrete(labels = c(
    "WT" = expression(italic("wt")), 
    "TRP1-TSRswap" = expression(italic("trp1-tsr-swap")), 
    "TRP1-delTSR" = expression(italic("trp1-") * Delta * italic("tsr")),
    "TRP1-TSRpoint" = expression(italic("trp1-tsr-point"))
  )) +
 # geom_text( x = 1, y = 700, label = paste0("Kruskal-Wallis p = ",round(tsr_p_val, 2)),
  #           hjust = 0, size = 5, color = "black") + # Add p-value at the bottom left
  # Add pairwise comparison significance (geom_signif)
  geom_signif(
    comparisons = significance_data1$comparisons,
    annotations = significance_data1$stars,  # Use the stars from the significance data
    y_position = significance_data1$y,  # Adjust y_position for each comparison as needed
    tip_length = 0.03  # Adjust length of the "tips" of the lines
  ) + 
  theme_bw(base_size = 17) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 15))


ooc_pl_TSR
cvdPlot(ooc_pl_TSR)
######################

# Calculate the fraction of infected mosquitoes
fraction_inf_tags <- data_inf_tags %>%
  group_by(construct) %>%
  summarize(
    infected_frac = mean(oocyst_no > 0, na.rm = TRUE), # Fraction where oocyst_no > 0
    n = n() # Number of samples
  )

# Perform Kruskal-Wallis test
kruskal_tags <- kruskal.test(oocyst_no ~ construct, data = data_inf_tags)
tags_p_val <- round(kruskal_tags$p.value, digits = 3)

#
################
ooc_pl_tags <- ggplot(data_inf_tags, aes(x = construct, y = oocyst_no, color = construct)) + 
  geom_jitter(width = 0.1, height = 0, size = 5, alpha = 0.5) + # Avoid overplotting by jittering dots
  geom_boxplot(alpha = 0.3) + # Add boxplot
  labs(y = "Number of oocysts", x = "") + # Y-axis label
  scale_color_manual(values = colors_tags) + # Set colors
  guides(color = "none") +
  geom_text(data = fraction_inf_tags, 
            aes(x = construct, y = -50, label = paste0("n = ", n)),
            inherit.aes = FALSE) +
  geom_text(data = fraction_inf_tags %>% 
              filter(construct == "WT"), 
            aes(x = construct, y = 850, label = paste0( "infected: ", round(infected_frac, 2))),
            inherit.aes = FALSE) +# Add n and infected fraction above jitter points
  geom_text(data = fraction_inf_tags%>% 
              filter(construct != "WT"), 
            aes(x = construct, y = 850, label =  round(infected_frac, 2)),
            inherit.aes = FALSE)+
  #annotate("text", x = 0.5, y = 700, label = paste0("Kruskal-Wallis p = ", tags_p_val),
  #         hjust = 0, size = 4.5, color = "black") +
  scale_x_discrete(labels = c(
    "WT" = expression(italic("wt")), 
    "TRP1-TMD-GFP" = expression(italic("trp1-tmd-gfp")), 
    "TRP1-flag2" = expression(italic("trp1-flag-gfp"))))+
  # geom_signif(
  #   comparisons = significance_data1$comparisons,
  #   annotations = significance_data1$stars,  # Use the stars from the significance data
  #   y_position = significance_data1$y,  # Adjust y_position for each comparison as needed
  #   tip_length = 0.03  # Adjust length of the "tips" of the lines
  # ) + 
  theme_bw( base_size = 17) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(size = 15))
ooc_pl_tags
cvdPlot(ooc_pl_tags)

######## sporozoite ratio tissue ######  

raw_data_ratio <- read_tsv("C:/Users/Computer/Desktop/Results/tissue_ratio.txt", col_types = cols(.default="c")) %>%
     select( # load data, rename columns to match everything
    construct = "clone",
    purified = "purified (Accudenz)",
    protease_inhibitor = "protease inhibitor",
    origin_tissue = "Tissue",
    wb_sample = "WB sample", 
    sporo_mosquito = "Sporozoite per Mosquito",
    "date",
    "cage_no",
    d_p_i = "day p.i.",
    "counts") %>% 
  mutate(origin_tissue = case_when(origin_tissue == "Salivary" ~ "Salivary Gland", # change tissue names to be uniform 
                                   origin_tissue == "salivary" ~ "Salivary Gland",
                                   origin_tissue == "salivary gland" ~ "Salivary Gland",
                                   origin_tissue == "Salivary gland" ~ "Salivary Gland",
                                   TRUE ~ origin_tissue)) %>% 
  mutate(
    # Replace commas with dots for decimal representation
    sporo_mosquito = str_replace_all(sporo_mosquito, ",", "."),
    # Remove any remaining periods (for thousand separators)
 #   sporo_mosquito = str_replace_all(sporo_mosquito, "\\.", "")
  ) %>%
  # Convert sporo_mosquito to numeric after replacements
  mutate(sporo_mosquito = as.numeric(sporo_mosquito)) %>% 
  filter(!(construct == "TRP1-TMD-GFP" & date == "27.06.2024")) 


data_ratio <- raw_data_ratio %>%
  mutate(construct = ifelse(construct == "784", "TRP1-flag2", construct)) %>% # rename the flag2 clone 784 
  filter(!construct == "741-X") %>% # remove data of 741-X since its wt-like 
  select(construct, 
         origin_tissue,
         sporo_mosquito,
         d_p_i,
         date, 
         cage_no) %>% 
  rowwise() %>% 
  mutate(sporo_mosquito = str_replace(sporo_mosquito, ",", ".")) %>% 
  filter(!sporo_mosquito %in% c("not counted","N/A")) %>% 
  mutate_at(.vars = c("sporo_mosquito"),
            .funs = as.numeric) %>% 
  mutate(origin_tissue = factor(origin_tissue, levels = c("Midgut", "Hemolymph", "Salivary Gland"))) %>% 
  group_by(origin_tissue) %>% 
  mutate(construct = str_remove(construct, " cage\\s?\\d+$")) %>%   # Removes " cage " and trailing numbers 
  mutate(construct = str_remove(construct, " "))  %>% # removes any spaces 
  mutate(d_p_i = as.factor(d_p_i)) %>% 
  ungroup() %>% 
  filter(!(construct == "TRP1-delTSR" & cage_no == "3" & d_p_i == "15"))




data_ratio_midgut <- data_ratio %>% 
  filter(origin_tissue == "Midgut") %>% 
  select(construct, 
         sporo_mosquito_midgut = "sporo_mosquito",
         d_p_i_midgut = "d_p_i",
        date,
        cage_no)

data_ratio_other <- data_ratio %>% 
  filter(!origin_tissue == "Midgut") %>% 
  select(
        construct, 
         sporo_mosquito_other = "sporo_mosquito",
         d_p_i_other = "d_p_i",
         date,
         origin_tissue)

data_w_ratio <- data_ratio_other %>% 
  left_join(data_ratio_midgut, by=c("construct", "date"))  %>% 
  filter(!is.na(sporo_mosquito_midgut), !is.na(sporo_mosquito_other)) %>% 
  rowwise()%>% 
  mutate(sporo_mosquito_midgut = str_replace(sporo_mosquito_midgut, ",", ".")) %>% 
  mutate(sporo_mosquito_other = str_replace(sporo_mosquito_other, ",", ".")) 


data_w_ratio_num <- data_w_ratio %>%  
  mutate_at(.vars = c("sporo_mosquito_other", "sporo_mosquito_midgut", "d_p_i_other"),
            .funs = as.numeric) %>%
  rowwise() %>% 
  mutate(ratio = sporo_mosquito_other/sporo_mosquito_midgut) %>% 
  
  group_by(construct, origin_tissue, cage_no) %>%
  
  mutate(ratio_mean_cage_no = mean(ratio)) %>% 
  ungroup() %>%
  select(construct, origin_tissue, ratio_mean_cage_no, cage_no) %>%
  mutate(ratio=ratio_mean_cage_no) %>%
  #dplyr::rename(ratio="ratio_mean_cage_no") %>%
  
  unique() %>%
  
  group_by(construct, origin_tissue) %>% 
  mutate(ratio_mean = mean(ratio_mean_cage_no)) %>% 
  ungroup()

ratio_tsr <- data_w_ratio_num %>% 
  filter(construct %in% TSRorder) %>% 
  mutate(construct = factor(construct, levels = TSRorder))

ratio_tags <- data_w_ratio_num %>% 
  filter(construct %in% restorder) %>% 
  mutate(construct = factor(construct, levels = restorder))

means_ratio <- data_w_ratio_num %>% 
  group_by(construct, origin_tissue) %>% 
  select(construct, origin_tissue, ratio_mean) %>% 
  unique()
#############
wilcox_rat_tsr <- data.frame(
  comparison = c("TRP1-delTSR vs WT (Hemolymph)", "TRP1-TSRswap vs WT (Hemolymph)",
                 "TRP1-TSRpoint vs WT (Hemolymph)", "TRP1-delTSR vs WT (Salivary Gland)", 
                 "TRP1-TSRswap vs WT (Salivary Gland)", "TRP1-TSRpoint vs WT (Salivary Gland)"),
  p_value = c(
    wilcox.test(ratio_tsr%>%
                  filter(construct == "TRP1-delTSR", origin_tissue == "Hemolymph") %>%
                  pull(ratio),
                ratio_tsr %>%
                  filter(construct == "WT", origin_tissue == "Hemolymph") %>% 
                  pull(ratio))$p.value,
    wilcox.test(ratio_tsr %>% 
                  filter(construct == "TRP1-TSRswap", origin_tissue == "Hemolymph") %>% 
                  pull(ratio),
                ratio_tsr %>%
                  filter(construct == "WT", origin_tissue == "Hemolymph") %>%
                  pull(ratio))$p.value,
    wilcox.test(ratio_tsr %>% 
                  filter(construct == "TRP1-TSRpoint", origin_tissue == "Hemolymph") %>% 
                  pull(ratio),
                ratio_tsr %>%
                  filter(construct == "WT", origin_tissue == "Hemolymph") %>%
                  pull(ratio))$p.value,
    wilcox.test(ratio_tsr %>% 
                  filter(construct == "TRP1-delTSR", origin_tissue == "Salivary Gland") %>% 
                  pull(ratio),
                ratio_tsr %>%
                  filter(construct == "WT", origin_tissue == "Salivary Gland") %>%
                  pull(ratio))$p.value,
    wilcox.test(ratio_tsr %>% 
                  filter(construct == "TRP1-TSRswap", origin_tissue == "Salivary Gland") %>% 
                  pull(ratio),
                ratio_tsr %>%
                  filter(construct == "WT", origin_tissue == "Salivary Gland") %>%
                  pull(ratio))$p.value,
    wilcox.test(ratio_tsr %>% 
                  filter(construct == "TRP1-TSRpoint", origin_tissue == "Salivary Gland") %>% 
                  pull(ratio),
                ratio_tsr %>%
                  filter(construct == "WT", origin_tissue == "Salivary Gland") %>%
                  pull(ratio))$p.value
  ))


wilcox_rat_tsr$adjusted_p <- p.adjust(wilcox_rat_tsr$p_value, method = "BH")

wilcox_rat_tsr <- wilcox_rat_tsr %>%
  mutate(significance = case_when(
    adjusted_p < 0.001 ~ "***",
    adjusted_p < 0.01 ~ "**",
    adjusted_p < 0.05 ~ "*",
    TRUE ~ "ns"
  ))

# wilcoxon tests on tag construct data 
wilcox_rat_tags <- data.frame(
  comparison = c("TRP1-TMD-GFP vs WT (Hemolymph)", "TRP1-flag2 vs WT (Hemolymph)",
                 "TRP1-TMD-GFP vs WT (Salivary Gland)", "TRP1-flag2 vs WT (Salivary Gland)"),
  p_value = c(
    wilcox.test(ratio_tags%>%
                  filter(construct == "TRP1-TMD-GFP", origin_tissue == "Hemolymph") %>%
                  pull(ratio),
                ratio_tags %>%
                  filter(construct == "WT", origin_tissue == "Hemolymph") %>% 
                  pull(ratio))$p.value,
    wilcox.test(ratio_tags%>%
                  filter(construct == "TRP1-flag2", origin_tissue == "Hemolymph") %>%
                  pull(ratio),
                ratio_tags %>%
                  filter(construct == "WT", origin_tissue == "Hemolymph") %>% 
                  pull(ratio))$p.value,
    wilcox.test(ratio_tags%>%
                  filter(construct == "TRP1-TMD-GFP", origin_tissue == "Salivary Gland") %>%
                  pull(ratio),
                ratio_tags %>%
                  filter(construct == "WT", origin_tissue == "Salivary Gland") %>% 
                  pull(ratio))$p.value,
    wilcox.test(ratio_tags %>% 
                  filter(construct == "TRP1-flag2", origin_tissue == "Salivary Gland") %>% 
                  pull(ratio),
                ratio_tags %>%
                  filter(construct == "WT", origin_tissue == "Salivary Gland") %>%
                  pull(ratio))$p.value
  ))


wilcox_rat_tags$adjusted_p <- p.adjust(wilcox_rat_tags$p_value, method = "BH")

wilcox_rat_tags <- wilcox_rat_tags %>%
  mutate(significance = case_when(
    adjusted_p < 0.001 ~ "***",
    adjusted_p < 0.01 ~ "**",
    adjusted_p < 0.05 ~ "*",
    TRUE ~ "ns"
  ))

############## Plots sporozoite numbers ratio ###########

# plot including the mean values as points TSR constructs
pl_rat_mean_tsr <- ggplot(ratio_tsr, aes(x = construct, y = ratio_mean_cage_no, color = construct)) +
  facet_grid( ~ origin_tissue, scales = "free_y") +
  geom_point(size = 2, alpha = 0.3) +
  scale_color_manual(values = colors_TSR) +
  labs(y = "ratio (HLS/MGS or SGS/MGS)") +
  geom_point(data = means_ratio %>% mutate(cage_no=NA) %>%
               filter(construct %in% TSRorder),
             aes(x = construct, y = ratio_mean, color = construct), 
             size = 4, alpha = 1, 
             shape = 18)+ 
  geom_signif(
    data = ratio_tsr %>% filter(origin_tissue == "Hemolymph"), # Apply to Hemolymph only
    comparisons = list(c("TRP1-delTSR", "WT"), c("TRP1-TSRswap", "WT"), c("TRP1-TSRpoint", "WT")),
    annotations = wilcox_rat_tsr$significance[1:3], # Adjust indices as needed
    y_position = c(0.2, 0.22, 0.24)
  ) +
  # Add geom_signif layer for another facet, e.g., Saliva
  geom_signif(
    data = ratio_tsr %>% filter(origin_tissue == "Salivary Gland"), # Apply to Saliva only
    comparisons = list(c("TRP1-delTSR", "WT"), c("TRP1-TSRswap", "WT"), c("TRP1-TSRpoint", "WT")),
    annotations = wilcox_rat_tsr$significance[4:6], # Adjust indices for Saliva comparisons
    y_position = c(0.2, 0.22, 0.24)
  ) +
  theme_bw( base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.y = element_text(size = 13))
# plot without signif, too small n for testing.
# plot to seperate plots for better scale labels 
pl_rat_mean_tsr

pl_rat_mean_tsr_hl <- ggplot(ratio_tsr %>% 
                               filter(origin_tissue == "Hemolymph"), 
                             aes(x = construct, y = ratio, color = construct)) +
  geom_point(size = 4, alpha = 0.3) +
  scale_color_manual(values = colors_TSR) +
  labs(y = "Ratio HLS/MGS", 
       x = "") +
  guides(color = "none") +
  scale_x_discrete(labels = c(
    "WT" = expression(italic("wt")), 
    "TRP1-delTSR" = expression(italic("trp1-") * Delta * italic("tsr")), 
    "TRP1-TSRswap" = expression(italic("trp1-tsr-swap")),
    "TRP1-TSRpoint" = expression(italic("trp1-tsr-point"))
  ))+
  geom_point(data = means_ratio %>% 
               filter(construct %in% TSRorder,
                      origin_tissue == "Hemolymph"),
             aes(x = construct, y = ratio_mean, color = construct), 
             size = 6, alpha = 1, shape = 18)+ 
  theme_bw( base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.y = element_text(size = 15))
pl_rat_mean_tsr_hl


pl_rat_mean_tsr_sg <- ggplot(ratio_tsr %>% 
                               filter(origin_tissue == "Salivary Gland"), 
                             aes(x = construct, y = ratio, color = construct)) +
  geom_point(size = 4, alpha = 0.3) +
  scale_color_manual(values = colors_TSR) +
  labs(y = "Ratio SGS/MGS",
       x = "") +
  guides(color = "none") +
  scale_x_discrete(labels = c(
    "WT" = expression(italic("wt")), 
    "TRP1-delTSR" = expression(italic("trp1-") * Delta * italic("tsr")), 
    "TRP1-TSRswap" = expression(italic("trp1-tsr-swap")),
    "TRP1-TSRpoint" = expression(italic("trp1-tsr-point"))
  ))+
  geom_point(data = means_ratio %>% 
               filter(construct %in% TSRorder,
                      origin_tissue == "Salivary Gland"),
             aes(x = construct, y = ratio_mean, color = construct), 
             size = 6, alpha = 1, shape = 18)+ 
  theme_bw( base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.y = element_text(size = 15))
pl_rat_mean_tsr_sg

# plot with only datapoints Tag constructs 
pl_rat_tags <-  ggplot(ratio_tags, aes(x = construct, y = ratio, color = construct)) +
  facet_grid( ~ origin_tissue) +
  labs(y = "ratio (HLS/MGS or SGS/MGS)") +
    geom_jitter(size = 2, height = 0, width = 0.1) +
  scale_color_manual(values = colors_tags) +
  theme_bw( base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.y = element_text(size = 13))

pl_rat_tags
# plot without signif, too small n for testing.
# plot to seperate plots for better scale labels 

pl_rat_mean_tags_hl <- ggplot(ratio_tags %>% 
                               filter(origin_tissue == "Hemolymph"), 
                             aes(x = construct, y = ratio, color = construct)) +
  geom_point(size = 4, alpha = 0.3) +
  scale_color_manual(values = colors_tags) +
  labs(y = "Ratio HLS/MGS",
       x = "") +
  guides(color = "none") +
  scale_x_discrete(labels = c(
    "WT" = expression(italic("wt")), 
    "TRP1-TMD-GFP" = expression(italic("trp1-tmd-gfp")), 
    "TRP1-flag2" = expression(italic("trp1-flag-gfp"))
  ))+
  geom_point(data = means_ratio %>% 
               filter(construct %in% restorder,
                      origin_tissue == "Hemolymph"),
             aes(x = construct, y = ratio_mean, color = construct), 
             size = 6, alpha = 1, shape = 18)+ 
  theme_bw( base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.y = element_text(size = 15))
pl_rat_mean_tags_hl


pl_rat_mean_tags_sg <- ggplot(ratio_tags %>% 
                               filter(origin_tissue == "Salivary Gland"), 
                             aes(x = construct, y = ratio, color = construct)) +
  geom_point(size = 4, alpha = 0.3) +
  scale_color_manual(values = colors_tags) +
  labs(y = "Ratio SGS/MGS",
       x = "") +
  guides(color = "none") +
  scale_x_discrete(labels = c(
    "WT" = expression(italic("wt")), 
    "TRP1-TMD-GFP" = expression(italic("trp1-tmd-gfp")), 
    "TRP1-flag2" = expression(italic("trp1-flag-gfp"))
  ))+
  geom_point(data = means_ratio %>% 
               filter(construct %in% restorder,
                      origin_tissue == "Salivary Gland"),
             aes(x = construct, y = ratio_mean, color = construct), 
             size = 6, alpha = 1, shape = 18)+ 
  theme_bw( base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.y = element_text(size = 15))
pl_rat_mean_tags_sg


# plot with mean values and significance testing included (Wilcoxon)
pl_rat_mean_tags_sig <- ggplot(ratio_tags, aes(x = construct, y = ratio, color = construct)) +
  facet_grid(~ origin_tissue) +
  labs(y = "ratio (HLS/MGS or SGS/MGS)") +
    geom_point(size = 2, alpha = 0.3) +
  scale_color_manual(values = colors_tags) +
  geom_point(data = means_ratio %>% 
               filter(construct %in% restorder),
             aes(x = construct, y = ratio_mean, color = construct), 
             size = 4, alpha = 1, shape = 18) + 
  # Add geom_signif layer for Hemolymph facet only
  geom_signif(
    data = ratio_tags %>% filter(origin_tissue == "Hemolymph"), # Apply to Hemolymph only
    comparisons = list(c("TRP1-TMD-GFP", "WT"), c("TRP1-flag2", "WT")),
    annotations = wilcox_rat_tags$significance[1:2], # Adjust indices as needed
    y_position = c(0.25, 0.3)
  ) +
  # Add geom_signif layer for another facet, e.g., Saliva
  geom_signif(
    data = ratio_tags %>% filter(origin_tissue == "Salivary Gland"), # Apply to Saliva only
    comparisons = list(c("TRP1-TMD-GFP", "WT"), c("TRP1-flag2", "WT")),
    annotations = wilcox_rat_tags$significance[3:4], # Adjust indices for Saliva comparisons
    y_position = c(0.25, 0.3)
  ) +
  theme_bw( base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.y = element_text(size = 13))

####### data absoulte sporo no ####

 
data_ratio_abs <- data_ratio %>% 
  group_by(construct, origin_tissue, cage_no) %>% 
  mutate(sporo_mosquito_mean = mean(sporo_mosquito))



means_abs <- data_ratio_abs %>% 
  filter(!is.na(sporo_mosquito)) %>% 
  group_by(construct, origin_tissue) %>% 
  summarise(mean_val = mean(sporo_mosquito_mean)) %>% 
  ungroup() 

pl_tsr_abs <- ggplot(data_ratio_abs %>% 
                       filter(construct %in% TSRorder) %>% 
                       mutate(construct = factor(construct, levels = TSRorder)), 
                              aes(x = construct, y = sporo_mosquito_mean, color = construct)) +
  facet_grid( ~ origin_tissue, scales = "free_y",
              labeller = as_labeller(c("Midgut" = "Midgut",
                                       "Hemolymph" = "Haemolymph",
                                       "Salivary Gland" = "Salivary gland"))) +
  geom_point(size = 4, alpha = 0.3) +
  scale_color_manual(values = colors_TSR) +
  labs(y = "Sporozoite count",
       x = "") +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x)))  +  # Log scale for y-axis
  guides(color = "none") +
  scale_x_discrete(labels = c(
    "WT" = expression(italic("wt")), 
    "TRP1-delTSR" = expression(italic("trp1-") * Delta * italic("tsr")), 
    "TRP1-TSRswap" = expression(italic("trp1-tsr-swap")),
    "TRP1-TSRpoint" = expression(italic("trp1-tsr-point"))
  ))+
    geom_point(data = means_abs %>%
               filter(construct %in% TSRorder) %>% 
                 mutate(construct = factor(construct, levels = TSRorder)),
             aes(x = construct, y = mean_val, color = construct),
             size = 6, alpha = 1, shape = 18)+
  theme_bw( base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.y = element_text(size = 15))
pl_tsr_abs


pl_tags_abs <- ggplot(data_ratio_abs %>% 
                       filter(construct %in% restorder) %>% 
                       mutate(construct = factor(construct, levels = restorder)), 
                     aes(x = construct, y = sporo_mosquito_mean, color = construct)) +
  facet_grid( ~ origin_tissue, scales = "free_y",
              labeller = as_labeller(c("Midgut" = "Midgut",
                                       "Hemolymph" = "Haemolymph",
                                       "Salivary Gland" = "Salivary gland")))+
  geom_point(size = 4, alpha = 0.3) +
  scale_color_manual(values = colors_tags) +
  labs(y = "Sporozoite count",
       x = "") +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x)))  +  # Log scale for y-axis
  guides(color = "none") +
  scale_x_discrete(labels = c(
    "WT" = expression(italic("wt")), 
    "TRP1-TMD-GFP" = expression(italic("trp1-tmd-gfp")), 
    "TRP1-flag2" = expression(italic("trp1-flag-gfp"))
  ))+
  geom_point(data = means_abs %>%
               filter(construct %in% restorder) %>% 
               mutate(construct = factor(construct, levels = restorder)),
             aes(x = construct, y = mean_val, color = construct),
             size = 6, alpha = 1, shape = 18)+
  theme_bw( base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.y = element_text(size = 15))
pl_tags_abs

############# sporozoite numbers over time ######

early<- c(10:16)
late <- c(17:26)

# data_ratio_newdates <- data_ratio %>% 
#   mutate(newdate = case_when(
#     d_p_i %in% early ~ 1,
#     d_p_i %in% late ~ 2,
#     TRUE ~ NA
#   )) %>% 
#   mutate(datename = case_when(
#     newdate == 1 ~ "15 - 17",
#     newdate == 2 ~ "19 - 25",
#     TRUE ~ "no"
#   ))
# 
# 
# ggplot(data_ratio_newdates %>% 
#          filter(construct %in% c("TRP1-delTSR", "WT")), aes(x = d_p_i, y = sporo_mosquito, color = construct))+
#   geom_point()+
#   labs(y = "sporozoite numbers per mosquito", x = "days post infectio") +
#     facet_grid(~ origin_tissue, scales = "free_x")+
#   theme_bw( base_size = 14)+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
# 
# 
# 



# plot with only datapoints TSR constructs 
pl_rat_tsr <-  ggplot(ratio_tsr, aes(x = construct, y = ratio, color = construct)) +
  facet_grid( ~ origin_tissue) +
  geom_jitter(size = 2, height = 0, width = 0.1) +
  labs(y = "ratio (HLS/MGS or SGS/MGS)") +
  scale_color_manual(values = colors_TSR) +
  theme_bw( base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))






blue_pastel_palette2 <- c("#ff6f69", "#ff8b94", "#ffaaa5", "#ffd3b6", "#c1dced", "#a8d6e6", "#68a8c4", "#619fbf", "#5a96ba")


# pl_tsr_dpi <- ggplot(data_ratio %>% 
#                    filter(construct %in% TSRorder), aes(x = construct, y = sporo_mosquito, color = d_p_i)) +
#   facet_grid( ~ origin_tissue) +
#   geom_jitter(size = 2, height = 0, width = 0.1) +
#   scale_color_manual(values = blue_pastel_palette2) +
#   theme_bw( base_size = 14) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
# 
# pl_tag_dpi <-  ggplot(data_ratio %>% 
#                         filter(construct %in% restorder), aes(x = construct, y = sporo_mosquito, color = d_p_i)) +
#   facet_grid( ~ origin_tissue) +
#   geom_jitter(size = 2, height = 0, width = 0.1) +
#   scale_color_manual(values = blue_pastel_palette2) +
#   theme_bw( base_size = 14) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
# 
# 
# pl_tsr_dpi_time <- ggplot(data_ratio %>% 
#                             filter(construct %in% TSRorder), aes(x))
# 

#################################################################################


col_v <- c("Midgut" =  "#D55E00", # Orange
           "Hemolymph" = "#009E73", # Teal
           "Salivary Gland" =  "#0072B2") # Blue)    

# Calculate mean values by group
df_dpi_means <- data_ratio %>% 
  group_by(construct, d_p_i, origin_tissue) %>% 
  summarise(mean_sporo_nr = mean(sporo_mosquito, na.rm = TRUE)) %>% 
  ungroup()

# Join mean values back to original data
df_dpi <- data_ratio %>% 
  group_by(construct, d_p_i, origin_tissue) %>% 
  mutate( sd_value = sd(sporo_mosquito, na.rm = TRUE), # Standard deviation
          se_value = sd(sporo_mosquito, na.rm = TRUE) / sqrt(n())) %>%  # Standard error
  ungroup() %>% 
  left_join(df_dpi_means, by = c("construct", "d_p_i", "origin_tissue"))

shape_val <- c(
  "Midgut" =  15, # square
  "Hemolymph" = 17, # triangle
  "Salivary Gland" =  19 # circle
  )
   "#5F1C76"
color_val <- c(
  "Midgut" = "#1A4369",
  "Hemolymph" = "#a739bf",
  "Salivary Gland" = "#6FA8DC"
)

sporo_dev <- ggplot(df_dpi %>% 
         filter(construct %in% TSRorder), aes(x = d_p_i, y = mean_sporo_nr, shape = origin_tissue, color = origin_tissue)) +
  geom_point(size = 2) +
  labs(y = "sporozoite numbers per mosquito", x = "days post infection") +
    geom_line(aes(group = interaction(origin_tissue, construct))) +  # Group by both origin_tissue and construct
  facet_wrap(facets = "construct", scales = "free") +
  geom_errorbar(aes(ymin = mean_sporo_nr - se_value, ymax = mean_sporo_nr + se_value), width = 0.1) + 
  scale_y_log10() +  # Log scale for y-axis
  scale_shape_manual(values = shape_val) +
  scale_color_manual(values = color_val) +
  theme_bw( base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

cvdPlot(sporo_dev)

meandata_time <- data_ratio %>% 
  mutate(day_fm = case_when(
    origin_tissue == "Midgut" & construct != "TRP1-TMD-GFP" & d_p_i %in% 10:16  ~ 1,
    origin_tissue == "Midgut" & construct != "TRP1-TMD-GFP" & d_p_i %in% 17:27 ~ 2,
    origin_tissue == "Midgut" & construct == "TRP1-TMD-GFP" & d_p_i %in% 10:17  ~ 1,
    origin_tissue == "Midgut" & construct == "TRP1-TMD-GFP" & d_p_i %in% 18:27 ~ 2,
    origin_tissue == "Hemolymph" ~ 1,
    origin_tissue == "Salivary Gland" ~ 2,
    TRUE ~ NA_real_ )) %>%  # Assign NA if none of the conditions are met
  group_by(construct, origin_tissue, day_fm) %>% 
  summarize(
    mean_value = mean(sporo_mosquito, na.rm = TRUE),
    sd_value = sd(sporo_mosquito, na.rm = TRUE), # Standard deviation
    se_value = sd(sporo_mosquito, na.rm = TRUE) / sqrt(n()) # Standard error
  ) %>% 
  ungroup() %>% 
  mutate(day_fm = factor(day_fm)) %>% 
  mutate(data_type = "mean")


data_ratio1 <- data_ratio  %>% 
  mutate(day_fm = case_when(
  origin_tissue == "Midgut" & d_p_i %in% 10:16 ~ 1,
  origin_tissue == "Midgut" & d_p_i %in% 17:27 ~ 2,
  origin_tissue == "Hemolymph" ~ 1,
  origin_tissue == "Salivary Gland" ~ 2,
  TRUE ~ NA_real_ )) %>% 
  mutate(data_type = "singlevalue")

tsr_time <- ggplot(meandata_time %>% 
         filter(construct %in% TSRorder), aes(x = day_fm, y = mean_value,#
                                              color = origin_tissue, shape = data_type)) +
  facet_grid(~ construct) + 
  geom_point(size = 4, alpha = 1) +  
  # geom_errorbar(aes(ymin = mean_value - se_value, ymax = mean_value + se_value), width = 0.1) + # Error bars (standard error)
  geom_point(data = data_ratio1 %>% 
               filter(construct %in% TSRorder), 
             aes(x = day_fm, y = sporo_mosquito,
                 color = origin_tissue, shape = data_type), 
             size = 2, alpha = 0.3) + 
  scale_shape_manual(values = c("mean" = 18, "singlevalue" = 16),
                     labels = c("Mean", "Replicate")) +
  scale_color_manual(values = color_val) +
  theme_bw( base_size = 14) +
  scale_y_log10() +  # Log scale for y-axis
  labs(y = "Mean Sporozoites per Mosquito", x = "Day Post Infection", 
       shape = "Data Type", color = "Origin Tissue") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_x_discrete(breaks = c(1, 2), labels = c("15 - 17", "19 - 25"))

tags_time <- ggplot(meandata_time %>% 
         filter(construct %in% restorder), aes(x = day_fm, y = mean_value,#
                                              color = origin_tissue, shape = data_type)) +
  facet_grid(~ construct) + 
  geom_point(size = 4, alpha = 1) +  
  # geom_errorbar(aes(ymin = mean_value - se_value, ymax = mean_value + se_value), width = 0.1) + # Error bars (standard error)
  geom_point(data = data_ratio1 %>% 
               filter(construct %in% restorder), 
             aes(x = day_fm, y = sporo_mosquito,
                 color = origin_tissue, shape = data_type), 
             size = 2, alpha = 0.3) + 
  scale_shape_manual(values = c("mean" = 18, "singlevalue" = 16),
                     labels = c("Mean", "Replicate")) +
  scale_color_manual(values = color_val) +
  theme_bw( base_size = 14) +
  scale_y_log10() +  # Log scale for y-axis
  labs(y = "Mean Sporozoites per Mosquito", x = "Day Post Infection", 
       shape = "Data Type", color = "Origin Tissue") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_x_discrete(breaks = c(1, 2), labels = c("15 - 17", "19 - 25"))


# Export figures as PDF

pdf(file = "D:/Documents/UNI/MASTER/THESIS/Figures Thesis/Exported14012025/abs_no_tags.pdf",
    title = "abs_no_tags",
    height = 6,
    width = 8)
print(pl_tags_abs)
dev.off()

pdf(file = "D:/Documents/UNI/MASTER/THESIS/Figures Thesis/Exported14012025/abs_no_tsr.pdf",
    title = "abs_no_tsr",
    height = 6,
    width = 8)
print(pl_tsr_abs)
dev.off()

pdf(file = "D:/Documents/UNI/MASTER/THESIS/Figures Thesis/Exported14012025/ooc_tag.pdf",
    title = "ooc_tag",
    height = 6,
    width = 8)
print(ooc_pl_tags)
dev.off()

pdf(file = "D:/Documents/UNI/MASTER/THESIS/Figures Thesis/Exported14012025/ooc_TSR_stat.pdf",
    title = "ooc_TSR_stat",
    height = 6,
    width = 8)
print(ooc_pl_TSR)
dev.off()

pdf(file = "D:/Documents/UNI/MASTER/THESIS/Figures Thesis/Exported14012025/ratio_hl_tsr.pdf",
    title = "ratio_hl_tags",
    height = 5,
    width = 6)
print(pl_rat_mean_tags_hl)
dev.off()

pdf(file = "D:/Documents/UNI/MASTER/THESIS/Figures Thesis/Exported14012025/ratio_hl_tsr.pdf",
    title = "ratio_hl_tsr",
    height = 5,
    width = 6)
print(pl_rat_mean_tsr_hl)
dev.off()

pdf(file = "D:/Documents/UNI/MASTER/THESIS/Figures Thesis/Exported14012025/ratio_sg_tags.pdf",
    title = "ratio_sg_tags",
    height = 5,
    width = 6)
print(pl_rat_mean_tags_sg)
dev.off()

pdf(file = "D:/Documents/UNI/MASTER/THESIS/Figures Thesis/Exported14012025/ratio_sg_tsr.pdf",
    title = "ratio_sg_tsr",
    height = 5,
    width = 6)
print(pl_rat_mean_tsr_sg)
dev.off()


pdf(file = "D:/Documents/UNI/MASTER/THESIS/Figures Thesis/Exported14012025/sporo_dev.pdf",
    title = "sporo_dev",
    height = 6,
    width = 8)
print(sporo_dev)
dev.off()

pl_rat_mean_tsr
pl_rat_mean_tags_hl
pl_rat_tsr
pl_rat_mean_tags_sg
pl_rat_mean_tags_sig
pl_tags_abs
ooc_pl_tags
ooc_pl_TSR
