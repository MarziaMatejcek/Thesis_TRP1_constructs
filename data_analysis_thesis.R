library(tidyverse)
library(ggrepel)
library(scales)

# import raw data and cleanup of column names to avoid problems with " , /, ." later on
raw_data <- read_tsv("C:/Users/Computer/Desktop/Results/gliding_assays.txt") %>%
  select(
    file_name = "File name", 
    "construct",
    "gliding",
    waving_flipping = "waving/flipping",
    patch_gliding = "patch-gliding",
    non_motile = "non motile",
    "floating",
    "wiggling",
    origin_tissue = "origin tissue",
    "total",
    "date",        
    d_p_i = "d. p. i."     
  ) %>% 
  filter(!is.na(file_name)) #drop the rows that have only 0 and no filename 


plot_data <- raw_data %>% # drop the columns not needed for plots
  select(-file_name,
         -date,
         -d_p_i,
         -total)

## plot flag2 data 

flag2_data <- raw_data %>% # create dataframe for only flag2 and WT data (since flag2 has nothing to do with TSR)
  filter(construct %in% c("TRP1flag2", "PbANKA"),
         !is.na(gliding),
         !str_detect(file_name, "clone741-X")) 


rep_flag2 <- flag2_data %>% # include information that there's 2 replicates for flag2
   group_by(
     construct, origin_tissue, date # group by columns needed for that
   ) %>% 
   tally %>% 
   group_by( 
     construct, origin_tissue
   ) %>%
  mutate(
    rep = seq_len(length(date)) 
    )

mt <- c("gliding", # vector of motility types in the right order
        "patch_gliding", 
        "waving_flipping",
        "wiggling",
        "non_motile",
        "floating" )

flag2_plot <- left_join(flag2_data, rep_flag2, by = c("construct", "origin_tissue", "date")) %>% 
  pivot_longer(cols = mt, # to get the total number of sporos for each experiment long dataframe needed
               names_to = "motility_type",
               values_to = "count") %>% 
  group_by(origin_tissue, construct, motility_type) %>% 
  summarise(summe = sum(count)) %>% 
  mutate(motility_type = factor(motility_type, levels = mt))

# color vector but maybe needs to be changed

mt_colors <- c("#009E73", # gliding 
              "#0072B2", # patchgliding
              "#56B4E9", # waving flipping
              "#F0E442", # wiggling
              "#D55E00", # nonmotile
              "#CC79A7" # floating
              )

flag2_pl <- ggplot(flag2_plot, aes(x = construct, y = summe, fill = motility_type)) +
  facet_wrap(~ origin_tissue, ncol = 2) + 
  geom_col(position = "fill") + 
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = seq(0, 100, 25)) + 
  labs(y = "motility type [%]") + 
  scale_fill_manual(breaks = mt, 
                      labels = c("gliding",
                                 "patch-gliding",
                                 "waving/flipping",
                                 "wiggling",
                                 "non motile",
                                 "floating"),
                      name = "motility type",
                      values = mt_colors) + 
  geom_text_repel(data = flag2_plot %>% filter(summe > 0), 
            aes(label = summe), 
            position = position_fill(vjust = 0.5), 
            direction = "x", 
            color = "black") +
  theme_bw()

# plot all TSR-related constructs
  
tsr_data <- raw_data %>% # create dataframe for only flag2 and WT data (since flag2 has nothing to do with TSR)
  filter(construct %in% c("TRP1_delTSR", "PbANKA", "TRP1_TSRswap", "TRP1_TSRpoint"),
         !is.na(gliding)) %>% 
  pivot_longer(cols = mt,
               names_to = "motility_type",
               values_to = "count") %>% 
  group_by(origin_tissue, construct, motility_type) %>% 
  summarise(summe = sum(count)) %>% 
  mutate(motility_type = factor(motility_type, levels = mt)) %>% 
  mutate(mot_type_simp = case_when(
    motility_type %in% c()
  ))


ggplot(tsr_data, aes(x = construct, y = summe, fill = motility_type)) +
  facet_wrap(~ origin_tissue, ncol = 2, scales = "free") + 
  geom_col(position = "fill") + 
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = seq(0, 100, 25)) + 
  labs(y = "motility type [%]") + 
  scale_fill_manual(breaks = mt, 
                    labels = c("gliding",
                               "patch-gliding",
                               "waving/flipping",
                               "wiggling",
                               "non motile",
                               "floating"),
                    name = "motility type",
                    values = mt_colors) + 
  geom_text(data = tsr_data %>% filter(summe > 0), 
                  aes(label = summe), 
                  position = position_fill(vjust = 0.5), 
                 # direction = "x", 
                  color = "black") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

################# midgut infection rate ##########

raw_data_inf <- read_tsv("C:/Users/Computer/Desktop/Results/infection_rate.txt") %>%
  select(
    construct = "clone", # rename columns to match the ones from the other dataframe
    midgut = "MG",
    oocyst_no = "Oocyst",
    "rep",
    "date") %>% 
  mutate(construct = ifelse(construct == "784", "TRP1-flag2", construct)) %>% # rename the flag2 clone 784 
  filter(!construct == "741-X") # remove data of 741-X since its wt-like 

count_inf <- raw_data_inf %>% 
  group_by(construct) %>% 
  tally()

blue_shades <- c("#1E90FF", # DodgerBlue
                 "#00BFFF", # DeepSkyBlue
                 "#4682B4", # SteelBlue
                 "#5F9EA0", # CadetBlue
                 "#6495ED", # CornflowerBlue
                 "#87CEFA"  # LightSkyBlue
)

ggplot(raw_data_inf, aes(x = construct, y = oocyst_no, color = construct)) + 
  geom_jitter(width = 0.1, height = 0) +
  labs(y = "number of oocysts") +
  scale_color_manual(values = blue_shades) +
  geom_boxplot(alpha = 0) +
  geom_text(data = count_inf, aes(x = construct, label = n)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(raw_data_inf, aes(x = construct, y = oocyst_no, color = construct)) + 
  geom_jitter(width = 0.1, height = 0) +
  labs(y = "number of oocysts") +
  scale_color_manual(values = blue_shades) +
  geom_boxplot(alpha = 0) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

########## sporozoite ratio tissue ######  

raw_data_ratio <- read_tsv("C:/Users/Computer/Desktop/Results/tissue_ratio.txt", col_types = cols(.default="c")) %>%
  select(
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
  mutate(origin_tissue = case_when(origin_tissue == "Salivary" ~ "Salivary Gland",
                                   origin_tissue == "salivary" ~ "Salivary Gland",
                                   origin_tissue == "salivary gland" ~ "Salivary Gland",
                                   origin_tissue == "Salivary gland" ~ "Salivary Gland",
                                   TRUE ~ origin_tissue))

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
  mutate(d_p_i = as.numeric(d_p_i))



  data_ratio_midgut <- data_ratio %>% 
    filter(origin_tissue == "Midgut") %>% 
    select(construct, 
           sporo_mosquito_midgut = "sporo_mosquito",
           d_p_i_midgut = "d_p_i",
           date)
  
  data_ratio_other <- data_ratio %>% 
    filter(!origin_tissue == "Midgut") %>% 
    select(construct, 
           sporo_mosquito_other = "sporo_mosquito",
           d_p_i_other = "d_p_i",
           date,
           origin_tissue)
  
  data_w_ratio <- data_ratio_midgut %>% 
    left_join(data_ratio_other) %>% 
    filter(!is.na(sporo_mosquito_midgut), !is.na(sporo_mosquito_other)) %>% 
    rowwise() %>% 
    mutate(sporo_mosquito_midgut = str_replace(sporo_mosquito_midgut, ",", ".")) %>% 
    mutate(sporo_mosquito_other = str_replace(sporo_mosquito_other, ",", ".")) 
  
  
  data_w_ratio_num <- data_w_ratio %>%  
    mutate_at(.vars = c("sporo_mosquito_other", "sporo_mosquito_midgut", "d_p_i_other"),
              .funs = as.numeric) %>%
    rowwise() %>% 
    mutate(ratio = sporo_mosquito_other/sporo_mosquito_midgut)


hemo <- data_w_ratio_num %>% 
  filter(origin_tissue == "Hemolymph")

sali <- data_w_ratio_num %>% 
  filter(!origin_tissue == "Hemolymph") 

  ggplot(hemo, aes(x = construct, y = ratio, color = d_p_i_other))+ 
    geom_point() +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  
  ggplot(sali, aes(x = construct, y = ratio,, color = d_p_i_other))+ 
    geom_point() +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  df1 <- data_ratio %>% 
    filter(construct == "WT")
  
  colors1 <- c("#D55E00", # Orange
              "#009E73", # Teal
              "#0072B2") # Blue
  
  colors2 <- c("#E69F00", # Orange-yellow
              "#56B4E9", # Sky blue
              "#999999") # Gray
  
  ggplot(data_ratio, aes(x = origin_tissue, y = sporo_mosquito))  +
    geom_point() +
    facet_wrap(~ construct) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) #+ 
    scale_y_continuous(breaks = seq(0, 200000, ), labels = seq(0, 200000, 1000000))
  

col_v <- c("Midgut" =  "#D55E00", # Orange
          "Hemolymph" = "#009E73", # Teal
          " Salivary Gland" =  "#0072B2") # Blue)    
    
    
#sporo_dev <-
ggplot(data_ratio, aes(x = d_p_i, y = sporo_mosquito, color = origin_tissue))+
  geom_point() +
  geom_line(aes(group = origin_tissue)) +
  facet_wrap(facets = "construct")+
  scale_color_manual(values = col_v) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# separate plots for constructs (days vary too much)

# delTSR 

delTSR_meandata <- data_ratio %>% 
  filter(construct == "TRP1-delTSR") %>% 
  mutate(day_fm = case_when(
    d_p_i %in% c(seq(10, 17, 1)) ~ 1,
    d_p_i %in% c(seq(18, 27, 1)) ~ 2,
    TRUE ~ 400
  )) %>% 
  group_by(origin_tissue, day_fm) %>% 
  mutate(mean_value = mean(sporo_mosquito)) 


delTSR_meandata2 <- data_ratio %>% 
  filter(construct == "TRP1-delTSR") %>% 
  mutate(day_fm = case_when(
    d_p_i %in% c(seq(10, 17, 1)) ~ 1,
    d_p_i %in% c(seq(18, 27, 1)) ~ 2,
    TRUE ~ 400
  )) %>% 
  group_by(origin_tissue, day_fm) %>% 
  summarize(
    mean_value = mean(sporo_mosquito, na.rm = TRUE),
    sd_value = sd(sporo_mosquito, na.rm = TRUE), # Standard deviation
    se_value = sd(sporo_mosquito, na.rm = TRUE) / sqrt(n()) # Standard error
  ) %>% 
  ungroup() %>% 
  mutate(day_fm = factor(day_fm))


ggplot(delTSR_meandata2, aes(x = day_fm, y = mean_value, color = origin_tissue)) +
  geom_point() +
  #geom_line() +
  geom_errorbar(aes(ymin = mean_value - se_value, ymax = mean_value + se_value), width = 0.1) + # Error bars (standard error)
  theme_bw() +
  labs(y = "Mean Sporozoites per Mosquito", x = "Day Post Infection") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
  scale_x_discrete(breaks = c(1, 2), labels = c("15 - 17", "19 - 25")) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05)))  # Adjust y-axis

#deltsr_sporo_dev <-
ggplot(data_ratio %>% 
         filter(construct == "TRP1-delTSR",
                cage_no == "2"), aes(x = d_p_i, y = sporo_mosquito, 
                                                 color = origin_tissue, 
                                                # shape = cage_no
                                                 )) + 
  geom_point() +
  scale_color_manual(values = col_v) +
  geom_line(aes(group = origin_tissue)) +
  theme_bw()

ggplot(data_ratio %>% 
         filter(construct == "TRP1-delTSR",
                cage_no == "1"), aes(x = d_p_i, y = sporo_mosquito, 
                                     color = origin_tissue, 
                                     # shape = cage_no
                )) + 
  geom_point() +
  scale_color_manual(values = col_v) +
  geom_line(aes(group = origin_tissue)) +
  theme_bw()


######################### bite back and i. v. data #####################
 raw_bite_iv <- read_tsv("Desktop/Results/bite_back_iv_data_clean.txt", col_types = cols(.default="c")) %>%
      select(!inf) %>% 
      mutate_at(.vars =c(paste0("field",seq(1,5)), "mean"),
                .funs = str_replace,
                pattern = ",", 
                replacement = ".")

i_v_mice <- c(481,482, 483, 471)  # best practice: no hard coding, better: additional input table (import)
bb_mice <- c(491, 492, 493, 472, 2891, 2892, 2901, 2902)

bite_iv <- raw_bite_iv %>% 
  mutate(inf = case_when(
    mouse %in% i_v_mice ~ "i_v",
    mouse %in% bb_mice ~ "biteback",
    TRUE ~ NA_character_
  ))%>%
  mutate(mean = case_when(
    mean == "#DIV/0!" ~ NA,
    #is.na(as.numeric(mean)) ~ 3,
    TRUE ~ as.numeric(mean)
    )) %>% 
  mutate(date_1 = as.Date(date, format = "%d.%m.%Y")) %>% 
  arrange(date_1) %>% 
  group_by(construct) %>% 
  arrange(date_1) %>% 
  mutate(day_p_i = dense_rank(date_1)) %>% 
  ungroup() %>% 
  group_by(construct, inf) %>% 
  mutate(group = cur_group_id())
  


pl_df_iv_bb <- bite_iv %>% 
 filter(counted == "rate") %>% 
  select(date_1, construct, mouse, mean, inf) %>%
  mutate(mean = as.numeric(mean)) %>% 
  mutate(mean_perc =  percent(mean, accuracy = 0.01))
  

#pl_df_iv_bb_wt <- bite_iv %>% 
 # filter(counted == "rate", construct == "PbANKA") %>% 
  #select(date, construct, mouse, mean, inf) %>%
  #mutate(mean = as.numeric(mean)) %>% 
  #mutate(mean_perc =  percent(mean, accuracy = 0.01))


#pl_df_iv_bb_tsrpoint <- bite_iv %>% 
 # filter(counted == "rate", construct == "TSRpoint") %>% 
#  select(date, construct, mouse, mean, inf) %>%
 # mutate(mean = as.numeric(mean)) %>% 
#  mutate(mean_perc =  percent(mean, accuracy = 0.01))




facetlab = c("1" = "bite back", 
             "2" = "i. v.",
             "3" = "bite back")


line_plot_data <- bite_iv %>% 
         filter(counted == "rate",
                !is.na(mean))

ggplot(line_plot_data, aes(x = day_p_i, y = mean, color = construct)) +
  geom_point(
    na.rm = TRUE
    ) +
  geom_line(
    #data=
    aes(group = mouse)
    ) +
  facet_wrap(ncol = 4, facets = "mouse", # labeller =as_labeller(facetlab)
             ) +
  theme_bw() + 
  ylab("parasitemia") +
  xlab("day post infection") + 
  #scale_y_continuous(breaks = seq(0, 5, 10), labels = seq(0, 5, 10)) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1))


ggplot(line_plot_data, aes(x = day_p_i, y = mean, color = mouse)) +
  geom_point(
    na.rm = TRUE
  ) +
  geom_line(
    #data=
    aes(group = mouse)
  ) +
  
  facet_grid(construct~inf)+
  
  # facet_wrap(ncol = 2, facets = "group", # labeller =as_labeller(facetlab)
  # ) +
  theme_bw() + 
  ylab("parasitemia") +
  xlab("day post infection") + 
  #scale_y_continuous(breaks = seq(0, 5, 10), labels = seq(0, 5, 10)) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1))


##############################

blues <- c("#004166", "#00517f" ,"#006299", "#0072b2", "#0082cc", "#0093e5", "#00a3ff")

redsTSRdel <- c("#1e070f", "#330c19", "#471124",  "#5c162e",  "#711b38",   "#852043",  "#9a254d")
redsTSRswap <- c("#893c00", "#a24700", "#bc5300", "#d55e00", "#ef6900", "#ff7609", "#ff8423")
redsTSRpoint <- c("#a66c2e", "#ba7934", "#c9853e", "#cf9252", "#d59f66", "#daab7a", "#e0b88e")

lilacabdel <- c( "#3d1365", "#4a177a", "#571b90", "#641fa5", "#7123ba",  "#7e27d0",  "#8b37d9")
lilacabswap <- c("#723a91",  "#8041a3", "#8f48b5", "#9a5abd", "#a56cc4",  "#b07ecc",  "#bb91d3")

greensflag2 <- c("#0d483c", "#115e4d", "#15735f", "#198971", "#1d9f83", "#21b495", "#25caa6")
greensGFP <- c("#3e691c", "#4a7d21", "#569127", "#62a52c", "#6eb931", "#7aca3a", "#88d04e")

"#a99623"
"#bea827"
"#d3bb2b"
"#d8c240"
"#dcc955"
"#e1d06a"
"#e5d680"








