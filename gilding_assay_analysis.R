library(tidyverse)
library(ggrepel)
library(scales)
#install.packages("colorBlindness")
library(colorBlindness)
library(ggsignif)

# import raw data and cleanup of column names to avoid problems with " , /, ." later on
raw_data <- read_tsv("C:/Users/Computer/Desktop/Results/gliding_assays.txt") %>%
  select(
    file_name = "File name", # fix category names 
    "construct",
    "gliding",
    waving_flipping = "waving/flipping", 
    patch_gliding = "patch-gliding",
    non_motile = "non motile",
    "floating",
    flexing = "wiggling",
    origin_tissue = "origin tissue",
    "total",
    "date",        
    d_p_i = "d. p. i."     
  ) %>% 
  filter(!is.na(gliding), !is.na(file_name))  %>% #drop the rows that have only 0 and no filename and NAs 
  mutate(origin_tissue = case_when(
           origin_tissue == "hemolymph" ~ "Hemolymph", # fix tissue names 
           origin_tissue == "salivary gland" ~ "Salivary gland",
           origin_tissue == "midgut" ~ "Midgut",
           TRUE ~ "NA"
         ))

# list of analysed construct in order for the plots
order_tsr <- c("PbANKA", "TRP1_delTSR", "TRP1_TSRswap", "TRP1_TSRpoint" ) 

mt <- c("gliding", # vector of motility types in the right order
        "patch_gliding", 
        "waving_flipping",
        "flexing",
        "non_motile",
        "floating" )

non_productive <- c("patch_gliding", # vector including all non-productive motility types 
                    "waving_flipping",
                    "flexing",
                    "non_motile")

simp <- c("gliding", # vector with simplified motility types
          "non_productive",
          "floating")

longdf <- raw_data %>%  # create a new long df 
  pivot_longer(cols = mt, # columns to transform into long data frame (mt is vector with all motility types)
              names_to = "motility_type",
              values_to = "count1") %>%  # where the new values go 
  mutate(motility_type_simp = case_when( # add column for simplified motility types
    motility_type %in% non_productive ~ "non_productive",
    TRUE ~ motility_type
    )) %>% 
  filter(complete.cases(.)) # remove all rows containing any NAs 
# use long dataframe for further plots etc. 


total_n <- longdf %>%  # to display total n per tissue type and construct sum up all entries for all motility types 
  group_by(origin_tissue, construct) %>% 
  summarize(n = sum(count1), .groups = "drop") # why groups drop?

##### plot flag2 data #####

flag2_data <- longdf %>% # create dataframe for only flag2 and WT data (since flag2 has nothing to do with TSR)
  filter(construct %in% c("TRP1flag2", "PbANKA"),
          !str_detect(file_name, "clone741-X")) # remove clone thats wildtype-like


rep_flag2 <- flag2_data %>% # include information that there's 2 replicates for flag2
  group_by(
    construct, origin_tissue, date # group by columns needed for that
  ) %>% 
  tally %>% # ?
  group_by( 
    construct, origin_tissue
  ) %>%
  mutate(
    rep = seq_len(length(date)) # ? 
  )


flag2_plot_data <- left_join(flag2_data, rep_flag2, by = c("construct", "origin_tissue", "date")) %>% # data for the plot
  group_by(origin_tissue, construct, motility_type) %>% # to get total number of counts for each type (sum for all videos of one kind)
  summarise(summe = sum(count1)) %>% # sum up the counts within one group
  mutate(motility_type = factor(motility_type, levels = mt)) # bring the types in order for the plot

# color vectors for all motility types 

blue_pastel_palette <- c("#68a8c4", # Gliding 
                         "#a8d6e6", # patch-gliding
                         "#c1dced", # waving/flipping
                         "#ffd3b6", # flexing
                         "#ffaaa5", # non motile
                         "#ff8b94") # floating

bluesimp <- c("#68a8c4", # Gliding
              "#ffd3b6", # Non productive
              "#ff8b94") # floating


flag2_plot_data_simp <- left_join(flag2_data, rep_flag2, by = c("construct", "origin_tissue", "date")) %>% 
  group_by(origin_tissue, construct, motility_type_simp) %>% # to get total number of counts for each type (sum for all videos of one kind)
  summarise(summe = sum(count1)) %>% # sum up the counts within one group
  mutate(motility_type_simp = factor(motility_type_simp, levels = simp)) # bring the types in order for the plot

# total n of examined sporozoites in wt and flag2 groups
total_n_flag2 <- total_n %>% 
  filter(construct %in% c("TRP1flag2", "PbANKA"))

# plot the simplified categories 
simp_plot_flag2 <- 
  ggplot(flag2_plot_data_simp, aes(x = construct, y = summe, fill = motility_type_simp)) +
    facet_wrap(~ origin_tissue, ncol = 2, 
               labeller = as_labeller(c("Hemolymph" = "Haemolymph",
                                        "Salivary gland" = "Salivary gland"))) +
  geom_col(position = "fill", width = 0.5, color = "black", size = 0.3) +
    scale_y_continuous(breaks = seq(0, 1, 0.25), labels = seq(0, 100, 25)) + 
    labs(y = "Motility pattern [%]", x = "") + 
    scale_fill_manual(breaks = simp, 
                      labels = c("Gliding",
                                 "Non productive",
                                 "Floating"),
                      name = "Motility pattern",
                      values = bluesimp) + 
  scale_x_discrete(labels = c("PbANKA" = expression(italic("wt")), 
                              "TRP1flag2" = expression(italic("trp1-flag-gfp")))) +
   # geom_text(data = flag2_plot_data_simp %>% filter(summe > 0), # absolute count for each category, removed for results, kept in supplementary 
   #                 aes(label = summe), 
   #                position = position_fill(vjust = 0.5), 
   #                direction = "x", 
   #                 color = "black") +
    geom_text(data = total_n_flag2,
              aes(x = construct, y = 1, label = paste0("n = ",n)),  # Position at the top of each bar
              vjust = -0.5,  # Adjust vertical position of the text
              color = "black",
              inherit.aes = FALSE) +  # Avoid inheriting aesthetics from the main plot
    theme_bw(base_size = 14)

cvdPlot(simp_plot_flag2) # checks plot for colorblind readibility 

# plot the elaborated categories 
elab_plot_flag2 <- 
ggplot(flag2_plot_data, aes(x = construct, y = summe, fill = motility_type)) +
  facet_wrap(~ origin_tissue, ncol = 2,
             labeller = as_labeller(c("Hemolymph" = "Haemolymph",
                                      "Salivary gland" = "Salivary gland"))) +
  geom_col(position = "fill", width = 0.5, color = "black", size = 0.3) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = seq(0, 100, 25)) + 
  labs(y = "Motility pattern [%]", x = "") + 
  scale_fill_manual(breaks = mt, 
                    labels = c("Gliding",
                               "Patch-gliding",
                               "Waving/flipping",
                               "Flexing",
                               "Non motile",
                               "Floating"),
                    name = "Motility pattern",
                    values = blue_pastel_palette) + 
  scale_x_discrete(labels = c("PbANKA" = expression(italic("wt")), 
                              "TRP1flag2" = expression(italic("trp1-flag-gfp")))) +
  # geom_text(data = flag2_plot_data %>% filter(summe > 0), # absolute count for each category, removed for results, kept in supplementary 
  #           aes(label = summe), 
  #           position = position_fill(vjust = 0.5), 
  #           #direction = "x", 
  #           color = "black") +
  geom_text(data = total_n_flag2,
            aes(x = construct, y = 1, label = paste0("n = ",n)),  # Position at the top of each bar
            vjust = -0.5,  # Adjust vertical position of the text
            color = "black",
            inherit.aes = FALSE) +  # Avoid inheriting aesthetics from the main plot
  theme_bw(
    base_size = 14)


cvdPlot(elab_plot_flag2) # checks plot for colorblind readibility 


####### plot all TSR-related constructs #######
# create dataframe for the TSR constructs 
tsr_data <- longdf %>% 
  filter(construct %in% c("TRP1_delTSR", "PbANKA", "TRP1_TSRswap", "TRP1_TSRpoint")) %>% 
  group_by(origin_tissue, construct, motility_type, motility_type_simp) %>% 
  mutate(count1 = as.numeric(count1)) %>%
  summarise(summe = sum(count1)) %>% 
  mutate(motility_type = factor(motility_type, levels = mt),
         construct = factor(construct, levels = order_tsr)) # bring the types in order for the plot

# filter number of sporozoites in tsr groups and wt
total_n_tsr <- total_n %>% 
  filter(construct %in% c("TRP1_delTSR", "PbANKA", "TRP1_TSRswap", "TRP1_TSRpoint")) 


# Now apply the group_by and summarise 
tsr_data_simp <- tsr_data %>% 
  group_by(origin_tissue, construct, motility_type_simp) %>% 
  summarise(summe_simp = sum(summe, na.rm = TRUE)) %>%  # Add na.rm = TRUE to ignore NAs
  mutate(construct = factor(construct, levels = order_tsr),
         motility_type_simp = factor(motility_type_simp, levels = simp))

####### plots TSR #########
# plot for tsr with simplified motility categories 
simp_plot_tsr <- 
  ggplot(tsr_data_simp , aes(x = construct, y = summe_simp, fill = motility_type_simp)) +
  facet_wrap(~ origin_tissue, ncol = 2, scales = "free_x",
             labeller = as_labeller(c("Hemolymph" = "Haemolymph",
                                      "Salivary gland" = "Salivary gland"))) +
  geom_col(position = "fill", width = 0.5, color = "black", size = 0.3) +
 # geom_col(position = "fill", width = 0.5) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = seq(0, 100, 25)) + 
  labs(y = "Motility pattern [%]",
       x = "") + 
  scale_fill_manual(breaks = simp, 
                    labels = c("Gliding",
                               "Non productive",
                               "Floating"),
                    name = "Motility pattern",
                    values = bluesimp) + 
  scale_x_discrete(labels = c(
    "PbANKA" = expression(italic("wt")), 
    "TRP1_TSRswap" = expression(italic("trp1-tsr-swap")), 
    "TRP1_delTSR" = expression(italic("trp1-") * Delta * italic("tsr")),
    "TRP1_TSRpoint" = expression(italic("trp1-tsr-point"))
  ))+
 # geom_text(data = tsr_data_simp %>% filter(summe_simp > 0),  # n for all groups (removed for results kept for supplementa)
        #    aes(label = summe_simp), 
         #   position = position_fill(vjust = 0.5), 
            #direction = "x", 
        #    color = "black") +
  geom_text(data = total_n_tsr,
            aes(x = construct, y = 1, label = paste0("n = ",n)),  # Position at the top of each bar
            vjust = -0.5,  # Adjust vertical position of the text
            color = "black",
            inherit.aes = FALSE) +  # Avoid inheriting aesthetics from the main plot
  theme_bw(
    base_size = 14
  )+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) #+  # Ensure label tilt with proper adjustments


simp_plot_tsr

# Plot with facet_grid and ensuring labels tilt
# plot elaborate categories for tsr constructs 
tsr_mot_elab <- ggplot(tsr_data, aes(x = construct, y = summe, fill = motility_type)) +
  facet_grid(~ origin_tissue, scales = "free_x",
             labeller = as_labeller(c("Hemolymph" = "Haemolymph",
                                    "Salivary gland" = "Salivary gland"))) +  # Free the x-axis scale in facet_grid
  geom_col(position = "fill", width = 0.5, color = "black", size = 0.3) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = seq(0, 100, 25)) + 
  labs(y = "Motility pattern [%]",
       x = "") + 
   scale_x_discrete(labels = c(
    "PbANKA" = expression(italic("wt")), 
    "TRP1_TSRswap" = expression(italic("trp1-tsr-swap")), 
    "TRP1_delTSR" = expression(italic("trp1-") * Delta * italic("tsr")),
    "TRP1_TSRpoint" = expression(italic("trp1-tsr-point"))
  ))+
  scale_fill_manual(breaks = mt, 
                    labels = c("Gliding",
                               "Patch-gliding",
                               "Waving/flipping",
                               "Flexing",
                               "Non motile",
                               "Floating"),
                    name = "Motility pattern",
                    values = blue_pastel_palette) + 
 # geom_text(data = tsr_data %>% filter(summe > 0), # absolute count in tsr constructs
#            aes(label = summe),
#            position = position_fill(vjust = 0.5),
#            color = "black") +
  geom_text(data = total_n_tsr,
            aes(x = construct, y = 1, label = paste0("n = ",n)),  # Position at the top of each bar
            vjust = -0.5,  # Adjust vertical position of the text
            color = "black",
            inherit.aes = FALSE) +  # Avoid inheriting aesthetics from the main plot
  theme_bw(base_size = 14)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) #+  # Ensure label tilt with proper adjustments
  

# Export figures as PDF

pdf(file = "D:/Documents/UNI/MASTER/THESIS/Figures Thesis/Exported14012025/mot_tsr_elab.pdf",
    title = "mot_tsr_elab",
    height = 7,
    width = 8)
print(tsr_mot_elab)
dev.off()

pdf(file = "D:/Documents/UNI/MASTER/THESIS/Figures Thesis/Exported14012025/mot_tsr_simp.pdf",
    title = "mot_tsr_simp",
    height = 7,
    width = 8)
print(simp_plot_tsr)
dev.off()

pdf(file = "D:/Documents/UNI/MASTER/THESIS/Figures Thesis/Exported14012025/mot_flag_simp.pdf",
    title = "mot_flag_simp",
    height = 7,
    width = 8)
print(simp_plot_flag2)
dev.off()

pdf(file = "D:/Documents/UNI/MASTER/THESIS/Figures Thesis/Exported14012025/mot_flag_elab.pdf",
    title = "mot_flag_elab",
    height = 7,
    width = 8)
print(elab_plot_flag2)
dev.off()
