#!/usr/bin/env Rscript

#------------------------------------------------------------------------------
# From the average intensity data per z-level, an average per individual/well
# was calculated and extreme outliers (3x interquartile range) were removed
# from the statistical analysis. 
#
# each plate is the same experiment but at different time points
# one plate (output excel) = one batch, random effect
# DAPI, Cy3, FITC = log(intensity)
# Control, 12.5, 50, 100, 200 = group (treatment)
#
# Linear mixed-effect model (LMM) works well for observations that are NOT
# independent (clustered / grouped), in this case, repeated measurements in
# different wells across a row.
# Linear regression assumes all observations are independent.
# Experiment_ID = batch

#------------------------------------------------------------------------------


suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(readxl)
  library(ggplot2)
  library(tidyr)
  library(rstatix)
  library(data.table)
  library(stringr)
  library(DescTools)
  library(patchwork)
  library(DHARMa)    # lmm diagnostics
  library(car)       # lmm Anova
  library(emmeans)   # lmm Post-Hoc
  library(lme4)      # lmm
})


#---- Global Publication Theme ------------------------------------------------


theme_set(theme_bw(base_size = 14))
pub_theme <- theme(
  plot.title     = element_text(face = "bold", size = 18, hjust = 0.5),
  plot.subtitle  = element_text(size = 12, hjust = 0.5),
  axis.title     = element_text(face = "bold", size = 16),
  axis.text      = element_text(color = "black", size = 13),
  axis.text.x    = element_text(angle = 45, hjust = 1),
  strip.text     = element_text(face = "bold", size = 14),
  legend.title   = element_text(face = "bold", size = 14),
  legend.text    = element_text(size = 13)
)
control_group_name <- "Control"
color_palette <- c("DAPI" = "#3a86ff", "Cy3" = "#ff006e", "FITC" = "#38b000")


#---- Import files ------------------------------------------------------------


# 1. READ RAW DATA (and tag with filename as Experiment_ID)
file.list <- list.files(pattern = '*.xls')

# Remove file extension
names(file.list) <- tools::file_path_sans_ext(file.list)

rawdf <- rbindlist(lapply(file.list, read_excel), idcol = "Experiment_ID")

# The size of a daphnia neonate is >10000 pixels, anything below is removed
# avoid unintended analysis of artefacts.
rows <- rawdf$`1_area` >= 10000
df_filtered <- rawdf[rows, ]

cat(sprintf(
  " > Removed observations w/ weak signal: %s of %s (%.2f%%)\n",
  sum(!rows),
  scales::comma(length(rows)),
  sum(!rows) / length(rows) * 100
))

# rename columns
RD <- df_filtered[, c('Experiment_ID', 'well_name', 'biggestobject_1_mean', 'biggestobject_2_mean', 'biggestobject_3_mean')]
colnames(RD) <- c('Experiment_ID', 'Well', 'Cy3', 'DAPI', 'FITC')


#------------------------------------------------------------------------------
# 2. EXTRACT ROW LETTER 
# Wells named "A01", "B03" etc => "A", "B"
RD <- RD %>%
  mutate(Row = substr(Well, 1, 1))

# Get average per rows (A01, A02, etc) for each experiment
# x <- RD[Experiment_ID == "20260109-DM-MPLX-CCCP-24h_measures" & Well == "A03" & Row == "A"]
# mean(x$Cy3, na.rm = TRUE)
Well_Summaries <- RD %>%
  group_by(Experiment_ID, Well, Row) %>%
  summarise(across(c(DAPI, Cy3, FITC), \(x) mean(x, na.rm = TRUE)), .groups = "drop")

# 3. READ METADATA AND MERGE
# Load the CSV you created
plate_maps <- read.csv("Plate_Layouts.csv")

# Check what IDs exist in your R data vs your CSV
unique(Well_Summaries$Experiment_ID)
unique(plate_maps$Experiment_ID)

# Test if they intersect
any(unique(Well_Summaries$Experiment_ID) %in% unique(plate_maps$Experiment_ID))

# Merge the biological data with the experimental design
RD_final <- Well_Summaries %>%
  left_join(plate_maps, by = c("Experiment_ID", "Row")) %>%
  filter(group != '' & group != 'Blank' & !is.na(group))

#----
# why is Blank removed?

#----


#==== refactor ====
#groups <- unique(RD_final$group) # includes all in order

# 2. Filter the data to only include these groups, and set the factor levels
#RD_final <- RD_final %>%
#  #filter(group %in% groups) %>%
#  mutate(group = factor(group, levels = groups))

#====

# --- MANUAL GROUP SELECTION & ORDERING ---

# 1. Define exactly which groups to keep, in the EXACT order you want them plotted.
# -> Replace these example strings with your actual group names.
# -> CRITICAL: Your control_group_name ("DMSO 0.2%") MUST be in this list!
desired_groups <- c(
  "Control",         # Usually put your vehicle control first
  "12.5",      # Replace with your actual text
  "50",
  "100",
  "200"   # Add as many as you need, dropping the ones you don't
)

# 2. Filter the data to only include these groups, and set the factor levels
RD_final <- RD_final %>%
  filter(group %in% desired_groups) %>%
  mutate(group = factor(group, levels = desired_groups))

# --- 3. PREPARE DATA STREAMS ---

# A. Long Format
data_long_raw <- RD_final %>%
  pivot_longer(cols = c("DAPI", "Cy3", "FITC"), names_to = "channel", values_to = "intensity") %>%
  mutate(channel = factor(channel, levels = c("DAPI", "Cy3", "FITC"))) %>%
  mutate(obs_id = row_number())

# B. Identify Extreme Outliers
outliers_found <- data_long_raw %>%
  group_by(channel, group) %>%
  identify_outliers(intensity) %>%
  filter(is.extreme == TRUE) 

# C. Clean Data (For Statistics & Boxplot)
data_clean <- data_long_raw %>%
  anti_join(outliers_found, by = "obs_id")


#---- LMM STATISTICS MODULE (Mixed-Effects) -----------------------------------

lmm_stats_list <- list()
lmm_models <- list() # Store models for DHARMa diagnostics
anova_results_list <- list() # Store overall ANOVA results

# Find the numeric factor level of the control group for Dunnett's
control_index <- which(levels(data_clean$group) == control_group_name)

for (chan in c("DAPI", "Cy3", "FITC")) {
  chan_data <- data_clean %>%
    filter(channel == chan) %>%
    drop_na(intensity, group, Experiment_ID)
  
  # 1. Fit Log-Normal Mixed Model & Store it
  # We use log(intensity) with a standard lmer to handle fluorescence
  # distribution
  fit_log_mixed <- lmer(
    log(intensity) ~ group + (1 | Experiment_ID),
    #log(intensity) ~ group + batch,
    data = chan_data,
    REML = TRUE
  )

  #==== Sanity check ====
  # expect boxes around zero with similar spread, 
  # no batch/plate should be off-zero
  #======================
  chan_data$resid <- residuals(fit_log_mixed)
  boxplot(resid ~ Experiment_ID, data = chan_data)

  #==== linear model ====
  data <- chan_data
  data$Experiment_ID <- as.character(as.numeric(as.factor(data$Experiment_ID)))
  res <- lm(log(intensity) ~ group * Experiment_ID, data)
  plot(log(intensity) ~ group, data)
  #======================


  lmm_models[[chan]] <- fit_log_mixed
  
  # Calculate and print the overall ANOVA (Analysis of Deviance)
  # Note: For lmmer, car::Anova defaults to Wald Chisq tests instead of F tests
  model_anova <- car::Anova(fit_log_mixed, type = "II")

  #==== F test is better for smaller sample sizes? Like 3 batches? ====
  #model_anova <- car::Anova(fit_log_mixed, type = "II", test.statistic = "F")
  #====
  
  cat(paste("\n--- ANOVA for", chan, "---\n"))
  print(model_anova)
  
  # Format and store the ANOVA table
  anova_df <- as.data.frame(model_anova)
  anova_df$channel <- chan
  anova_df$Term <- rownames(anova_df)
  anova_df <- anova_df %>% select(channel, Term, everything())
  anova_results_list[[length(anova_results_list) + 1]] <- anova_df
  
  # 2. Run Dunnett's Post-Hoc via emmeans
  chan_means <- emmeans(fit_log_mixed, ~ group, type = "response")
  dunnett_res <- contrast(chan_means, method = "trt.vs.ctrl", ref = control_index)

  #==== validation ============================================================
  # Because you used type = "response" on a log-scale model, the contrasts
  # come back as ratios, not differences:
  # 1.50 = 50% increase; 0.50 = 50% decrease
  #============================================================================
  
  # 3. Convert results to dataframe and parse out the treatment names
  res_df <- as.data.frame(dunnett_res)
  
  # Clean up the contrast column to just get the target group name
  res_df <- res_df %>%
    mutate(
      group = str_remove(contrast, paste0(" / ", control_group_name)),

      # Fallback just in case
      group = str_remove(group, paste0(" - ", control_group_name)), 
      
      channel = chan,
      p_val = p.value, # <- redundent? 2 columns with pval
      label = case_when(
        p_val <= 0.001 ~ "***", 
        p_val <= 0.01  ~ "**", 
        p_val <= 0.05  ~ "*", 
        TRUE           ~ ""
      )
    ) %>%
    # Keep only significant labels to map to plots
    filter(label != "")
  
  lmm_stats_list[[length(lmm_stats_list) + 1]] <- res_df
}

# Bind all channels into the master significance dataframe for the plots
final_sig_df <- if(length(lmm_stats_list) > 0) bind_rows(lmm_stats_list) else data.frame()
final_anova_df <- bind_rows(anova_results_list)
plot_subtitle <- "Significance: Log-Normal LMM (Mixed Effects) + Dunnett's Test"


#---- DHARMa DIAGNOSTIC PLOTS -------------------------------------------------
# 1. Open a high-resolution PNG file
png(
  filename = "FigureS1_DHARMa_Diagnostics.png", 
  width = 12, height = 15, units = "in", res = 300
)

# 2. Set up grid canvas: Increased top margin for extra headroom
par(mfrow = c(3, 3), mar = c(4, 4, 4.5, 2) + 0.1)

# 3. Loop through saved models and plot
for (chan in c("DAPI", "Cy3", "FITC")) {
  
  # Fetch the correct data for this specific channel
  chan_data <- data_clean %>%
    filter(channel == chan) %>%
    drop_na(intensity, group)
  
  # Pull the stored model for the current channel
  fit_log <- lmm_models[[chan]]
  
  # Simulate residuals
  sim_res <- simulateResiduals(fittedModel = fit_log, plot = FALSE)
  
  # Print outliers to console for each channel
  # DHARMa outliers are defined as scaled residuals of exactly 0 or 1
  outlier_indices <- which(sim_res$scaledResiduals == 0 | sim_res$scaledResiduals == 1)
  
  if(length(outlier_indices) > 0) {
    warning(paste(
      "In", chan, "- Residuals were outside the simulation envelope (DHARMa outliers):",
      "\n",
      paste(capture.output(print(chan_data[outlier_indices, ])), collapse = "\n")
    ))  
  } else {
    cat(paste("No DHARMa outliers detected for", chan, "\n"))
  }
  
  # Plot 1: QQ Plot (Left Column)
  plotQQunif(sim_res)
  mtext(paste(chan, "- QQ Plot"), side = 3, line = 2.5, font = 2, adj = 0, cex = 1.2)
  
  # Plot 2: Residuals vs Categorical Predictor (Middle Column)
  plotResiduals(sim_res, form = chan_data$group)
  mtext(paste(chan, "- Residuals"), side = 3, line = 2.5, font = 2, adj = 0, cex = 1.2)
  
  # Plot 3: Dispersion Test (Right Column)
  testDispersion(sim_res)
  mtext(paste(chan, "- Dispersion Test"), side = 3, line = 2.5, font = 2, adj = 0, cex = 1.2)
}

# 4. Close the device and reset canvas
dev.off()
par(mfrow = c(1, 1))
print("Diagnostic plots successfully saved as 'FigureS1_DHARMa_Diagnostics.png'")


#---- PLOT GENERATION (Normalized for Batch Effects) --------------------------
# 1. Normalize data to Vehicle Control within each experiment
plot_data <- data_clean %>%
  group_by(Experiment_ID, channel) %>%
  mutate(
    # Find the mean of the control group for this specific experiment and channel
    control_mean = mean(intensity[group == control_group_name], na.rm = TRUE),
    # Calculate fold change
    fold_change = intensity / control_mean
  ) %>%
  ungroup()

# Define Facet Names for Scientific Labels
facet_names <- c(
  "Cy3"  = "MXR system (Rhodamine B)",
  "DAPI" = "Mitochondrial Membrane Potential\n(BioTracker Blue)",
  "FITC" = "Esterase Activity (Calcein AM)"
)

# 1. Define Highlight Settings
alpha_highlight <- c("DAPI" = 0.9, "Cy3" = 0.2, "FITC" = 0.2)
# 2. Prepare Data for Labels (Using Fold Change max values)
max_y_values <- plot_data %>%
  group_by(channel) %>%
  summarise(max_y = max(fold_change, na.rm = TRUE), .groups = 'drop')

# Survival text positioning (Calculated directly from raw metadata)
survival_text <- plate_maps %>%
  # 1. Filter to keep only the desired groups
  filter(group %in% desired_groups) %>%
  # 2. Set factor levels so the order matches the plot exactly
  mutate(group = factor(group, levels = desired_groups)) %>%
  # 3. Calculate the true mean survival across ALL initial wells
  group_by(group) %>%
  summarise(mean_survival = mean(survival, na.rm = TRUE), .groups = 'drop') %>%
  # 4. Create the text label, rounding to 0 decimal places (e.g., "S: 98%")
  mutate(label = paste0("S: ", round(mean_survival * 100, 0), "%")) %>%
  # 5. Cross with max Y values for plotting height
  crossing(max_y_values) %>% 
  mutate(y_pos_surv = max_y * 1.12)

# Map max Y coordinates to the lmmM significance df so the stars hover correctly
if(nrow(final_sig_df) > 0) {
  # Using base R merge to bypass factor/character join mismatches
  final_sig_df <- merge(final_sig_df, max_y_values, by = "channel", all.x = TRUE)
  final_sig_df$y_pos_raw <- final_sig_df$max_y * 1.05
}

# 3. DEFINE THE PLOTTING FUNCTION
create_plot <- function(target_channel, show_y_axis = TRUE) {
  
  # Filter Data
  sub_plot  <- plot_data %>% filter(channel == target_channel)
  sub_stats <- final_sig_df %>% filter(channel == target_channel)
  sub_surv  <- survival_text %>% filter(channel == target_channel)
  max_y     <- max(sub_plot$fold_change, na.rm = TRUE)
  
  # --- LOGIC: Target FITC, DAPI or Cy3 ---
  is_target  <- (target_channel == "DAPI") 
  
  # Styling Parameters
  font_style <- ifelse(is_target, "bold", "plain")
  surv_font_style <- ifelse(is_target, "bold.italic", "italic")
  frame_width <- ifelse(is_target, 2.0, 0.5)
  tick_width <- ifelse(is_target, 1.2, 0.5)
  box_linewidth <- ifelse(is_target, 1.2, 0.5)
  current_alpha <- alpha_highlight[[target_channel]]
  current_title <- facet_names[[target_channel]]
  
  current_title <- facet_names[[target_channel]]
  
  p <- ggplot() +
    # Add a reference line at 1.0 (Control Baseline)
    geom_hline(yintercept = 1.0, linetype = "dashed", color = "gray60", linewidth = 0.8) +
    
    # Boxplot (using fold_change)
    geom_boxplot(data = sub_plot, 
                 aes(x = group, y = fold_change), 
                 fill = color_palette[target_channel], 
                 alpha = current_alpha, 
                 linewidth = box_linewidth,
                 outlier.shape = NA) +
    
    # Jittered Points (Colored by Experiment_ID)
    geom_jitter(data = sub_plot, 
                aes(x = group, y = fold_change, color = Experiment_ID, shape = Experiment_ID), 
                width = 0.2, 
                size = 2.5,
                alpha = current_alpha) + 
    
    # High-contrast colorblind-friendly palette for the dots
    scale_color_viridis_d(option = "magma", end = 0.8) +
    
    # Text: Stats Stars (Using the merged y_pos_raw)
    geom_text(data = sub_stats, aes(x = group, y = y_pos_raw, label = label), 
              fontface = "bold", size = 8) + 
    
    # Text: Survival Rates
    geom_text(data = sub_surv, aes(x = group, y = max_y * 1.15, label = label), 
              fontface = "italic", size = 3.5) +
    
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
    
    labs(title = current_title, 
         x = "µg/L", 
         y = if(show_y_axis) "Fold Change (vs Control)" else "",
         color = "Replicate Run:",
         shape = "Replicate Run:") + # Consolidate legends
    
    # Global Theme 
    pub_theme + 
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 14, face = font_style),
      axis.title.x = element_text(face = font_style, size = 14), 
      axis.title.y = if(show_y_axis) element_text(face = font_style, size = 14) else element_blank(),
      axis.text.x  = element_text(angle = 45, hjust = 1, face = font_style, color="black", size=12),
      axis.text.y  = element_text(face = font_style, color="black", size=12),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = frame_width),
      axis.ticks = element_line(color = "black", linewidth = tick_width),
      axis.ticks.length = unit(0.2, "cm")
    )
  
  return(p)
}
# 4. GENERATE & COMBINE
p_cy3  <- create_plot("Cy3", show_y_axis = TRUE)
p_dapi <- create_plot("DAPI", show_y_axis = FALSE)
p_fitc <- create_plot("FITC", show_y_axis = FALSE)

# Combine with Patchwork (Fixed Operator Precedence!)
p1_combined <- (p_cy3 + p_dapi + p_fitc) +
  plot_layout(ncol = 3, guides = "collect") +
  plot_annotation(
    title = "Merged Replicate Analysis of CCCP exposure (Normalized)",
    subtitle = paste(plot_subtitle, "| S: Survival Rate"),
    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
  ) & 
  theme(legend.position = "bottom")

print(p1_combined)

# 5. Save Plot 
ggsave("Figure1_BoxPlot_FoldChange.png", p1_combined, width = 15, height = 8, dpi = 300)


#---- XPORT STATISTICAL TABLES & DATA -----------------------------------------
cat(" > Exporting statistical tables...\n")

# 1. Save Omnibus ANOVA Results (Analysis of Deviance from the Mixed Model)
if(nrow(final_anova_df) > 0) {
  write.csv(final_anova_df, "Stats_LMM_ANOVA_Results.csv", row.names = FALSE)
}

# 2. Save Dunnett's Post-Hoc Results (P-values and significance stars)
if(nrow(final_sig_df) > 0) {
  
  # Safely select and rename columns regardless of what emmeans named the test statistics
  export_sig_df <- final_sig_df %>%
    rename(Treatment = group, 
           P_Value = p_val,
           Significance = label) %>%
    # any_of() prevents crashes if "estimate" is actually called "ratio", etc.
    select(channel, Treatment, P_Value, Significance, 
           any_of(c("ratio", "estimate", "SE", "df", "t.ratio", "z.ratio")))
  
  write.csv(export_sig_df, "Stats_LMM_Dunnett_Results.csv", row.names = FALSE)
}

# 3. Save the Normalized Data (Highly useful if you ever need to plot in GraphPad Prism)
write.csv(plot_data, "Processed_Normalized_Plot_Data.csv", row.names = FALSE)

# 4. Save Outliers List (For your audit trail / methods section)
if(nrow(outliers_found) > 0) {
  write.csv(outliers_found, "Stats_Removed_Outliers.csv", row.names = FALSE)
}

#---- Done --------------------------------------------------------------------