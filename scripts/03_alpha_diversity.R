################################################################################
#          _______  _____  _____ ___ ____  _   _ __  __ _____ ____  
#         | ____\ \/ / _ \|  ___|_ _/ ___|| | | |  \/  | ____|  _ \ 
#         |  _|  \  / | | | |_   | |\___ \| |_| | |\/| |  _| | | | |
#         | |___ /  \ |_| |  _|  | | ___) |  _  | |  | | |___| |_| |
#         |_____/_/\_\___/|_|   |___|____/|_| |_|_|  |_|_____|____/ 
#
################################################################################
#
# R script to analyze apha diversity results
#
# Arthur Escalas September 2020
# arthur.escalas@gmail.com
################################################################################


################################ LOAD DATA #####################################


for (keep_F1 in c(TRUE, FALSE)) {
  
  if (keep_F1 == TRUE) {
    today <- "with_F1"
    dir_save <- paste0(dir_alpha_div, today, "/")
    dir.create(dir_save)
    df_data <- read.csv(paste0(dir_data, 
                               "data_alpha_diversity.csv")) 
  } else {
    today <- "without_F1"
    dir_save <- paste0(dir_alpha_div, today, "/")
    dir.create(dir_save)
    df_data <- read.csv(paste0(dir_data, 
                               "data_alpha_diversity.csv")) %>% 
      filter(treatment_tank != "F1")
  }


# names of the indexes ----

index_alpha <- c("taxo_q0", "taxo_q1","phylo_q0", "phylo_q1")
names(index_alpha) <- index_alpha



###################### SUMMARIZE DIVERSITY ESTIMATES ###########################


# create a directory for the plots 

dir_plot <- paste0(dir_save, "plots_alpha_diversity/")
dir.create(dir_plot)


# ---------------- Compare the different compartments --------------------------


png(paste0(dir_plot, "plot_alpha_all_compartments_phylo_q0q1.png"), height = 15, 
    width = 25,
    unit = "cm", res = 200)
  par(mfrow = c(1,2), mar = c(3,4,1,1), mgp = c(2.5,0.5,0), las = 1, 
      font.lab = 2, oma = c(2,1,1,1))
  boxplot(df_data[, "phylo_q0"] ~ df_data$sample_type, main = "", xlab = "", 
          ylab = "Phylogenetic richness", drop = TRUE, outline = FALSE)
  boxplot(df_data[, "phylo_q1"] ~ df_data$sample_type, main = "", xlab = "", 
          ylab = "Phylogenetic diversity", drop = TRUE, outline = FALSE)
  mtext("Compartment of the mesocosm", 1, outer = TRUE, font = 2, line = 0)
dev.off()


# summarize diversity values

out <- df_data %>% group_by(sample_type) %>% 
  summarize_at(vars(contains("_q")), .funs = list(~mean(.), ~sd(.))) %>% data.frame()

write.csv(out, paste0(dir_save, "table_meansd_alpha_per_compartment.csv"), row.names = FALSE)


# ---------------------- Compare water samples ---------------------------------

dat <- df_data %>% filter(sample_type == "Water")

# estimate average and sd diversity at d0 ----
summ_d0 <- dat %>% filter(sampling_day == 0) %>% 
  group_by(treatment_tank) %>% 
  dplyr::summarize_at(index_alpha, mean) %>% 
  dplyr::summarize_at(index_alpha, list(avg = mean, sd = sd)) %>% data.frame()

# estimate average and sd diversity at d8 ----
summ_d8 <- dat %>% filter(sampling_day == 8) %>% 
  group_by(treatment, treatment_tank) %>% 
  dplyr::summarize_at(index_alpha, mean) %>% 
  group_by(treatment) %>%
  dplyr::summarize_at(index_alpha, list(avg = mean, sd = sd)) %>% data.frame()

# combine all data to plot ----
out <- rbind(c(treatment = "D0", summ_d0), summ_d8)

write.csv(out, paste0(dir_save, "table_meansd_alpha_water.csv"), row.names = FALSE)


# ---------------------- Compare SEDIMENT samples ------------------------------

# all sediments ----

dat <- df_data %>% filter(sample_type == "Sediment" & 
                            position_in_the_tank == "tank" & 
                            sampling_day == 8)


# estimate average and sd diversity at d0 ----
summ_d0 <- df_data %>% filter(sample_type == "Sediment" & 
                                sampling_day == 0) %>% 
  # group_by(treatment_tank) %>% 
  dplyr::summarize_at(index_alpha, list(avg = mean, sd = sd)) %>% data.frame()

# estimate average and sd diversity at d8 ----
summ_d8 <- dat %>% filter(sampling_day == 8) %>% 
  group_by(treatment, treatment_tank) %>% 
  dplyr::summarize_at(index_alpha, mean) %>% 
  group_by(treatment) %>%
  dplyr::summarize_at(index_alpha, list(avg = mean, sd = sd)) %>% data.frame()

# combine all data to plot ----
out <- rbind(c(treatment = "D0", summ_d0), summ_d8)

write.csv(out, paste0(dir_save, "table_meansd_alpha_sediment.csv"), row.names = FALSE)


# different positions in tank at D8 ----

dat <- df_data %>% filter(sample_type == "Sediment" & sampling_day == 8 & 
                            position_in_the_tank != "tank")

# summarize diversity values

summ_d0 <- df_data %>% filter(sample_type == "Sediment" & 
                                sampling_day == 0 & tank_id == "Bucket") %>% 
  dplyr::summarize_at(index_alpha, list(avg = mean, sd = sd)) %>% data.frame()

summ_d8 <- dat %>% filter(position_in_the_tank == "open_shelter") %>% 
  group_by( treatment) %>% 
  summarize_at(vars(index_alpha), .funs = list(avg = mean, sd = sd)) %>% data.frame()

# combine all data to plot ----
out <- rbind(c(#position_in_the_tank = NA, 
               treatment = "D0", summ_d0), summ_d8)

write.csv(out, paste0(dir_save, "table_meansd_alpha_sediment_per_tank_position.csv"), row.names = FALSE)


# ------------------------- Compare TURF samples --------------------------------

dat <- df_data %>% filter(sample_type == "Turf")

# estimate average and sd diversity at d0 ----
summ_d0 <- dat %>% filter(sampling_day == 0) %>% 
  # group_by(treatment) %>% 
  dplyr::summarize_at(index_alpha, list(avg = mean, sd = sd)) %>% data.frame()

# estimate average and sd diversity at d8 ----
summ_d8 <- dat %>% filter(sampling_day == 8) %>% 
  group_by(treatment, treatment_tank) %>% 
  dplyr::summarize_at(index_alpha, mean) %>% 
  group_by(treatment) %>%
  dplyr::summarize_at(index_alpha, list(avg = mean, sd = sd)) %>% data.frame()

# combine all data to plot ----
out <- rbind(c(treatment = "D0", summ_d0), summ_d8)

write.csv(out, paste0(dir_save, "table_meansd_alpha_turf.csv"), row.names = FALSE)




############################## TEST DIFFERENCES ################################

dir_stat <- paste0(dir_save, "results_stats/")
dir.create(dir_stat)

# ============================ MIXED MODELS ====================================

ls_lmm <- list()

# --------------------------- define and fit models ----------------------------

#  compartments ----

model_list <- list(mod1 = " ~ sample_type + (1| sample_type : treatment_tank)")

input_data <- df_data
fits  <- lapply(index_alpha, function(idx) {
  fit_lmm(input_data, model_list, nm_var_Y = idx, reml = TRUE)
}) %>% setNames(index_alpha)

ls_lmm$compartments <- fits


#  water ----

model_list <- list(mod1 = " ~ treatment + (1| treatment : treatment_tank )")

input_data <- df_data %>% filter(sample_type == "Water" & sampling_day == 8)

fits  <- lapply(index_alpha, function(idx) {
  fit_lmm(input_data, model_list, nm_var_Y = idx, reml = TRUE)
}) %>% setNames(index_alpha)

ls_lmm$water <- fits


#  Sediment ----

model_list <- list(mod1 = " ~ treatment + (1| treatment : treatment_tank )")

input_data <- df_data %>% filter(sample_type == "Sediment" & sampling_day == 8 &
                                   position_in_the_tank == "tank")

fits  <- lapply(index_alpha, function(idx) {
  fit_lmm(input_data, model_list, nm_var_Y = idx, reml = TRUE)
}) %>% setNames(index_alpha)

ls_lmm$sediment <- fits


#  Turf ----

model_list <- list(mod1 = " ~ treatment + (1| treatment : treatment_tank )")

input_data <- df_data %>% filter(sample_type == "Turf" & sampling_day == 8)

fits  <- lapply(index_alpha, function(idx) {
  fit_lmm(input_data, model_list, nm_var_Y = idx, reml = TRUE)
}) %>% setNames(index_alpha)

ls_lmm$turf <- fits


#  Open shelters between treatments ----

model_list <- list(mod1 = " ~ treatment + (1| treatment : treatment_tank )")

input_data <- df_data %>% filter(sample_type == "Sediment" & sampling_day == 8 &
                                   position_in_the_tank == "open_shelter")

fits  <- lapply(index_alpha, function(idx) {
  fit_lmm(input_data, model_list, nm_var_Y = idx, reml = TRUE)
}) %>% setNames(index_alpha)

ls_lmm$open_shelter <- fits


#  Open vs closed shelters ----

model_list <- list(mod1 = " ~ treatment * position_in_the_tank + (1| treatment : treatment_tank )")

input_data <- df_data %>% filter(sample_type == "Sediment" & sampling_day == 8 &
                                   position_in_the_tank != "tank")

fits  <- lapply(index_alpha, function(idx) {
  fit_lmm(input_data, model_list, nm_var_Y = idx, reml = TRUE)
}) %>% setNames(index_alpha)

ls_lmm$open_vs_closed_shelter <- fits

saveRDS(ls_lmm, file = paste0(dir_stat, "res_lmm.rds"))


# ------------------------- get the results of anova ---------------------------

df_anova <- lapply(ls_lmm, function(X) {
  lapply(X, function(x) {
    out <- x$anova$mod1[,5:6] %>% data.frame()
    names(out) <- c("F-value", "p-value")
    round(out, 4) %>% rownames_to_column("factor")
  }) %>% reformat_as_df(new_var_name = "yvar")
}) %>% reformat_as_df(new_var_name = "compartment") %>% 
  separate(yvar, into = c("facet","order"), sep = "_", remove = FALSE) %>% 
  dplyr::select(compartment, everything())

write.csv(df_anova, file = paste0(dir_stat, "table_res_lmm_ANOVA.csv"), row.names = FALSE)



############################## FINAL PLOTS #####################################


# load the data for plot ----

nms_data <- c("water","sediment","turf","sediment_per_tank_position")
ls_dat <- lapply(nms_data, function(x) {
  read.csv(paste0(dir_save, "table_meansd_alpha_", x, ".csv"))
}) %>% setNames(c("Water", "Sediment", "Turf", "Open_shelter"))

nms_index_to_plot <- c("Taxonomic richness", "Taxonomic diversity", "Taxonomic diversity",
                       "Phylogenetic richness", "Phylogenetic diversity", "Phylogenetic diversity") %>% 
  setNames(index_alpha)


for (idx in index_alpha) {
  
  nm_yvar <- nms_index_to_plot[idx]
  var_avg <- paste0(idx, "_avg")
  var_sd <- paste0(idx, "_sd")

  
  png(paste0(dir_plot, "plot_alpha_diversity_", idx,".png"),
      height = 15, width = 20, res = 400, unit = "cm")
  
  par(mfrow = c(2,2), oma = c(4,1,1,0), xpd = NA, mar = c(3,1,1,1))
  
  tmp <- do.call(rbind, ls_dat)
  ylim <- round(setrge(c(min(tmp[, var_avg]) - max(tmp[, var_sd]),
                         max(tmp[, var_avg]) + max(tmp[, var_sd])), 20),0)
  
  for (nm_dat in names(ls_dat)) {
    
    dat <- ls_dat[[nm_dat]]
    if (min(ylim) < 0) {ylim[1] <- 0}
    cols <- c(Control = "white", Fish = "black")
    cols_border <- c(Control = "black", Fish = "black")
    
    if (nm_dat %in% names(ls_dat)[3:4]) {
      nmx <- nms_treatment
    } else {
      nmx <-  ""
    }
    niceplot(limX=c(0.75,1.75),limY=ylim, marg=c(1,3,1,1), cexY = 0.7, cexX = .9,
             tick=-0.4,   nlab=4, labX=c(1,1.5), nmlabX=nmx, nmX = "",
             nmY = "", lineX = 1,cexXt = 0.8,
             lasX=1,lasY=1, cexYt = 0.7, lineYt = 2)
    text(1.25, max(ylim) + 0.05*diff(ylim), labels = gsub("_", " ", nm_dat), 
         font = 4, cex = 1.2)
    mask_d0 <- dat$treatment == "D0"
    segments(x0 =0.75, y0 = dat[mask_d0, var_avg], x1 = 1.75, 
             y1 = dat[mask_d0, var_avg],lty = 1, col = "grey")
    segments(x0 =0.75, y0 = dat[mask_d0, var_avg] - dat[mask_d0, var_sd], x1 = 1.75, 
             y1 = dat[mask_d0, var_avg] - dat[mask_d0, var_sd],lty = 3, col = "grey")
    segments(x0 =0.75, y0 = dat[mask_d0, var_avg] + dat[mask_d0, var_sd], x1 = 1.75, 
             y1 = dat[mask_d0, var_avg] + dat[mask_d0, var_sd],lty = 3, col = "grey")
    
    meansey(c(1,1.5), dat[!mask_d0, var_avg], dat[!mask_d0, var_sd], pchp = 22, 
            colp = cols_border, bgp = cols, cexp = 2, colb = cols_border, 
            lgt = 0.01)
    # pvalues
    tmp <- df_anova %>% filter(yvar == idx & compartment == tolower(nm_dat))
    text(1.25, max(ylim) - 0.1 * diff(ylim), font = 3,
         labels = paste0("F-value = ", round(tmp$`F-value`, 1),
                         "   p-value = ", round(tmp["p-value"], 3)))  
  }
  mtext("Experimental treatment", 1, outer = TRUE, font = 2, line = 2)
  mtext(nm_yvar, 2, outer = TRUE, font = 2, line = -1)
  
  dev.off()
  
  
}

}
