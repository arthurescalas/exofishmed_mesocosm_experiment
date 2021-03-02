################################################################################
#  _______  _____  _____ ___ ____  _   _ __  __ _____ ____  
# | ____\ \/ / _ \|  ___|_ _/ ___|| | | |  \/  | ____|  _ \ 
# |  _|  \  / | | | |_   | |\___ \| |_| | |\/| |  _| | | | |
# | |___ /  \ |_| |  _|  | | ___) |  _  | |  | | |___| |_| |
# |_____/_/\_\___/|_|   |___|____/|_| |_|_|  |_|_____|____/ 
#
################################################################################
#
# R script to analyse nutrient dynamics in macroalgae, feces and fishes
#
# Arthur Escalas October 2020
# arthur.escalas@gmail.com
################################################################################


########################### LOAD DATA #######################################

for (keep_F1 in c(TRUE, FALSE)) {
  
  if (keep_F1 == TRUE) {
    dir_save <- paste0(dir_nutrients_dynamics, "with_F1/")
    dir.create(dir_save)
    df_data <- read.csv(paste0(dir_data, 
                               "data_nutrients_algae.csv")) 
  } else {
    dir_save <- paste0(dir_nutrients_dynamics, "without_F1/")
    dir.create(dir_save)  
    df_data <- read.csv(paste0(dir_data, 
                               "data_nutrients_algae.csv")) %>% 
      filter(treatment_tank != "F1")
  }
  

df_data$group <- factor(df_data$group, levels = c("Environmental control D0",
                                      "Control D8", "Fish D8"))

vars_nut <- c(names(df_data)[grep("pct_", names(df_data))], "NP")

# aggregate data per tank at D8 ----
# remove the environmental control at D0

df_agg <- df_data %>% filter(group != "Environmental control D0") %>% 
  group_by(taxonomy, group, treatment, treatment_tank) %>% 
  summarise_at(vars_nut, mean) %>% data.frame()

# compute mean and average values in each group ----

df_summ <- df_agg %>% group_by(taxonomy, treatment) %>% 
  summarise_at(vars_nut, list(~mean(., na.rm = TRUE), ~sd(., na.rm = TRUE)))

write.csv(df_summ, file = paste0(dir_save, 
                                  "table_avg_sd_nutrient_per_species.csv"),
          row.names = FALSE)



# objects for analyses ----

nms_for_plot <-  gsub(" ", "\n", levels(df_data$group))
nms_species <- df_data$taxonomy %>% unique() %>% as.character()
nms_species_for_plot <- gsub("_", " ", nms_species)
nms_species_for_plot <- gsub("Cystoseira", "Treptacantha", nms_species_for_plot)


########################### MIXED MODELS #######################################

# define the models ----
# mod1 include only treatment effect and does not account for time effect
# model have a random effect of tank within treatment

model_list <- list(mod1 = " ~ treatment + (1| treatment : treatment_tank )")


# fit the lmm ----
lmm_fits <-  lapply(nms_species, function(nm_sp) {
  lapply(vars_nut, function(v) {
  dat <- df_data %>% filter(taxonomy == nm_sp & group != "Environmental control D0")
    fit_lmm(dat, model_list, nm_var_Y = v, reml = TRUE)
  }) %>% setNames(vars_nut)
}) %>% setNames(nms_species)

saveRDS(lmm_fits, file = paste0(dir_save, "res_lmm.rds"))


# ANOVA on fixed effects ----

df_anova <- lapply(lmm_fits, function(X) {
  lapply(X, function(x) {
    x$anova$mod1 %>% data.frame() %>% round(4)
  }) %>% reformat_as_df(new_var_name = "yvar")
}) %>% reformat_as_df(new_var_name = "species")

write.csv(df_anova, file = paste0(dir_save, 
                                  "table_lmm_ANOVA.csv"),
          row.names = FALSE)


########################### PLOT THE DATA ######################################

cols_treatment <- c("grey", "white", "black") %>% 
  setNames(levels(df_data$group))
cols_border <- c("black","black"," grey40") %>% 
  setNames(levels(df_data$group))

nms_vars <- c("pct_N", "pct_P", "NP")
nms_vars_toplot <- c("% N", "% P", "N:P") %>% setNames(nms_vars)

# estimate avarrages and SD across tanks ----

df_data_plot <- df_data

avg_d0 <- df_data_plot %>% filter(sampling_day == 0) %>% group_by(taxonomy) %>% 
  dplyr::summarize_at(nms_vars, funs(mean(., na.rm = TRUE))) %>% data.frame() %>% 
  column_to_rownames("taxonomy")
sd_d0 <- df_data_plot %>% filter(sampling_day == 0) %>% group_by(taxonomy) %>% 
  dplyr::summarize_at(nms_vars, funs(sd(., na.rm = TRUE))) %>% data.frame() %>% 
  column_to_rownames("taxonomy")

tmp <- cbind(avg_d0, sd_d0 %>% setNames(paste0("sd_", names(sd_d0))))
write.csv(tmp, file = paste0(dir_save, "table_algae_content_d0.csv"))


df_plot <- df_data_plot %>% filter(sampling_day == 8)
tmp <- df_plot %>% group_by(treatment, treatment_tank, taxonomy) %>% 
  dplyr::summarize_at(nms_vars, funs(mean(., na.rm = TRUE)))
avgs <- tmp %>% group_by(treatment, taxonomy) %>% 
  summarize_at(nms_vars, mean, na.rm = TRUE) %>% data.frame()
sds <- tmp %>% group_by(treatment, taxonomy) %>% 
  summarize_at(nms_vars, sd, na.rm = TRUE) %>% data.frame()


cols <- c(control = "white", fish = "black")
cols_border <- c(control = "black", fish = "black")

png(paste0(dir_save, "plot_NP_concentration_d8.png"),
    height = 15, width = 20, res = 400, unit = "cm")
par(mfrow = c(3,3), oma = c(4,1,2,1), xpd = NA, mar = c(3,2,1,1))

  for (var in nms_vars) {
    # ylim <- round(setrge(df_data_plot[, var], 0.1),1)
    
    if (var == "pct_N") {
      nms_main <- nms_species_for_plot %>% setNames(nms_species)
      ylim <- c(0.5,2.5)
      nmsx <- rep("",2)
    } 
    if (var == "pct_P") {
      nms_main <-  rep("", length(nms_species)) %>% setNames(nms_species)
      ylim <- c(0.02,0.1)
      nmsx <- rep("",2)
    }
    if (var == "NP") {
      nms_main <-  rep("", length(nms_species)) %>% setNames(nms_species)
      ylim <- c(10,80)
      nmsx <- nms_treatment
    }
    
    for (nm_sp in nms_species[-3]) {
    
      niceplot(limX=c(0.75,1.75),limY=ylim, marg=c(1,3,1,1), cexY = 0.7, cexX = .9,
               tick=-0.4,   nlab=4, labX=c(1,1.5), nmlabX=nmsx, nmX = "",
               nmY = nms_vars_toplot[var], lineX = 1,cexXt = 0.8, 
               lasX=1,lasY=1, cexYt = 0.7, lineYt = 2)
      text(1.25, max(ylim) + 0.15*diff(ylim), labels = nms_main[nm_sp], font = 4, cex = 1.2)
      segments(x0 =0.75, y0 = avg_d0[nm_sp, var] %>% unlist(), x1 = 1.75, 
               y1 = avg_d0[nm_sp, var] %>% unlist(),lty = 1, col = "grey")
      segments(x0 =0.75, y0 = avg_d0[nm_sp, var] %>% unlist() + sd_d0[nm_sp, var] %>% unlist(), x1 = 1.75, 
               y1 = avg_d0[nm_sp, var] %>% unlist() + sd_d0[nm_sp, var] %>% unlist(),lty = 3, col = "grey")
      segments(x0 =0.75, y0 = avg_d0[nm_sp, var] %>% unlist() - sd_d0[nm_sp, var] %>% unlist(), x1 = 1.75, 
               y1 = avg_d0[nm_sp, var] %>% unlist() - sd_d0[nm_sp, var] %>% unlist(),lty = 3, col = "grey")
      
      dat_avg <- avgs[avgs$taxonomy == nm_sp,var]
      dat_sd <- sds[sds$taxonomy == nm_sp,var]
      meansey(c(1,1.5), dat_avg %>% t(), dat_sd %>% t(), pchp = 22, 
              colp = cols_border, bgp = cols, cexp = 2, colb = cols_border, 
              lgt = 0.01)
      
      # pvalues
      tmp <- df_anova %>% filter(yvar == var & species == nm_sp)
      text(1.25, max(ylim) - 0.1 * diff(ylim), font = 3,
           labels = paste0("F-value = ", round(tmp$F.value, 1),
                           "   p-value = ", round(tmp["Pr..F."], 3)))  
    }
  }

mtext("Experimental treatment", 1, outer = TRUE, font = 2, line = 2)

dev.off()

}




######################## ESTIMATE NUTRIENT FLUXES ##############################


# --------------------- Fluxes of excreted nutrients ---------------------------

df_excretion <- read.csv(paste0(dir_data, "data_excretion.csv"))


df_excretion <- df_excretion %>% mutate(flux_NH4_min = NH4 / excretion_duration_min,
                                        flux_PO4_min = PO4 / excretion_duration_min,
                                        flux_NH4_hour = flux_NH4_min * 60,
                                        flux_PO4_hour = flux_PO4_min * 60,
                                        flux_NH4_day = flux_NH4_min * 60 * 24,
                                        flux_PO4_day = flux_PO4_min * 60 * 24,
                                        NP = flux_NH4_day / flux_PO4_day) 
avg_fluxes <- df_excretion %>% summarize(avg_flux_NH4_day = mean(flux_NH4_day),
                                         avg_flux_PO4_day = mean(flux_PO4_day),
                                         avg_NP = mean(NP),
                                         sd_flux_NH4_day = sd(flux_NH4_day),
                                         sd_flux_PO4_day = sd(flux_PO4_day),
                                         sd_NP = sd(NP),
                                         avg_body_mass = mean(body_mass, na.rm = TRUE))
# flux per day for 3 fishes
avg_fluxes * 3

# total flux  for 3 fishes for 8 days
avg_fluxes * 3 * 8





# ------------------------ Fluxes of egested nutrients -------------------------


# load gut trait data to estimate the volume of the gut ----
# experimental fishes are 22-39 g and 109-132 mm 
# 122.2+-6.85 mm
# 31.4+-6.4

df_gut_traits <- read.csv(paste0(dir_data, "data_gut_traits.csv"))
gut_vol <- df_gut_traits$VolumeD

# average gut volume (mL ou cm3)
# Volume D: diametre moyen
# VolumeA: aire

mean(gut_vol)
sd(gut_vol)
min(gut_vol) 
max(gut_vol) 

# quantity of feces for 3 fishes 
#density in mg.ml = g.l and dry weight

df_density <- read.csv(paste0(dir_data, "data_feces_density.csv"))

density <- df_density$density

mean(density) 
sd(density)
min(density)
max(density) 


# compute feces quantity in mg

# in one gut
mean(density) * mean(gut_vol) 
min(density) * min(gut_vol) 
max(density) * max(gut_vol) 

# in 3 guts
mean(density) * mean(gut_vol) * 3
min(density) * min(gut_vol) * 3
max(density) * max(gut_vol) * 3

# in 3 guts over the course of the experiment
mean(density) * mean(gut_vol) * 3 * 8
min(density) * min(gut_vol) * 3 * 8
max(density) * max(gut_vol) * 3 * 8




# Load data for nutrient content of feces  ----

df_nfns_exp <- read.csv(paste0(dir_data, "data_nutrients_feces.csv"))

df_nfns_exp %>% summarise(avg_C = mean(pct_C_mean),
                           sd_C = sd(pct_C_mean),
                           avg_N = mean(pct_N_mean),
                           sd_N = sd(pct_N_mean))



