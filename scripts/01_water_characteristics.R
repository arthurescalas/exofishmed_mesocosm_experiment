################################################################################
#  _______  _____  _____ ___ ____  _   _ __  __ _____ ____  
# | ____\ \/ / _ \|  ___|_ _/ ___|| | | |  \/  | ____|  _ \ 
# |  _|  \  / | | | |_   | |\___ \| |_| | |\/| |  _| | | | |
# | |___ /  \ |_| |  _|  | | ___) |  _  | |  | | |___| |_| |
# |_____/_/\_\___/|_|   |___|____/|_| |_|_|  |_|_____|____/ 
#
################################################################################
#
# R script to analyze the non µbiome data describing the characteristics of water
# nutrient, pigments, CMF
#
# Arthur Escalas September 2020
# arthur.escalas@gmail.com
################################################################################


for (keep_F1 in c(TRUE, FALSE)) {

if (keep_F1 == TRUE) {
  dir_save <- paste0(dir_water_characteristics, "with_F1/")
  dir.create(dir_save)
} else {
  dir_save <- paste0(dir_water_characteristics, "without_F1/")
  dir.create(dir_save)  
}


########################### LOAD CLEANED DATA ##################################


df_data <- read.csv(file = paste0(dir_data, "data_water_characteristics.csv"))

nms_vars <- c("NO2_NO3", "NH4", "PO4", "NP", "Chla","total_bact_mL")
nms_vars_toplot <- c("NO2 + NO3", "NH4", "N:P", "Chlorophyll a", 
                     "PO4", "Total bacterial abundance") %>% setNames(nms_vars)

# metadata variables ----

vars_meta <- c("sample_id", "sampling_id", "sampling_day", "sampling_date",
               "sample_type", "tank_id", "treatment", "treatment_tank")

# remove F1 samples if needed ----

if (keep_F1 == FALSE) {
  df_data <- df_data %>% filter(treatment_tank != "F1")
} 



####################### COMPUTE SOME SUMMARY STATS #############################

# aggregate at the tank level ----

df_agg <- df_data %>% group_by(sampling_day, treatment, treatment_tank) %>% 
  summarize_at(nms_vars, mean, na.rm = TRUE) %>% data.frame()

# compute average at D0 for reference
df_agg %>% filter(sampling_day == 0) %>% 
  summarize_at(nms_vars, c(mean, sd), na.rm = TRUE) %>% data.frame()

# average and SD per treatment at D0 and D8 ----

avg_agg <- df_agg %>% filter(sampling_day == 8) %>% group_by(treatment) %>% 
  summarize_at(nms_vars, mean, na.rm = TRUE) %>% data.frame()


out <- data.frame(rbind(avg_agg %>% as.matrix(),
      c("pct_change", apply(avg_agg[,-1], 2, function(x) { (x[2] - x[1]) / x[1] *100 })),
      c("fold_change", apply(avg_agg[,-1], 2, function(x) { x[2] / x[1] }))))

write.csv(out, file = paste0(dir_save, "table_average_fold_change.csv"), row.names = FALSE)


############################## TEST STATS ######################################

ls_data <- list(d0 = df_data %>% filter(sampling_day == 0),
                d8 = df_data %>% filter(sampling_day == 8))

# ---------------------------- MIXED MODELS ------------------------------------

# define the models ----
# mod1 include only treatment effect and does not account for time effect
# model has a random effect of tank within treatment

model_list <- list(mod1 = " ~ treatment + (1| treatment : treatment_tank )")

# fit the lmm ----
lmm_fits <- lapply(ls_data, function(input_data) {
  out <- lapply(nms_vars, function(v) {
    fit_lmm(input_data, model_list, nm_var_Y = v, reml = TRUE)
  }) %>% setNames(nms_vars)
  out
})

saveRDS(lmm_fits, file = paste0(dir_save, "res_lmm.rds"))


# ---------------------- KRUSKALL WALLIS TEST TESTS ----------------------------


df_res_kw <- lapply(ls_data, function(input_data) {
  
  tmp <- input_data %>% group_by(sampling_day, treatment, treatment_tank) %>% 
    summarize_at(nms_vars, mean)
  res <- do.call(rbind, lapply(nms_vars, function(v) {
    out <- do.call(c, kruskal.test(formula(paste0(v, "~ treatment")), data = tmp)[1:3]) %>% 
      data.frame() %>% t()
    row.names(out) <- v
    out
  }) %>% setNames(nms_vars)) %>% data.frame() %>% setNames(c("Chi2", "df","p-value")) %>% 
    rownames_to_column("y_var")
  res
}) %>% reformat_as_df("sampling_day")

write.csv(df_res_kw, file = paste0(dir_save, "table_res_KW.csv"), row.names = FALSE)


# ---------------------------- PLOT THE DATA -----------------------------------


png(paste0(dir_save, "plot_water_variables_d0_vs_d8.png"),
    height = 15, width = 25, res = 400, unit = "cm")
par(mfrow = c(2,3), mar = c(3,3,1,1), mgp = c(2,0.5,0), oma = c(2,2,1,1), las = 1)

for (var in nms_vars) {
  
  xlim <- ylim <- c(0, round(max(setrge(df_data[, var], 10)),2))
  avgs <- df_data %>% group_by(treatment, sampling_day, tank_id, treatment_tank) %>% 
    dplyr::summarize_at(nms_vars, funs(mean(., na.rm = TRUE)))
  sds <- df_data %>% group_by(treatment, sampling_day, tank_id, treatment_tank) %>% 
    dplyr::summarize_at(nms_vars, funs(sd(., na.rm = TRUE)))
  mask <- table(avgs$tank_id) == 1
  avgs <- avgs %>% filter(tank_id != names(mask[mask]))
  sds <- sds %>% filter(tank_id != names(mask[mask]))
  
  pchs <- c(control = 21, fish = 22)
  cols <- c(control = "dodgerblue", fish = "forestgreen")
  plot(NA, main  = nms_vars_toplot[var], xlab = "", 
       ylab = "", xlim = xlim, ylim = ylim)
  abline(a = 0, b = 1, lty = 2)
  for (tk in as.numeric(unique(avgs$tank_id))) {
    mask <-  as.numeric(avgs$tank_id) == tk
    gp <- unique(avgs$treatment[mask]) %>% as.character()
    text(avgs[mask,var] %>% t(), col = "grey",
         labels = unique(avgs$treatment_tank[mask]), adj = c(0.5,-0.5))
    meansexy(avgs[mask,var] %>% t(), sds[mask,var] %>% t(), pchp = pchs[gp], colp = cols[gp], 
             bgp = cols[gp], cexp = 1, colb = cols[gp], lgt = 0.01, lwd_seg = 1)

  }
}
mtext(text = "Value of the variable at day 0", 1, outer = TRUE, font = 2)
mtext(text = "Value of the variable at day 8", 2, outer = TRUE, las = 0, font = 2)
dev.off()



# -------------------------- ONLY VARIABLES AT D8 ------------------------------

# sumarize lmm results ----

df_anova <- lapply(lmm_fits$d8, function(x) {
    x$anova$mod1 %>% data.frame()
}) %>% reformat_as_df(new_var_name = "yvar")


# Prepare data for plot ----

nms_vars_toplot <- c("NO2 + NO3 (µg/L)", "NH4 (µg/L)","PO4 (µg/L)", "N:P",
                     "Chlorophyll a (µg/L)", 
                     "Total bacterial abundance\n(cells/mL, x 1 000 000)") %>% setNames(nms_vars)

df_data_plot <- df_data
df_data_plot$total_bact_mL <- df_data_plot$total_bact_mL / 10^6

avg_d0 <- df_data_plot %>% filter(sampling_day == 0) %>% group_by(tank_id)%>% 
  dplyr::summarize_at(nms_vars, funs(mean(., na.rm = TRUE))) %>% 
  summarize_at(nms_vars, mean, na.rm = TRUE) %>% data.frame()
sd_d0 <- df_data_plot %>% filter(sampling_day == 0) %>% group_by(tank_id)%>% 
  dplyr::summarize_at(nms_vars, funs(mean(., na.rm = TRUE))) %>% 
  summarize_at(nms_vars, sd, na.rm = TRUE) %>% data.frame()

df_plot <- df_data_plot %>% filter(sampling_day == 8)
tmp <- df_plot %>% group_by(treatment, tank_id, treatment_tank) %>% 
  dplyr::summarize_at(nms_vars, funs(mean(., na.rm = TRUE)))
avgs <- tmp %>% group_by(treatment) %>% 
  summarize_at(nms_vars, mean, na.rm = TRUE) %>% data.frame()
sds <- tmp %>% group_by(treatment) %>% 
  summarize_at(nms_vars, sd, na.rm = TRUE) %>% data.frame()

# make the plot ----

png(paste0(dir_save, "plot_water_variables_d8.png"),
    height = 15, width = 15, res = 400, unit = "cm")
par(mfrow = c(3,2), oma = c(4,1,1,0), xpd = TRUE, mgp = c(2.5,0.5,0), mar = c(2,3,1,1))

for (var in nms_vars) {
  ylim <- c(0, round(max(setrge(df_data_plot[, var], 20)),2))
  pchs <- c(control = 21, fish = 22)
  cols <- c(control = "white", fish = "black")
  cols_border <- c(control = "black", fish = "black")
  
  if (var %in% c("Chla", "total_bact_mL")) {
    nmx <- nms_treatment
  } else {
    nmx <-  ""
  }
  niceplot(limX=c(0.75,1.75),limY=ylim, marg=c(1,4,1,1), cexY = 0.7, cexX = .9,
           tick=-0.4,   nlab=4, labX=c(1,1.5), nmlabX=(nmx), nmX = "",
           nmY = nms_vars_toplot[var], lineX = 1,cexXt = 0.8,
           lasX=1,lasY=1, cexYt = 0.7, lineYt = 2)
  segments(x0 =0.75, y0 = avg_d0[var] %>% unlist(), x1 = 1.75, 
           y1 = avg_d0[var] %>% unlist(),lty = 1, col = "grey")
  segments(x0 =0.75, y0 = avg_d0[var] %>% unlist() + sd_d0[var] %>% unlist(), x1 = 1.75, 
           y1 = avg_d0[var] %>% unlist() + sd_d0[var] %>% unlist(),lty = 3, col = "grey")
  segments(x0 =0.75, y0 = avg_d0[var] %>% unlist() - sd_d0[var] %>% unlist(), x1 = 1.75, 
           y1 = avg_d0[var] %>% unlist() - sd_d0[var] %>% unlist(),lty = 3, col = "grey")
  
  tmp <- df_anova[,] %>% filter(yvar == var)
  text(1.25, max(ylim) - 0.1 * diff(ylim), font = 3,
       labels = paste0("F-value = ", round(tmp$F.value, 1),
                       "   p-value = ", round(tmp["Pr..F."], 3)))
  
  # data points
  meansey(c(1,1.5), avgs[,var] %>% t(), sds[,var] %>% t(), pchp = 22, 
             colp = cols_border, bgp = cols, cexp = 1.5, colb = cols_border, 
             lgt = 0.01)
}

mtext("Experimental treatment", 1, outer = TRUE, font = 2, line = 2)

dev.off()


}## eo for keep_F1



#######################  SUMMARIZE RESULTS OF STATS ############################


# -------------------------- KRUSKALL ------------------------------

nms_folds <- c("with_F1","without_F1")

res_ls <- lapply(nms_folds, function(X) {
  a <- read.csv(paste0(dir_water_characteristics, X, "/table_res_KW.csv"))
  names(a) <- paste0(names(a), "_", X)
  a
}) %>% setNames(nms_folds)

tmp <- cbind(res_ls$with_F1, res_ls$without_F1[, c(2,4)]) 
tmp <- tmp[, c(5,1,2,4,6,7)]
tmp[, 3:6] <- apply(tmp[, 3:6], 2, function(x) round(x, 3))

write.csv(tmp, file = paste0(dir_water_characteristics, "table_comparison_KW_F1.csv"),
          row.names = FALSE)


# -------------------------- MIXED MODELS ------------------------------

res_ls <- lapply(nms_folds, function(X) {
  a <- readRDS(paste0(dir_water_characteristics, X, "/res_lmm.rds"))
  a$d8
}) %>% setNames(nms_folds)

# ANOVA on fixed effects ----

df_anova <- lapply(res_ls, function(X) {
  lapply(X, function(x) {
    x$anova$mod1 %>% data.frame()
  }) %>% reformat_as_df(new_var_name = "yvar")
}) %>% reformat_as_df(new_var_name = "F1")

write.csv(df_anova, file = paste0(dir_water_characteristics, 
                                  "table_comparison_lmm_F1.csv"),
          row.names = FALSE)


