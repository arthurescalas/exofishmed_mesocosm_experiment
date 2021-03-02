################################################################################
#  _______  _____  _____ ___ ____  _   _ __  __ _____ ____  
# | ____\ \/ / _ \|  ___|_ _/ ___|| | | |  \/  | ____|  _ \ 
# |  _|  \  / | | | |_   | |\___ \| |_| | |\/| |  _| | | | |
# | |___ /  \ |_| |  _|  | | ___) |  _  | |  | | |___| |_| |
# |_____/_/\_\___/|_|   |___|____/|_| |_|_|  |_|_____|____/ 
#
################################################################################
#
# R script to globally set up project wp4_proj_001_bioinfo_pipeline
#
# Arthur Escalas June 2019
# arthur.escalas@gmail.com
################################################################################


# ==============================================================================
# cleaning memory
# ==============================================================================

rm(list = ls()) 

# ==============================================================================
# Define the directories of the project
# create objects corresponding to the directories and create directories locally
# ==============================================================================

# ------------------------------------------------------------------------------
# identify the project name, find where we are and create local directories

dir_project <- getwd()
proj_name   <- unlist(strsplit(dir_project, 
                               split = 'marbec_exofishmed', 
                               fixed = TRUE))[2]
dir_project_global <- unlist(strsplit(dir_project, 
                                      split = proj_name, 
                                      fixed = TRUE))[1]

dir_scripts   <- paste0(dir_project, "/scripts/")
dir_functions <- paste0(dir_scripts, "/functions/")

dir_data     <- paste0(dir_project, "/data/")

dir_analyses <- paste0(dir_project, "/analyses/")
dir_water_characteristics    <- paste0(dir_analyses, "01_water_characteristics/")
dir_nutrients_dynamics       <- paste0(dir_analyses, "02_nutrients_dynamics/")
dir_alpha_div   <- paste0(dir_analyses, "03_alpha_diversity/")
dir_composition <- paste0(dir_analyses, "04_analyze_composition/")


dir.create(dir_scripts)

dir.create(dir_analyses)

dir.create(dir_water_characteristics)
dir.create(dir_nutrients_dynamics)
dir.create(dir_alpha_div)
dir.create(dir_composition)


# ------------------------------------------------------------------------------
# Soucre material from the utility_function project

# source utility functions
for (f in list.files(dir_functions, full.names = T)) { source(f) }

# source package list
source(paste0(dir_scripts, "packages.R"))
install_if_not_there(cran_packages, type = "CRAN")
install_if_not_there(bioc_packages, type = "bioconductor")
sapply(c(cran_packages, bioc_packages), require, character.only = TRUE)




# ------------------------------------------------------------------------------
# Objects for analyses

df_col_pch <- data.frame(compartment = c("Water","Sediment","Turf", "tank",
                                         "closed_shelter", "open_shelter"),
                         pch_ctrl  = c(1,1,1,1,2,0),
                         pch_fish  = c(21,21,21,21,24,22),
                         color = c("dodgerblue", "darkorange4","firebrick1",
                                   "darkorange4", "darkorange4", "darkorange4")) %>% 
  column_to_rownames("compartment")


nms_treatment <- c(Control = "Control", Fish = "Fish")
nms_sample_type <- c("Water","Sediment","Turf","Algae", "Fish")

cols_treatment <- c(Environmental_control = "grey", Fish = "black", 
                    Control = "dodgerblue")
cols_sample_type <- c("#1E90FF","#8B4500","#FF3030","#228B22", "#000000") %>% 
  setNames(nms_sample_type)


pch_treatment <- c(Control = 21, Fish = 21, Environmental_control = 3)
pch_tank_position <- c(open_shelter = 22, closed_shelter = 24, 
                       tank = 21)

