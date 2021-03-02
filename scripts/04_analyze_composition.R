################################################################################
#          _______  _____  _____ ___ ____  _   _ __  __ _____ ____  
#         | ____\ \/ / _ \|  ___|_ _/ ___|| | | |  \/  | ____|  _ \ 
#         |  _|  \  / | | | |_   | |\___ \| |_| | |\/| |  _| | | | |
#         | |___ /  \ |_| |  _|  | | ___) |  _  | |  | | |___| |_| |
#         |_____/_/\_\___/|_|   |___|____/|_| |_|_|  |_|_____|____/ 
#
################################################################################
#
# R script to analyze composition of communities
#
# Arthur Escalas October 2020
# arthur.escalas@gmail.com
################################################################################


################################ LOAD DATA #####################################

for (keep_F1 in c(TRUE, FALSE)) {
  
  if (keep_F1 == TRUE) {
    today <- "with_F1"
    dir_save <- paste0(dir_composition, today, "/")
    dir.create(dir_save)
  } else {
    today <- "without_F1"
    dir_save <- paste0(dir_composition, today, "/")
    dir.create(dir_save)
  }
  

  # names of ranks and indices ----

nms_ranks_lefse <- c("Phylum", "Class", "Order", "Family", "Genus", "Species", "ASV")
nms_ranks_compo <- c("Phylum", "Class", "Family", "ASV")
nms_ranks_div   <- c("Phylum", "Family", "ASV")

nms_index <- c("taxo_q0", "taxo_q1", "phylo_q0", "phylo_q1")
names(nms_index) <- nms_index


# Load the phyloseq objects and the diversity estimates ----
ps_data <- readRDS(paste0(dir_data, "phyloseq_object_core_with_diversity.rds"))

# Load diversity estimates for all the ranks, taxo and phylo, alpha and beta
ls_hn_div <- readRDS(paste0(dir_data, "list_all_diversity_estimates_per_ranks.rds"))

# load the list of phyloseq objects at different taxomnomic ranks
ls_ps_data <- readRDS(paste0(dir_data, "list_glom_phyloseq_per_taxa_rank.rds"))


# remove the F1 samples if needed ----

if (keep_F1 == TRUE) {
  
  # identify the F1 samples
  metadata <- meta(ps_data)
  nms_f1_sples <- metadata %>% filter(treatment_tank == "F1") %>%
    dplyr::select(sample_id_fastq) %>% unlist()
  
  # remove these sampls from the diversity objects ----
  ls_hn_div <- lapply(ls_hn_div, function(a) {
    a$alpha_taxo <- a$alpha_taxo[! row.names(a$alpha_taxo) %in% nms_f1_sples,]
    a$beta_taxo <- lapply(a$beta_taxo, function(x) {
      x[! row.names(x) %in% nms_f1_sples, ! colnames(x) %in% nms_f1_sples]
    })
    a$alpha_phylo <- a$alpha_phylo[! row.names(a$alpha_phylo) %in% nms_f1_sples,]
    a$beta_phylo <- lapply(a$beta_phylo, function(x) {
      x[! row.names(x) %in% nms_f1_sples, ! colnames(x) %in% nms_f1_sples]
    })
    a
  })
  
  # list of phyloseq objects agglometrated at different taxonomic levels ----
  ls_ps_data <- lapply(ls_ps_data, function(x) {
    row.names(x@tax_table) <- taxa_names(x)
    subset_samples(x, ! sample_id_fastq %in% nms_f1_sples)
  })
  
  # phyloseq object with alpha taxonoic diversity at the OTU level
  ps_data <- subset_samples(ps_data, ! sample_id_fastq %in% nms_f1_sples)
} 


# extract data from phyloseq ----

# metadata
metadata <- meta(ps_data) %>% data.frame()
nms_vars <- names(metadata)[which(! names(metadata) %in% nms_index)]
for (i in nms_vars) {
  metadata[, i] <- factor(metadata[, i], levels = sort(unique(metadata[, i])))
}
#taxonomy
taxonomy <- ps_data@tax_table@.Data %>% data.frame()
ls_taxonomy <- lapply(ls_ps_data, tax_table)

# extract dissimilarity matrices ----

ls_diss <- lapply(ls_hn_div, function(rk) {
  out <- flatten(lapply(c("taxo", "phylo"), function(X) {
    dat <- rk[[paste0("beta_", X)]]
    # lapply(dat, as.dist)
  })) %>% setNames(nms_index)
  out
})

# transfom them into data.frame ----

ls_diss_df <- lapply(ls_diss, function(rk) {
  out <- lapply(rk, function(X) {
    dendextend::dist_long(as.dist(X))
  })
  out
})


############################# COMPARE SAMPLE TYPES #############################

# Make a directory to store the results ----

dir_out <- paste0(dir_save, "sample_types/")
dir.create(dir_out)

# create some objects for analyses and remove unecessary samples ----

metadata_tmp <- metadata
fact <- factor(metadata$sample_type)

ls_ps_data_tmp <- ls_ps_data
ps_data_tmp <- ps_data
ls_diss_tmp <- lapply(ls_diss, function(rk) { 
  lapply(rk, function(d) {
    as.dist(d)
  })
})

#============================= PCOA ORDINATION =================================

# Make the PCOA ----

ls_pcoa <- lapply(ls_diss_tmp, function(rk) { 
  lapply(rk, function(d) pcoa(D = as.dist(d)) )
})

saveRDS(ls_pcoa, file = paste0(dir_out, "list_pcoa.rds"))


#============================= COMPOSITION =================================

dir_plots <- paste0(dir_out, "tables_composition/")
dir.create(dir_plots)

# Transform data for plot: merge samples and tidyfy ----

ls_data_plot <- lapply(nms_ranks_compo, function(rk) {
    merge_samples(ls_ps_data_tmp[[rk]], "sample_type", sum) %>% 
    microbiome::transform("compositional")
}) %>% setNames(nms_ranks_compo)

# Save the composition tables -----

  lapply(nms_ranks_compo, function(rk) {
    
    # Save the table of composition ----
    dat <- ls_data_plot[[rk]]
    colnames(dat@tax_table) <- c("Domain", "Phylum", "Class",  "Order",
                                 "Family", "Genus", "Species",  "ASV", "other")
    otu <- dat %>% otu_table() %>% as.matrix() %>% t()
    tax <- dat %>% tax_table()
    out <- cbind(tax, round(otu * 100, 2))
    write.csv(cbind(tax, otu), row.names = FALSE,
              file = paste0(dir_plots, "table_composition_", rk, ".csv"))
  })


# =============== TEST DIFFERENCE BETWEEN GROUPS ===============================

# Run the tests ----

ls_res <- lapply(nms_ranks_div, function(rk) {
  lapply(ls_diss_tmp[[rk]], function(d) {
    tmp <- WdS.test(d, fact, nrep = 999, strata = NULL)
    out <- data.frame(p_value = tmp$p.value, statistic = tmp$statistic)
    return(round(out, 3))
  }) %>% reformat_as_df(new_var_name = "data_type")
}) %>% setNames(nms_ranks_div)

saveRDS(ls_res, file = paste0(dir_out, "res_wdstest.rds"))


#     format the results as a table ----

df_res <- reformat_as_df(ls_res, new_var_name = "rank")
df_res$q <- factor(str_split_fixed(df_res$data_type, "_", 2)[,2])
df_res$div_facet <- factor(str_split_fixed(df_res$data_type, "_", 2)[,1])

write.csv(df_res, file = paste0(dir_out, "table_results_wdstest.csv"), row.names =FALSE)


#============================= LEFSE ANALYSIS =================================

dir_lefse <- paste0(dir_out, "lefse/")
dir.create(dir_lefse)

# Run the LEFSE analysis ----

ls_res_lefse <- list()
for (rk in nms_ranks_lefse) {
  ls_res_lefse[[rk]] <- gimme_lefse(ps_data_tmp, variable = "sample_type", 
                                    taxrank = rk, 
                                    data_transform = "compositional")
}
saveRDS(ls_res_lefse, file = paste0(dir_lefse, "ls_res_lefse.rds"))


# Export the results ----

lapply(nms_ranks_lefse, function(rk) {
  
  # taxonomy ----
  
  tab_taxa <- ps_data@tax_table@.Data %>% data.frame()
  if (rk != "ASV") {
    tab_taxa <- do.call(rbind, lapply(split(tab_taxa, tab_taxa[, rk]), function(x) {
      x[1,]
    }))
    tab_taxa$feature <- tab_taxa[, rk]
  } else {
    tab_taxa$feature <- tab_taxa[, "ASV"]
  }
  
  # merge taxonomy and lefse results ----
  
  tmp <- ls_res_lefse[[rk]]$res_stats_lda
  tmp$feature <- gsub("\\.", "-", tmp$feature)
  out <- tmp %>%  
    filter(num_sign_diff != 0 & kw_pvalues < 0.05) %>% 
    left_join(tab_taxa, by = "feature") %>% 
    dplyr::select(tab_taxa %>% names(), "max_gp", everything()) %>% 
    arrange(Domain, Phylum, Class, Order, Family, Genus, Species, OTU)
  
  write.csv(out, file = paste0(dir_lefse, "table_lefse_", rk, ".csv"), row.names = FALSE)
})




############################# COMPARE WATER SAMPLES ###########################

# Make a directory to store the results ----

dir_out <- paste0(dir_save, "water_samples/")
dir.create(dir_out)

# create some objects for analyses and remove unecessary samples ----

mask <- metadata$sample_type == "Water" & metadata$sampling_day == 8  
metadata_tmp <- metadata[mask,]
fact <- factor(metadata_tmp$treatment, levels = sort(unique(metadata_tmp$treatment)))

ps_data_tmp <- subset_samples(ps_data, mask)

ls_ps_data_tmp <- lapply(ls_ps_data, function(X) {
  subset_samples(X, mask)
})

ls_diss_tmp <- lapply(ls_diss, function(rk) { 
  lapply(rk, function(d) {
    as.dist(d[mask, mask])
  })
})

#============================= PCOA ORDINATION =================================

# Make the PCOA ----

ls_pcoa <- lapply(ls_diss_tmp, function(rk) { 
  lapply(rk, function(d) {
    pcoa(D = d) 
  })
})

saveRDS(ls_pcoa, file = paste0(dir_out, "list_pcoa.rds"))


#============================= COMPOSITION =================================

dir_plots <- paste0(dir_out, "tables_composition/")
dir.create(dir_plots)

# Transform data for plot: merge samples and tidyfy ----

ls_data_plot <- lapply(nms_ranks_compo, function(rk) {
  merge_samples(ls_ps_data_tmp[[rk]], "treatment", sum) %>% 
    microbiome::transform("compositional")
}) %>% setNames(nms_ranks_compo)

# Make and save the stacked barplots -----

lapply(nms_ranks_compo, function(rk) {
  
  # Save the table of composition ----
  dat <- ls_data_plot[[rk]]     
  colnames(dat@tax_table) <- c("Domain", "Phylum", "Class",  "Order",  
                               "Family", "Genus", "Species",  "ASV", "other")
  otu <- dat %>% otu_table() %>% as.matrix() %>% t()
  tax <- dat %>% tax_table()
  out <- cbind(tax, round(otu * 100, 2))
  write.csv(cbind(tax, otu), row.names = FALSE,
            file = paste0(dir_plots, "table_composition_", rk, ".csv"))
})


# =============== TEST DIFFERENCE BETWEEN GROUPS ===============================

ls_res <- list()

for (rk in nms_ranks_div) {
  for (idx in names(ls_diss_tmp[[rk]])) {
    d <- ls_diss_tmp[[rk]][[idx]]
    tmp <- adonis(d ~ treatment / treatment_tank, metadata_tmp, 
                  strata = metadata_tmp$treatment)
    out <- data.frame(tmp$aov.tab)[, c(3,5,6)] %>% 
      rownames_to_column("factor") %>% 
      filter(factor != "Residuals")
    out$factor <- gsub("metadata_tmp\\$", "", as.character(out$factor))
    names(out) <- c("factor","F-value","R2", "p-value")
    ls_res[[rk]][[idx]] <- out
  }
}

saveRDS(ls_res, file = paste0(dir_out, "res_permanova.rds"))

#     format the results as a table ----

df_res <- lapply(ls_res, reformat_as_df, new_var_name = "data_type") %>% 
  reformat_as_df(new_var_name = "rank")
df_res$q <- factor(str_split_fixed(df_res$data_type, "_", 2)[,2])
df_res$div_facet <- factor(str_split_fixed(df_res$data_type, "_", 2)[,1])

write.csv(df_res, file = paste0(dir_out, "table_results_permanova.csv"), row.names =FALSE)


#============================= LEFSE ANALYSIS =================================

dir_lefse <- paste0(dir_out, "lefse/")
dir.create(dir_lefse)

# Run the LEFSE analysis ----

ls_res_lefse <- list()
for (rk in nms_ranks_lefse) {
  ls_res_lefse[[rk]] <- gimme_lefse_2_groups(ps_data_tmp, 
                                             variable = fact, 
                                             taxrank = rk, 
                                             data_transform = "compositional")
}
saveRDS(ls_res_lefse, file = paste0(dir_lefse, "ls_res_lefse.rds"))


# Export the results ----

lapply(nms_ranks_lefse, function(rk) {
  
  # taxonomy ----
  tab_taxa <- ps_data_tmp@tax_table@.Data %>% data.frame()
  if (rk != "ASV") {
    tab_taxa <- do.call(rbind, lapply(split(tab_taxa, tab_taxa[, rk]), function(x) {
      x[1,]
    }))
    tab_taxa$feature <- as.character(tab_taxa[, rk])
  } else {
    tab_taxa$feature <- as.character(tab_taxa[, "ASV"])
  }

  # merge taxonomy and lefse results ----
  tmp <- ls_res_lefse[[rk]]$res_stats_lda
  tmp$feature <- gsub("\\.", "-", tmp$feature)
  out <- tmp %>%  
    filter(kw_pvalues < 0.05) %>% 
    left_join(tab_taxa, by = "feature") %>% 
    dplyr::select(tab_taxa %>% names(), "max_gp", everything()) %>% 
    arrange(Domain, Phylum, Class, Order, Family, Genus, Species, OTU)
  
  write.csv(out, file = paste0(dir_lefse, "table_lefse_", rk, ".csv"), row.names = FALSE)
})

       
       
############################ COMPARE SEDIMENT SAMPLES ##########################

# Make a directory to store the results ----

dir_out <- paste0(dir_save, "sediment_samples/")
dir.create(dir_out)

# create some objects for analyses and remove unecessary samples ----

mask <- metadata$sample_type == "Sediment" & 
  metadata$sampling_day == 8 &
  metadata$treatment != "Environmental_control" &
  metadata$position_in_the_tank %in% c("tank")
metadata_tmp <- metadata[mask,]
fact <- factor(metadata_tmp$treatment, levels = sort(unique(metadata_tmp$treatment)))

ps_data_tmp <- subset_samples(ps_data, mask)

ls_ps_data_tmp <- lapply(ls_ps_data, function(X) {
  subset_samples(X, mask)
})

ls_diss_tmp <- lapply(ls_diss, function(rk) { 
  lapply(rk, function(d) {
    as.dist(d[mask, mask])
  })
})

#============================= PCOA ORDINATION =================================

# Make the PCOA ----

ls_pcoa <- lapply(ls_diss_tmp, function(rk) { 
  lapply(rk, function(d) {
    pcoa(D = d) 
  })
})

saveRDS(ls_pcoa, file = paste0(dir_out, "list_pcoa.rds"))


#============================= COMPOSITION =================================

dir_plots <- paste0(dir_out, "tables_composition/")
dir.create(dir_plots)

# Transform data for plot: merge samples and tidyfy ----

ls_data_plot <- lapply(ls_ps_data_tmp, function(rk) {
  merge_samples(rk, "treatment", sum) %>% 
    microbiome::transform("compositional")
})

# Make and save the stacked barplots -----

lapply(nms_ranks_compo, function(rk) {
  
  # Save the table of composition ----
  dat <- ls_data_plot[[rk]]     
  colnames(dat@tax_table) <- c("Domain", "Phylum", "Class",  "Order",  
                               "Family", "Genus", "Species",  "ASV", "other")
  otu <- dat %>% otu_table() %>% as.matrix() %>% t()
  tax <- dat %>% tax_table()
  out <- cbind(tax, round(otu * 100, 2))
  write.csv(cbind(tax, otu), row.names = FALSE,
            file = paste0(dir_plots, "table_composition_", rk, ".csv"))

})


# =============== TEST DIFFERENCE BETWEEN GROUPS ===============================

ls_res <- list()

for (rk in nms_ranks_div) {
  for (idx in names(ls_diss_tmp[[rk]])) {
    d <- ls_diss_tmp[[rk]][[idx]]
    tmp <- adonis(d ~ treatment / treatment_tank, metadata_tmp, 
                  strata = metadata_tmp$treatment)
    out <- data.frame(tmp$aov.tab)[, c(3,5,6)] %>% 
      rownames_to_column("factor") %>% 
      filter(factor != "Residuals")
    out$factor <- gsub("metadata_tmp\\$", "", as.character(out$factor))
    names(out) <- c("factor","F-value","R2", "p-value")
    ls_res[[rk]][[idx]] <- out
  }
}

saveRDS(ls_res, file = paste0(dir_out, "res_permanova.rds"))

#     format the results as a table ----

df_res <- lapply(ls_res, reformat_as_df, new_var_name = "data_type") %>% 
  reformat_as_df(new_var_name = "rank")
df_res$q <- factor(str_split_fixed(df_res$data_type, "_", 2)[,2])
df_res$div_facet <- factor(str_split_fixed(df_res$data_type, "_", 2)[,1])

write.csv(df_res, file = paste0(dir_out, "table_results_permanova.csv"), row.names =FALSE)


#============================= LEFSE ANALYSIS =================================

dir_lefse <- paste0(dir_out, "lefse/")
dir.create(dir_lefse)

# Run the LEFSE analysis ----

ls_res_lefse <- list()
for (rk in nms_ranks_lefse) {
  ls_res_lefse[[rk]] <- gimme_lefse_2_groups(ps_data_tmp,
                                             variable = fact, 
                                             taxrank = rk, 
                                             data_transform = "compositional")
}
saveRDS(ls_res_lefse, file = paste0(dir_lefse, "ls_res_lefse.rds"))


# Export the results ----

lapply(nms_ranks_lefse, function(rk) {
  
  # taxonomy ----
  
  tab_taxa <- ps_data_tmp@tax_table@.Data %>% data.frame()
  if (rk != "ASV") {
    tab_taxa <- do.call(rbind, lapply(split(tab_taxa, tab_taxa[, rk]), function(x) {
      x[1,]
    }))
    tab_taxa$feature <- as.character(tab_taxa[, rk])
  } else {
    tab_taxa$feature <- as.character(tab_taxa[, "ASV"])
  }
  
  # merge taxonomy and lefse results ----
  
  tmp <- ls_res_lefse[[rk]]$res_stats_lda
  tmp$feature <- gsub("\\.", "-", tmp$feature)
  out <- tmp %>%  
    filter(kw_pvalues < 0.05) %>% 
    left_join(tab_taxa, by = "feature") %>% 
    dplyr::select(tab_taxa %>% names(), "max_gp", everything()) %>% 
    arrange(Domain, Phylum, Class, Order, Family, Genus, Species, OTU)
  
  write.csv(out, file = paste0(dir_lefse, "table_lefse_", rk, ".csv"), row.names = FALSE)
  
})      





#################### COMPARE SEDIMENT CLOSED SHELTER  ##########################

# Make a directory to store the results ----

dir_out <- paste0(dir_save, "sediment_closed_shelter/")
dir.create(dir_out)

# create some objects for analyses and remove unecessary samples ----

mask <- metadata$sample_type == "Sediment" & 
  metadata$sampling_day == 8 &
  metadata$treatment != "Environmental_control" &
  metadata$position_in_the_tank %in% c("closed_shelter")
metadata_tmp <- metadata[mask,]
fact <- factor(metadata_tmp$treatment, levels = sort(unique(metadata_tmp$treatment)))

ps_data_tmp <- subset_samples(ps_data, mask)

ls_ps_data_tmp <- lapply(ls_ps_data, function(X) {
  X <- subset_samples(X, mask)
  X@sam_data$treatment_position <- fact
  X
})

ls_diss_tmp <- lapply(ls_diss, function(rk) { 
  lapply(rk, function(d) {
    as.dist(d[mask, mask])
  })
})

#============================= PCOA ORDINATION =================================

# Make the PCOA ----

ls_pcoa <- lapply(ls_diss_tmp, function(rk) { 
  lapply(rk, function(d) {
    pcoa(D = d) 
  })
})

saveRDS(ls_pcoa, file = paste0(dir_out, "list_pcoa.rds"))


#============================= COMPOSITION =================================

dir_plots <- paste0(dir_out, "tables_composition/")
dir.create(dir_plots)

# Transform data for plot: merge samples and tidyfy ----

ls_data_plot <- lapply(nms_ranks_compo, function(rk) {
  merge_samples(ls_ps_data_tmp[[rk]], "treatment_position", sum) %>% 
    microbiome::transform("compositional")
}) %>% setNames(nms_ranks_compo)

# Make and save the stacked barplots -----

lapply(nms_ranks_compo, function(rk) {
  
  # Save the table of composition ----
  dat <- ls_data_plot[[rk]]     
  colnames(dat@tax_table) <- c("Domain", "Phylum", "Class",  "Order",  
                               "Family", "Genus", "Species",  "ASV", "other")
  otu <- dat %>% otu_table() %>% as.matrix() %>% t()
  tax <- dat %>% tax_table()
  out <- cbind(tax, round(otu * 100, 2))
  write.csv(cbind(tax, otu), row.names = FALSE,
            file = paste0(dir_plots, "table_composition_", rk, ".csv"))
  
  
})


# =============== TEST DIFFERENCE BETWEEN GROUPS ===============================

ls_res <- list()

for (rk in nms_ranks_div) {
  for (idx in names(ls_diss_tmp[[rk]])) {
    d <- ls_diss_tmp[[rk]][[idx]]
    tmp <- adonis(d ~ treatment, metadata_tmp)
    out <- data.frame(tmp$aov.tab)[, c(3,5,6)] %>% 
      rownames_to_column("factor")
    out$factor <- gsub("metadata_tmp\\$", "", as.character(out$factor))
    names(out) <- c("factor","F-value","R2", "p-value")
    ls_res[[rk]][[idx]] <- out
  }
}

saveRDS(ls_res, file = paste0(dir_out, "res_permanova_treatment_treatmenttank.rds"))

#     format the results as a table ----

df_res <- lapply(ls_res, reformat_as_df, new_var_name = "data_type") %>% 
  reformat_as_df(new_var_name = "rank")
df_res$q <- factor(str_split_fixed(df_res$data_type, "_", 2)[,2])
df_res$div_facet <- factor(str_split_fixed(df_res$data_type, "_", 2)[,1])

write.csv(df_res, file = paste0(dir_out, "table_results_permanova.csv"), row.names =FALSE)





####################### COMPARE SEDIMENT OPEN SHELTERS #########################


# Make a directory to store the results ----

dir_out <- paste0(dir_save, "sediment_open_shelter/")
dir.create(dir_out)

# create some objects for analyses and remove unecessary samples ----

mask <- metadata$sample_type == "Sediment" & 
  metadata$sampling_day == 8 &
  metadata$treatment != "Environmental_control" &
  metadata$position_in_the_tank %in% c("open_shelter")
metadata_tmp <- metadata[mask,]
fact <- factor(metadata_tmp$treatment, levels = sort(unique(metadata_tmp$treatment)))

ps_data_tmp <- subset_samples(ps_data, mask)

ls_ps_data_tmp <- lapply(ls_ps_data, function(X) {
  X <- subset_samples(X, mask)
  X@sam_data$treatment_position <- fact
  X
})

ls_diss_tmp <- lapply(ls_diss, function(rk) { 
  lapply(rk, function(d) {
    as.dist(d[mask, mask])
  })
})

#============================= PCOA ORDINATION =================================

# Make the PCOA ----

ls_pcoa <- lapply(ls_diss_tmp, function(rk) { 
  lapply(rk, function(d) {
    pcoa(D = d) 
  })
})

saveRDS(ls_pcoa, file = paste0(dir_out, "list_pcoa.rds"))


#============================= COMPOSITION =================================

dir_plots <- paste0(dir_out, "tables_composition/")
dir.create(dir_plots)

# Transform data for plot: merge samples and tidyfy ----

ls_data_plot <- lapply(ls_ps_data_tmp, function(rk) {
  merge_samples(rk, "treatment_position", sum) %>% 
    microbiome::transform("compositional")
})

# Make and save the stacked barplots -----

lapply(nms_ranks_compo, function(rk) {
  
  # Save the table of composition ----
  dat <- ls_data_plot[[rk]]     
  colnames(dat@tax_table) <- c("Domain", "Phylum", "Class",  "Order",  
                               "Family", "Genus", "Species",  "ASV", "other")
  otu <- dat %>% otu_table() %>% as.matrix() %>% t()
  tax <- dat %>% tax_table()
  out <- cbind(tax, round(otu * 100, 2))
  write.csv(cbind(tax, otu), row.names = FALSE,
            file = paste0(dir_plots, "table_composition_", rk, ".csv"))
  
})


# =============== TEST DIFFERENCE BETWEEN GROUPS ===============================

ls_res <- list()

for (rk in nms_ranks_div) {
  for (idx in names(ls_diss_tmp[[rk]])) {
    d <- ls_diss_tmp[[rk]][[idx]]
    tmp <- adonis(d ~ treatment / treatment_tank, metadata_tmp, 
                  strata = metadata_tmp$treatment)
    out <- data.frame(tmp$aov.tab)[, c(3,5,6)] %>% 
      rownames_to_column("factor")
    out$factor <- gsub("metadata_tmp\\$", "", as.character(out$factor))
    names(out) <- c("factor","F-value","R2", "p-value")
    ls_res[[rk]][[idx]] <- out
  }
}

saveRDS(ls_res, file = paste0(dir_out, "res_permanova_treatment_treatmenttank.rds"))

#     format the results as a table ----

df_res <- lapply(ls_res, reformat_as_df, new_var_name = "data_type") %>% 
  reformat_as_df(new_var_name = "rank")
df_res$q <- factor(str_split_fixed(df_res$data_type, "_", 2)[,2])
df_res$div_facet <- factor(str_split_fixed(df_res$data_type, "_", 2)[,1])

write.csv(df_res, file = paste0(dir_out, "table_results_permanova.csv"), row.names =FALSE)


#============================= LEFSE ANALYSIS =================================

dir_lefse <- paste0(dir_out, "lefse/")
dir.create(dir_lefse)

# Run the LEFSE analysis ----

ls_res_lefse <- list()
for (rk in nms_ranks_lefse) {
  ls_res_lefse[[rk]] <- gimme_lefse_2_groups(ps_data_tmp,
                                    variable = fact, 
                                    taxrank = rk, 
                                    data_transform = "compositional")
}
saveRDS(ls_res_lefse, file = paste0(dir_lefse, "ls_res_lefse.rds"))


# Export the results ----

lapply(nms_ranks_lefse, function(rk) {
  
  # taxonomy ----
  
  tab_taxa <- ps_data_tmp@tax_table@.Data %>% data.frame()
  if (rk != "ASV") {
    tab_taxa <- do.call(rbind, lapply(split(tab_taxa, tab_taxa[, rk]), function(x) {
      x[1,]
    }))
    tab_taxa$feature <- as.character(tab_taxa[, rk])
  } else {
    tab_taxa$feature <- as.character(tab_taxa[, "ASV"])
  }
  
  # merge taxonomy and lefse results ----
  
  tmp <- ls_res_lefse[[rk]]$res_stats_lda
  tmp$feature <- gsub("\\.", "-", tmp$feature)
  out <- tmp %>%  
    filter(kw_pvalues < 0.05) %>% 
    left_join(tab_taxa, by = "feature") %>% 
    dplyr::select(tab_taxa %>% names(), "max_gp", everything()) %>% 
    arrange(Domain, Phylum, Class, Order, Family, Genus, Species, OTU)
  
  write.csv(out, file = paste0(dir_lefse, "table_lefse_", rk, ".csv"), row.names = FALSE)
})    




############################ COMPARE TURF SAMPLES #############################

# Make a directory to store the results ----

dir_out <- paste0(dir_save, "turf_samples/")
dir.create(dir_out)

# create some objects for analyses and remove unecessary samples ----
mask <- metadata$sample_type == "Turf" & 
  metadata$sampling_day == 8 &
  metadata$treatment != "Environmental_control"
metadata_tmp <- metadata[mask,]
fact <- factor(metadata_tmp$treatment, levels = sort(unique(metadata_tmp$treatment)))

ps_data_tmp <- subset_samples(ps_data, mask)

ls_ps_data_tmp <- lapply(ls_ps_data, function(X) {
  subset_samples(X, mask)
})

ls_diss_tmp <- lapply(ls_diss, function(rk) { 
  lapply(rk, function(d) {
    as.dist(d[mask, mask])
  })
})

#============================= PCOA ORDINATION =================================

# Make the PCOA ----

ls_pcoa <- lapply(ls_diss_tmp, function(rk) { 
  lapply(rk, function(d) {
    pcoa(D = d) 
  })
})

saveRDS(ls_pcoa, file = paste0(dir_out, "list_pcoa.rds"))


#============================= COMPOSITION =================================

dir_plots <- paste0(dir_out, "tables_composition/")
dir.create(dir_plots)

# Transform data for plot: merge samples and tidyfy ----

ls_data_plot <- lapply(ls_ps_data_tmp, function(rk) {
  merge_samples(rk, "treatment", sum) %>% 
    microbiome::transform("compositional")
})

# Make and save the stacked barplots -----

lapply(nms_ranks_compo, function(rk) {
  
  # Save the table of composition ----
  dat <- ls_data_plot[[rk]]     
  colnames(dat@tax_table) <- c("Domain", "Phylum", "Class",  "Order",  
                               "Family", "Genus", "Species",  "ASV", "other")
  otu <- dat %>% otu_table() %>% as.matrix() %>% t()
  tax <- dat %>% tax_table()
  out <- cbind(tax, round(otu * 100, 2))
  write.csv(cbind(tax, otu), row.names = FALSE,
            file = paste0(dir_plots, "table_composition_", rk, ".csv"))

})


# =============== TEST DIFFERENCE BETWEEN GROUPS ===============================

ls_res <- list()

for (rk in nms_ranks_div) {
  for (idx in names(ls_diss_tmp[[rk]])) {
    d <- ls_diss_tmp[[rk]][[idx]]
    tmp <- adonis(d ~ treatment / treatment_tank, metadata_tmp, 
                  strata = metadata_tmp$treatment)
    out <- data.frame(tmp$aov.tab)[, c(3,5,6)] %>% 
      rownames_to_column("factor") %>% 
      filter(factor != "Residuals")
    out$factor <- gsub("metadata_tmp\\$", "", as.character(out$factor))
    names(out) <- c("factor","F-value","R2", "p-value")
    ls_res[[rk]][[idx]] <- out
  }
}

saveRDS(ls_res, file = paste0(dir_out, "res_permanova.rds"))
ls_res <- readRDS(paste0(dir_out, "res_permanova.rds"))

#     format the results as a table ----

df_res <- lapply(ls_res, reformat_as_df, new_var_name = "data_type") %>% 
  reformat_as_df(new_var_name = "rank")
df_res$q <- factor(str_split_fixed(df_res$data_type, "_", 2)[,2])
df_res$div_facet <- factor(str_split_fixed(df_res$data_type, "_", 2)[,1])

write.csv(df_res, file = paste0(dir_out, "table_results_permanova.csv"), row.names =FALSE)


#============================= LEFSE ANALYSIS =================================

dir_lefse <- paste0(dir_out, "lefse/")
dir.create(dir_lefse)

# Run the LEFSE analysis ----

ls_res_lefse <- list()
for (rk in nms_ranks_lefse) {
  ls_res_lefse[[rk]] <- gimme_lefse_2_groups(ps_data_tmp,
                                             variable = fact, 
                                             taxrank = rk, 
                                             data_transform = "compositional")
}
saveRDS(ls_res_lefse, file = paste0(dir_lefse, "ls_res_lefse.rds"))


# Export the results ----

lapply(nms_ranks_lefse, function(rk) {
  
  # taxonomy ----
  
  tab_taxa <- ps_data_tmp@tax_table@.Data %>% data.frame()
  if (rk != "ASV") {
    tab_taxa <- do.call(rbind, lapply(split(tab_taxa, tab_taxa[, rk]), function(x) {
      x[1,]
    }))
    tab_taxa$feature <- as.character(tab_taxa[, rk])
  } else {
    tab_taxa$feature <- as.character(tab_taxa[, "ASV"])
  }
  
  # merge taxonomy and lefse results ----
  
  tmp <- ls_res_lefse[[rk]]$res_stats_lda
  tmp$feature <- gsub("\\.", "-", tmp$feature)
  out <- tmp %>%  
    filter(kw_pvalues < 0.05) %>% 
    left_join(tab_taxa, by = "feature") %>% 
    dplyr::select(tab_taxa %>% names(), "max_gp", everything()) %>% 
    arrange(Domain, Phylum, Class, Order, Family, Genus, Species, OTU)
  
  write.csv(out, file = paste0(dir_lefse, "table_lefse_", rk, ".csv"), row.names = FALSE)
  
})      




####################### SUMMARIZE THE PERMANOVA RESULTS ########################

# count the umber of times treatment and treatment_tank effects are significants

nms_folds <- list.files(dir_save)

ls_permanova <- list()
for (f in nms_folds) {
  file <- list.files(paste0(dir_save, f, "/"), pattern = "_permanova.csv")
  if (length(file) != 0) {
    if (file.exists(paste0(dir_save, f, "/",file))) {
      ls_permanova[[f]] <- read.csv(paste0(dir_save, f, "/",file))
    }
  }
}

df_summary_stats <- ls_permanova %>% reformat_as_df("compartment") %>% 
  dplyr::filter(!factor %in% c("Residuals","Total")) %>% 
  dplyr::filter(rank %in% nms_ranks_div & q %in% c("q0","q1")) %>% 
  arrange(factor, compartment, div_facet, rank, q)

df_summary_signif <- df_summary_stats %>% 
    group_by(compartment, factor) %>% 
    summarise(num_signif = sum(p.value < 0.05),
              avg_p = mean(p.value),
              avg_f = mean(F.value),
              avg_r2 = mean(R2)) %>% data.frame()

write.csv(df_summary_stats, row.names = FALSE,
          file = paste0(dir_save, "table_summary_results_permanova.csv"))
write.csv(df_summary_signif, row.names = FALSE,
          file = paste0(dir_save, "table_summary_results_significant_permanova.csv"))


####################### MAKE THE FINAL PCOA PLOT ###############################


# Load the PCOAs ----

nms_folds <- c("water_samples","sediment_open_shelter")

ls_pcoa <- list()
for (f in nms_folds) {
  file <- paste0(dir_save, f, "/list_pcoa.rds")
  ls_pcoa[[f]] <- readRDS(file)
}

# laod the lefse results at the OTU level ----

ls_lefse <- list()
for (f in nms_folds) {
  file <- paste0(dir_save, f, "/lefse/table_lefse_ASV.csv")
  ls_lefse[[f]] <- read.csv(file)
}




# Make the plots ----

dir_plots <- paste0(dir_save, "plots_pcoa/")
dir.create((dir_plots))

for (rk in c("Phylum","Class","Family","Genus")) {

  for (idx in c("taxo_q0","taxo_q1","phylo_q0","phylo_q1")) {

  png(paste0(dir_plots, "figure_PCOA_", rk, "_", idx, ".png"),
      height = 18, width = 30, unit = "cm", res = 400)
  
  layout(matrix(c(1,1,2,2,3,3,3,
                  4,4,5,5,6,6,6
                  ),2,7,TRUE))
  par(las = 1, mar = c(6,4,2,1), oma = c(2,2,2,1), mgp = c(2.5,.5,0), xpd = NA,
      font.lab = 2)
  
    # -----------------  water samples ------------------------
  
    # objects for plot 
  
    mask <- metadata$sample_type == "Water" & metadata$sampling_day == 8  
    metadata_tmp <- metadata[mask,]
    fact <- factor(metadata_tmp$treatment, levels = sort(unique(metadata_tmp$treatment)))
  
    col_border <- cols_sample_type["Water"]
    pchvec <- pch_treatment[as.character(metadata_tmp$treatment)]
    colvec <- as.character(metadata_tmp$treatment)
    colvec[as.character(metadata_tmp$treatment) == "Control"] <- "white"
    colvec <- gsub("Fish", col_border, colvec)
  
    # plot the PCOA 
    
    # dat <- ls_pcoa$water_samples[[rk]][[idx]]
    dat <- ls_pcoa$water_samples$ASV[[idx]]
    eig <- dat$values$Eigenvalues[dat$values$Eigenvalues > 0]
    expl_var <- round(eig / sum(eig) * 100, 1)
    lapply(list(c(1,2)), function(axes) {
      plot(dat$vectors[, axes], pch = pchvec, bg = colvec, col = col_border, 
           main = "Water", cex = 2, cex.lab = 1.2,cex.main = 1.5,
           xlab = paste0("Axis ", axes[1], " (", expl_var[axes[1]], "%)"),
           ylab = paste0("Axis ", axes[2], " (", expl_var[axes[2]], "%)"))
    })  
    lapply(unique(metadata_tmp$treatment_tank), function(x) {
      tmp <- dat$vectors[metadata_tmp$treatment_tank == x,1:2]
      s <- 1:(nrow(tmp)-1)
      segments(tmp[s,1], tmp[s,2], tmp[s+1,1], tmp[s+1,2], col = col_border)
      segments(tmp[nrow(tmp),1], tmp[nrow(tmp),2], tmp[1,1], tmp[1,2], col = col_border)
      
    })
    points(dat$vectors[, 1:2], pch = pchvec, bg = colvec, col = col_border, 
           cex = 2)
  
    # plot the abundance of families differntiating treatments by lefse 
    
    dat <- ls_lefse$water_samples %>% filter(max_gp == "Fish")
    mask_unknown <- grep(pattern = "Unknown", dat[, rk])
    if (length(mask_unknown) != 0) {
      dat <- dat[-mask_unknown, ]
    }
    
    # number of OTUs per taxa
    dat_plot_num <- do.call(rbind, lapply(split(dat, dat[, rk], drop = TRUE), function(X) {
      nrow(X) 
    })) %>% t()
    ylim <- c(0, round(max(setrge(dat_plot_num))))
    p <- barplot(dat_plot_num, col = paste0(col_border, 70),cex.lab = 1.2,
                 border = c("black"), ylab = "Number of ASVs", xaxt="n")
    labs <- colnames(dat_plot_num)
    text(cex=1, x= seq(1.2, by = 1.2, length.out = length(labs)) - 0.6 ,
         y= 0 - 0.05 * diff(ylim), labs, srt=45, adj = c(1,1))
    
    # abundance of these OTUs

    tmp <- ls_ps_data$ASV@otu_table %>% as.matrix()
    
    ps_data <- readRDS(paste0(dir_data, 
                              "phyloseq_object_cleaned_with_tree_rarefied.rds")) 
    toto <- subset_taxa(ps_data, ASV %in% row.names(tmp))  
    toto <- subset_samples(toto, sample_id_fastq %in% colnames(tmp))  
    otu <- toto@otu_table %>% as.matrix()
    
    dat_plot_ab <- do.call(rbind, lapply(split(dat, dat[, rk], drop = TRUE), function(X) {
      
      mask <- row.names(otu) %in% X$ASV
      do.call(c, tapply(colSums(otu[mask, metadata_tmp$sample_id_fastq]),
                        metadata_tmp$treatment,
                        function(x) {
                          c(avg = mean(x), sd = sd(x))
                        }))
    })) %>% t()
    
    toplot <- log(dat_plot_ab+1)
    toplot[toplot == "-Inf"] <- 0
    ylim <- c(0, round(max(setrge(toplot)),2))
    
    p <- barplot(toplot[c(1,3),], beside = TRUE, col = c("white",col_border),
                 cex.lab = 1.2, ylim = c(0,8),
                 border = c("black","black"), 
                 ylab = "log (number of sequences of ASVs)", xaxt="n")

    
    labs <- sapply(labs, function(x) c(x, "")) %>% as.vector()
    text(cex=1, x= p +0.6,
         y= 0 - 0.05 * diff(ylim), labs, srt=45, adj = c(1,1))

    # -----------------  Open shelter samples ------------------------
    
    # objects for plot 
    
    mask <- metadata$sample_type == "Sediment" & 
      metadata$sampling_day == 8 &
      metadata$treatment != "Environmental_control" &
      metadata$position_in_the_tank %in% c("open_shelter")
    metadata_tmp <- metadata[mask,]
    fact <- factor(metadata_tmp$treatment, levels = sort(unique(metadata_tmp$treatment)))
  
    col_border <- cols_sample_type["Sediment"]
    pchvec <- pch_treatment[as.character(metadata_tmp$treatment)]
    colvec <- as.character(metadata_tmp$treatment)
    colvec[as.character(metadata_tmp$treatment) == "Control"] <- "white"
    colvec <- gsub("Fish", col_border, colvec)
    
    # plot the PCOA 
    
    dat <- ls_pcoa$sediment_open_shelter[[rk]][[idx]]
    dat <- ls_pcoa$sediment_open_shelter$ASV[[idx]]
    eig <- dat$values$Eigenvalues[dat$values$Eigenvalues > 0]
    expl_var <- round(eig / sum(eig) * 100, 1)
    lapply(list(c(1,2)), function(axes) {
      plot(dat$vectors[, axes], pch = pchvec, bg = colvec, col = col_border, 
           main = "Sediment open shelter", cex = 2, cex.lab = 1.2,cex.main = 1.5,
           xlab = paste0("Axis ", axes[1], " (", expl_var[axes[1]], "%)"),
           ylab = paste0("Axis ", axes[2], " (", expl_var[axes[2]], "%)"))
    })  
    lapply(unique(metadata_tmp$treatment_tank), function(x) {
      tmp <- dat$vectors[metadata_tmp$treatment_tank == x,1:2]
      s <- 1:(nrow(tmp)-1)
      segments(tmp[s,1], tmp[s,2], tmp[s+1,1], tmp[s+1,2], col = col_border)
      segments(tmp[nrow(tmp),1], tmp[nrow(tmp),2], tmp[1,1], tmp[1,2], col = col_border)
      
    })
    points(dat$vectors[, 1:2], pch = pchvec, bg = colvec, col = col_border, 
           cex = 2)
    
    # plot the abundance of families differntiating treatments by lefse 
    
    dat <- ls_lefse$sediment_open_shelter %>% filter(max_gp == "Fish")
    mask_unknown <- grep(pattern = "Unknown", dat[, rk])
    if (length(mask_unknown) != 0) {
      dat <- dat[-mask_unknown, ]
    }
    
    # number of OTUs per taxa
    dat_plot_num <- do.call(rbind, lapply(split(dat, dat[, rk], drop = TRUE), function(X) {
      nrow(X) 
    })) %>% t()
    ylim <- c(0, round(max(setrge(dat_plot_num))))
    p <- barplot(dat_plot_num, col = paste0(col_border, 70),cex.lab = 1.2,
                 border = c("black"), ylab = "Number of ASVs", xaxt="n")
    labs <- colnames(dat_plot_num)
    text(cex=1, x= seq(1.2, by = 1.2, length.out = length(labs)) - 0.6 ,
         y= 0 - 0.05 * diff(ylim), labs, srt=45, adj = c(1,1))
    
    # abundance of these OTUs
    
    tmp <- ls_ps_data$ASV@otu_table %>% as.matrix()
    
    toto <- subset_taxa(ps_data, ASV %in% row.names(tmp))  
    toto <- subset_samples(toto, sample_id_fastq %in% colnames(tmp))  
    otu <- toto@otu_table %>% as.matrix()
    
    dat_plot_ab <- do.call(rbind, lapply(split(dat, dat[, rk], drop = TRUE), function(X) {
      
      mask <- row.names(otu) %in% X$ASV
      do.call(c, tapply(colSums(otu[mask, metadata_tmp$sample_id_fastq]),
                        metadata_tmp$treatment,
                        function(x) {
                          c(avg = mean(x), sd = sd(x))
                        }))
    })) %>% t()
    
    toplot <- log(dat_plot_ab+1)
    toplot[toplot == "-Inf"] <- 0
    ylim <- c(0, round(max(setrge(toplot)),2))
    
    p <- barplot(toplot[c(1,3),], beside = TRUE, col = c("white",col_border),
                 cex.lab = 1.2, ylim = c(0,8),
                 border = c("black","black"), 
                 ylab = "log (number of sequences of ASVs)", xaxt="n")

    labs <- sapply(labs, function(x) c(x, "")) %>% as.vector()
    text(cex=1, x= p +0.6,
         y= 0 - 0.05 * diff(ylim), labs, srt=45, adj = c(1,1))
   
    # ----
  dev.off()

  }
}

}
