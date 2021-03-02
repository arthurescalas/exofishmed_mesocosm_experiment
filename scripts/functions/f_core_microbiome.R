################################################################################
#  _______  _____  _____ ___ ____  _   _ __  __ _____ ____  
# | ____\ \/ / _ \|  ___|_ _/ ___|| | | |  \/  | ____|  _ \ 
# |  _|  \  / | | | |_   | |\___ \| |_| | |\/| |  _| | | | |
# | |___ /  \ |_| |  _|  | | ___) |  _  | |  | | |___| |_| |
# |_____/_/\_\___/|_|   |___|____/|_| |_|_|  |_|_____|____/ 
#
################################################################################
#
# Functions to make core microbiome analysis
#
# Arthur Escalas November 2019
# arthur.escalas@gmail.com
################################################################################



# Define the function that do the core analysis
get_core_taxa <- function(ps) {
  
  # remove the otus not detected in the dataset
  # ps <- subset_taxa(ps, rowSums(microbiome::abundances(ps)) != 0)
  
  # Abundance vs occurence 
  otu_tab       <- microbiome::abundances(ps) %>% t() 
  otu_aboccur   <- abuocc(otu_tab, panel = FALSE)
  otu_occurence <- otu_aboccur$plt.spc
  otu_avg_ab    <- colSums(otu_tab) / otu_occurence
  
  # transform spc.plt vector into table in order to calculate specific richness
  otu_richness <- data.frame(otu_aboccur$spc.plt)
  
  # Variance calculation
  otu_square    <- otu_tab^2
  otu_sumsquare <- data.frame(colSums(otu_square))
  otu_variance  <- otu_sumsquare / otu_occurence - otu_avg_ab^2
  
  # Estimate dispersion for each OTU
  otu_dispersion <- (otu_variance / otu_avg_ab) * otu_occurence
  
  # Calculate IC for Poisson distribution using Chi square distribution 
  #   (value and formula within Zar p574)
  otu_poisson_ic <- pois.exact(otu_occurence, conf.level = 0.95)
  
  # Make the final table describing OTU distribution
  otu_stats <- data.frame(otu_name   = colnames(otu_tab),
                          occurrence = otu_occurence,
                          avg_ab     = otu_avg_ab,
                          dispersion = otu_dispersion) %>% 
    bind_cols(otu_poisson_ic[, c("lower","upper","conf.level")]) %>% 
    setNames(c("otu_name", "occurrence", "avg_ab", "dispersion", "lower", 
               "upper","conf.level")) %>% 
    mutate(core   = dispersion > upper)#,
    #        Phylum = tax_table(ps)[otu_name %in% tax_table(ps)[,"OTU"], "Phylum"],
    #        Class  = tax_table(ps)[otu_name %in% tax_table(ps)[,"OTU"], "Class"],
    #        Order  = tax_table(ps)[otu_name %in% tax_table(ps)[,"OTU"], "Order"])
  
  # output
  return(otu_stats)
}

# And the function to plots it


draw_core_plot <- function(core_tab, title) {
  
  core_core_tab <- core_tab %>% filter(core == TRUE)
  plot(core_tab$occurrence, log(core_tab$dispersion), pch = 21, col = "grey", bg = "grey",
       cex = 0.5, ylab = "Dispersion index", xlab = "Occurrence", main = title)
  points(core_core_tab$occurrence, log(core_core_tab$dispersion), 
         col = "forestgreen", bg = "forestgreen", cex = 0.5, pch = 21)
  points(core_tab$occurrence, log(core_tab$lower), cex = 0.3, col = "blue", 
         bg = "blue", pch = 21)
  points(core_tab$occurrence, log(core_tab$upper), cex = 0.3, col = "red", 
         bg = "red", pch = 21)
  
  plot(core_tab$occurrence, log(core_tab$avg_ab), cex = 0.5, col = "grey", pch = 21, 
       bg = "grey", ylab = "Average abundance", xlab = "Occurrence", main = title)
  points(core_core_tab$occurrence, log(core_core_tab$avg_ab), 
         col = "forestgreen", bg = "forestgreen", cex = 0.5, pch = 21)
}


# Ratio of sequences from the Core / total_pool ----

get_core_ratio <- function(ps, ps_core) {
  ratio_core_global <- data.frame(prop_seq = round(sum(readcount(ps_core)) / 
                                                     sum(readcount(ps)) *100, 1),
                                  prop_otu = round(ntaxa(ps_core) / 
                                                     ntaxa(ps) * 100, 1),
                                  num_seq = sum(num_seq = readcount(ps_core)),
                                  num_otu = ntaxa(ps_core))
  ratio_core_sple <- data.frame(prop_seq = round(readcount(ps_core) / 
                                                   sample_sums(ps) *100, 1),
                                prop_otu = round(microbiome::alpha(ps_core, 
                                                                   index = "observed") %>% unlist() / 
                                                   ntaxa(ps) * 100, 1),
                                num_seq = readcount(ps_core),
                                num_otu = microbiome::alpha(ps_core, index = "observed") %>% unlist()
  )
  return(list(global = ratio_core_global, per_sample = ratio_core_sple))
}
