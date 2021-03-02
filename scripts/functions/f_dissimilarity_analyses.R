################################################################################
#  _______  _____  _____ ___ ____  _   _ __  __ _____ ____  
# | ____\ \/ / _ \|  ___|_ _/ ___|| | | |  \/  | ____|  _ \ 
# |  _|  \  / | | | |_   | |\___ \| |_| | |\/| |  _| | | | |
# | |___ /  \ |_| |  _|  | | ___) |  _  | |  | | |___| |_| |
# |_____/_/\_\___/|_|   |___|____/|_| |_|_|  |_|_____|____/ 
#
################################################################################
#
# Functions for dissmilarity and beta diversity analyses
#
# Arthur Escalas May 2020
# arthur.escalas@gmail.com
################################################################################

# get alpha and beta diversity
get_hn_div <- function(ps, dir_tree, do_beta = TRUE, do_phylo = FALSE, 
                       which_q = 0:2, path_to_phylog = NULL) {
  
  ps <- microbiome::transform(ps, transform = "compositional")
  
  if (do_phylo) {
    if (is.null(path_to_phylog)) {
      # transform the tree ----
      # 1°/ you need to make a .txt file in a newick format
      filename <- paste0(dir_tree, "phylo_tree.txt")
      write.tree(ps@phy_tree, file = filename)
      # 2°/ you need to read it as a character string
      tree <- paste(readLines(filename))
      # 3°/ newick2phylog{ade4} converts the character string into a phylog object
      tree_phylog <- newick2phylog(tree)
      #saveRDS(tree_phylog, file = paste0(dir_tree, "phylo_tree_phylog.rds"))
    } else {
      tree_phylog <- readRDS(path_to_phylog)
    }
  } else {
    tree_phylog <- NULL
  }
  
  # estimate HN diversity ----
  
  hn_div <- chao_alpha_beta(as(otu_table(ps), "matrix") %>% t(),
                            tree_phylog = tree_phylog, q = which_q,
                            beta = do_beta, run_example = F)
  hn_div
}

get_hill_numbers_dissimilarities <- function(phyloseq_object, dir_tmp) {
  
  # Transform data
  compo <- microbiome::transform(phyloseq_object, "compositional") 

  # -------------------- Estimate the dissmilarities  -------------------------
  
  # prepare a list of trees in phylog format
  filename <- paste0(dir_tmp, "temporary_phylo_tree.txt")
  write.tree(compo@phy_tree, file = filename)
  tree <- paste(readLines(filename))
  tree_phylog <- newick2phylog(tree)
  
  # estimate diversity
  hn_div <- chao_alpha_beta(as(otu_table(compo), "matrix") %>% t(),
                            tree_phylog = tree_phylog, q = c(0,1,2),
                            beta = TRUE, run_example = F)
  # extract only beta div
  diss_hn <- lapply(hn_div[c("beta_taxo", "beta_phylo")], function(X) { 
    lapply(X, as.dist)
  })
  diss_hn <- unlist(diss_hn, recursive = FALSE)
  names(diss_hn) <- gsub("\\.", "_", names(diss_hn))
  
  # -------------------- Return results in a list  -------------------------
  
  ls_diss <- list(taxo_q0 = diss_hn$beta_taxo_q0,
                  taxo_q1 = diss_hn$beta_taxo_q1,
                  taxo_q2 = diss_hn$beta_taxo_q2,
                  phylo_q0 = diss_hn$beta_phylo_q0,
                  phylo_q1 = diss_hn$beta_phylo_q1,
                  phylo_q2 = diss_hn$beta_phylo_q2)
  return(ls_diss)
}





# Function that estimate several dissmilarity indexes on a single phyloseq object

get_all_dissimilarities <- function(phyloseq_object, dir_tmp) {

  # Transform data
  compo <- microbiome::transform(phyloseq_object, "compositional") 
  clr <- microbiome::transform(phyloseq_object, "clr")
  
  # -------------------- Estimate the dissmilarities  -------------------------
  
  #     aitchison ----
  aitchison <-  vegdist(otu_table(clr) %>% t(), "euclidean") 
  
  #     bray-Curtis ----
  bray <-  vegdist(otu_table(compo) %>% t(), "bray") 
  
  #     hill numbers ----
  
  # prepare a list of trees in phylog format
  filename <- paste0(dir_tmp, "temporary_phylo_tree.txt")
  write.tree(compo@phy_tree, file = filename)
  tree <- paste(readLines(filename))
  tree_phylog <- newick2phylog(tree)
  
  # estimate diversity
  hn_div <- chao_alpha_beta(as(otu_table(compo), "matrix") %>% t(),
                            tree_phylog = tree_phylog, q = c(0,1,2),
                            beta = TRUE, run_example = F)
  # extract only beta div
  diss_hn <- lapply(hn_div[c("beta_taxo", "beta_phylo")], function(X) { 
    lapply(X, as.dist)
  })
  diss_hn <- unlist(diss_hn, recursive = FALSE)
  names(diss_hn) <- gsub("\\.", "_", names(diss_hn))
  
  #     unifrac ----
  
  UniFrac <- UniFrac(compo, weighted=FALSE, normalized=TRUE, parallel=FALSE, fast=TRUE)
  wUniFrac <- UniFrac(compo, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
  
  # -------------------- Return results in a list  -------------------------
  
  ls_diss <- list(aitchison = aitchison, braycurtis = bray,
                  taxo_q0 = diss_hn$beta_taxo_q0,
                  taxo_q1 = diss_hn$beta_taxo_q1,
                  taxo_q2 = diss_hn$beta_taxo_q2,
                  phylo_q0 = diss_hn$beta_phylo_q0,
                  phylo_q1 = diss_hn$beta_phylo_q1,
                  phylo_q2 = diss_hn$beta_phylo_q2, 
                  UniFrac  = UniFrac,
                  wUniFrac = wUniFrac)
return(ls_diss)
}


