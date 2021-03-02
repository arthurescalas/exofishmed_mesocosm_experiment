#### LEFSE functions


# Function that run the lefse analysis #########################################

# adapted from https://github.com/xia-lab/MicrobiomeAnalystR

# the function takes a phyloseq oject, a grouping variable, a taxonomic rank and the
# name of the taxa table transformation. 
# phyloseq_obj: your phyloseq object
# variable: the name of the grouping variable to make the comparison
#  (must be contained in the phyloseq_obj)
# taxrank: le niveau taxonomique auquel vous voulez faire l'analyse
# data_transform: method for transforming OTU abundance before analysis (thius
#  rely on the transform method from microbiome package: "clr" or "compositional")
#
# It does the following
# - groups the otu table according to the taxonomic rank provided (taxa are features)
# - tests for overall difference between groups using Kruskall Wallis test, for each features separately
# - tests for pairwise difference between groups using Dunn test, for each features separately
# - estimates an effect size (LDA_score) for each feature using a Linear Discrimiant Analysis (LDA)
# 
# The function then return a list of results:
# res_stats_lda: table of statistical tests and LDA results, taxa in rows and in colmuns:
#     feature: the name of the feature
#     kw_pvalues: pvalue associated with the KW test
#     kw_pvalues_fdr pvalue corrected using False Discovery Rate
#     num_sign_diff: the number of pairwise dunn tests that were significants
#     prop_sign_diff: the proportion of pairwise dunn tests that were signific
#     LDA_score: the effect size estimated using the LDAeffectsize() function
#     min_ab, max_ab, max_diff, the lowest, highest average group abundance and their difference
#     max_gp: name of the group in which the feature is the most abundant
#     res_dunn_test: a table of pvalues for each paiswise Dunn test
#     avg_group_ab: the average abundance of features in each groups
# taxa_level: lke niveau taxonomique
# grouping_factor: la variable de groupe


gimme_lefse <- function(phyloseq_obj, variable, taxrank = "OTU",
                        data_transform = "compositional") {
  
  library(MASS);
  library(microbiome);
  out <- list()
  
  # Prepare the data ==========================================================
  
  # remove unobserved OTUs ----
  vec <<- rowSums(phyloseq_obj@otu_table)
  ps_no0 <- phyloseq::subset_taxa(phyloseq_obj, vec != 0)
  
  # make the variable as a factor
  
  if (is.character(variable) & length(variable) == 1) {
    factor_group <- factor(ps_no0@sam_data[,variable] %>% unlist()) %>% 
      setNames(row.names(ps_no0@sam_data[,variable]))
  } else {
    factor_group <- variable
  }
  factor_levs <- levels(factor_group)
  
  # transform data if needed ---
  # aggregate data at the chosen taxonomic resolution
  # WARNING: it seems that making the analysis at other level than OTU
  # causse an error in the LDA function because of correlated taxa
  
  if (taxrank == "OTU") {
    ps_norm <- microbiome::transform(ps_no0, data_transform)
    data_to_analyze <- ps_norm@otu_table %>% t() %>% data.frame()
  } else {
    # make count table for the corresponding taxonomic rank
    ps_taxrank <- aggregate_taxa(ps_no0, level = taxrank) 
    ps_norm <- microbiome::transform(ps_taxrank, data_transform)
    data_to_analyze <- ps_norm@otu_table %>% t() %>% data.frame()
  }
  
  # Perform the statistical analysis ===========================================
  
  # Kruskal wallis W rank sum test to test difference of abundance ----
  #  of features between samples 
  
  kw_raw_pval <- apply(data_to_analyze, 2, function(x) {
    kruskal.test(x, factor_group)$p.value
  });
  kw_corr_pval <- p.adjust(kw_raw_pval, method = "fdr");
  
  # Dunne test to compare pairs of groups ----
  dunn_raw_pval <- apply(data_to_analyze, 2, function(x) {
    tmp <- dunn.test::dunn.test(x, factor_group)
    res <- tmp$P.adjusted %>% setNames(tmp$comparisons)
    res
  })
  dunn_corr_pval <- apply(dunn_raw_pval, 2, function(x) {
    p.adjust(x, method = "fdr")
  }) %>% t()
  colnames(dunn_corr_pval) <- gsub(" - ", "_vs_", rownames(dunn_raw_pval))
  
  # make a table of statistical test results
  # raw and correcte pvalues of KW test
  # corrected pvalues for each pairwise comparison between groups levels (1 column per pair)
  
  tab_res_stats <- data.frame(kw_pvalues = signif(kw_raw_pval,digits=5),
                              kw_pvalues_fdr = signif(kw_corr_pval,digits=5)) %>% 
    rownames_to_column(var = "feature") %>% 
    mutate(num_sign_diff = apply(dunn_corr_pval, 1, function(x) sum(x < 0.05)),
           prop_sign_diff = num_sign_diff / ncol(dunn_corr_pval))
  
  # reorganize data according to kruskal test pvalue and store them in the output
  
  ord.inx <- order(kw_corr_pval)
  
  data_to_analyze <- data_to_analyze[, ord.inx]
  out$res_stats_lda <- tab_res_stats[ord.inx, ]
  out$res_dunn_test <- dunn_corr_pval[ord.inx, ] %>% data.frame() %>% 
    rownames_to_column(var = "feature")
  
  # Perform the LINEAR DISCRIMINANT ANALYSIS ===================================
  
  set.seed(666);
  
  lda <-  suppressWarnings(lda(factor_group~., data = data_to_analyze, tol=1e-6))
  
  # Make a prediction using the LDA model in order to get sample scoordinates
  plda <- predict(object = lda,
                  newdata = data_to_analyze)
  
  # format the LDA results into a table
  tab_lda_results <- data.frame(t(lda$means)) %>% 
    rownames_to_column(var = "feature") %>% 
    mutate(min_ab = apply(t(lda$means),1,min),
           max_ab = apply(t(lda$means),1,max),
           max_diff = max_ab - min_ab,
           max_gp = factor_levs[apply(t(lda$means),1,which.max)])
    
  
  # Estimate the LDA score
  tab_lda_results$LDA_score <- LDAeffectsize(lda, data_to_analyze, factor_group)
  
  # reorganize columns
  res_LDA <- tab_lda_results %>% 
    dplyr::select(feature, LDA_score, min_ab, max_ab, max_diff, max_gp)
  
  out$res_stats_lda <- out$res_stats_lda %>% left_join(res_LDA, by = "feature")
  
  out$avg_group_ab <- t(lda$means) %>% data.frame() %>% 
    rownames_to_column(var = "feature") 
  
  # Plot some validation plots
  
  # par(mfrow = c(2,2), mar = c(3,3,1,1), mgp = c(1.5,0.5,0))
  # plot(plda$x[,1:2], col = factor_group, main = "Original group")
  # plot(plda$x[,1:2], col = plda$class, main = "LDA classification")
  # plot(tab_lda_results$LDA_score %>% sort, ylab = "LDA score", xlab = "rank")
  # plot(lda$scaling[,1], tab_lda_results$LDA_score, ylab = "LDA score", 
  #      xlab = "original LDA loadings")
  
  # Output the results =========================================================
  
  out$taxa_level <- taxrank;
  out$grouping_factor <- factor_group;
  
  return(out);
}


gimme_lefse_2_groups <- function(phyloseq_obj, variable, taxrank = "OTU",
                                data_transform = "compositional") {
  
  library(MASS);
  library(microbiome);
  out <- list()
  
  # Prepare the data ==========================================================
  
  # remove unobserved OTUs ----
  vec <<- rowSums(phyloseq_obj@otu_table)
  ps_no0 <- subset_taxa(phyloseq_obj, vec != 0)
  
  # make the variable as a factor
  if (is.character(variable) & length(variable) == 1) {
    factor_group <- factor(ps_no0@sam_data[,variable] %>% unlist()) %>% 
      setNames(row.names(ps_no0@sam_data[,variable]))
  } else {
    factor_group <- variable
  }
  factor_levs <- levels(factor_group)
  
  # transform data if needed ---
  # aggregate data at the chosen taxonomic resolution
  # WARNING: it seems that making the analysis at other level than OTU
  # causse an error in the LDA function because of correlated taxa
  
  if (taxrank == "OTU") {
    ps_norm <- microbiome::transform(ps_no0, data_transform)
    data_to_analyze <- ps_norm@otu_table %>% t() %>% data.frame()
  } else {
    # make count table for the corresponding taxonomic rank
    ps_taxrank <- aggregate_taxa(ps_no0, level = taxrank) 
    ps_norm <- microbiome::transform(ps_taxrank, data_transform)
    data_to_analyze <- ps_norm@otu_table %>% t() %>% data.frame()
  }
  
  # Perform the statistical analysis ===========================================
  
  # Wilcoxon rank sum test to test difference of abundance ----
  #  of features between samples 
  
  kw_raw_pval <- apply(data_to_analyze, 2, function(x) {
    wilcox.test(x ~ factor_group)$p.value
  });
  kw_corr_pval <- p.adjust(kw_raw_pval, method = "fdr");
  

  # make a table of statistical test results
  # raw and correcte pvalues of KW test
  # corrected pvalues for each pairwise comparison between groups levels (1 column per pair)
  
  tab_res_stats <- data.frame(kw_pvalues = signif(kw_raw_pval,digits=5),
                              kw_pvalues_fdr = signif(kw_corr_pval,digits=5)) %>% 
    rownames_to_column(var = "feature") 
  
  # reorganize data according to kruskal test pvalue and store them in the output
  
  ord.inx <- order(kw_corr_pval)
  
  data_to_analyze <- data_to_analyze[, ord.inx]
  out$res_stats_lda <- tab_res_stats[ord.inx, ]

  
  # Perform the LINEAR DISCRIMINANT ANALYSIS ===================================
  
  set.seed(666);
  
  lda <-  suppressWarnings(lda(factor_group~., data = data_to_analyze, tol=1e-6))
  
  # Make a prediction using the LDA model in order to get sample scoordinates
  plda <- predict(object = lda,
                  newdata = data_to_analyze)
  
  # format the LDA results into a table
  tab_lda_results <- data.frame(t(lda$means)) %>% 
    rownames_to_column(var = "feature") %>% 
    mutate(min_ab = apply(t(lda$means),1,min),
           max_ab = apply(t(lda$means),1,max),
           max_diff = max_ab - min_ab,
           max_gp = factor_levs[apply(t(lda$means),1,which.max)]) 
  
  # Estimate the LDA score
  tab_lda_results$LDA_score <- LDAeffectsize(lda, data_to_analyze, factor_group)
  
  # reorganize columns
  res_LDA <- tab_lda_results %>% 
    dplyr::select(feature, LDA_score, min_ab, max_ab, max_diff, max_gp)
  
  out$res_stats_lda <- out$res_stats_lda %>% left_join(res_LDA, by = "feature")
  
  out$avg_group_ab <- t(lda$means) %>% data.frame() %>% 
    rownames_to_column(var = "feature") 
  
  # Plot some validation plots
  
  # par(mfrow = c(2,2), mar = c(3,3,1,1), mgp = c(1.5,0.5,0))
  # plot(plda$x[,1:2], col = factor_group, main = "Original group")
  # plot(plda$x[,1:2], col = plda$class, main = "LDA classification")
  # plot(tab_lda_results$LDA_score %>% sort, ylab = "LDA score", xlab = "rank")
  # plot(lda$scaling[,1], tab_lda_results$LDA_score, ylab = "LDA score", 
  #      xlab = "original LDA loadings")
  
  # Output the results =========================================================
  
  out$taxa_level <- taxrank;
  out$grouping_factor <- factor_group;
  
  return(out);
}


# Function that estimate the LDA-based effect size #############################

# Yu lab method for calculation
# https://github.com/YuLab-SMU/MicrobiotaProcess/blob/423177f6ec1cfcd2f5838ae9408ec8c08c5baaf6/R/ml-method.R
# which used some code taken from segata bitbucket
# https://bitbucket.org/nsegata/lefse/src/default/lefse.py

LDAeffectsize <- function(ldares, data_input, factor_group) {
  
  compareclass <- combn(levels(factor_group), 2) %>% t()
  w <- ldares$scaling[,1]
  mm <- ldares$mean
  compareres <- list()
  
  for (p in seq_len(nrow(compareclass))){
    comb <- compareclass[p,]
    w.unit <- w 
    LD <- as.matrix(data_input) %*% as.matrix(w.unit)
    tmpp1 <- data_input[factor_group == comb[1], ]
    tmpp2 <- data_input[factor_group == comb[2], ]
    effect.size <- abs(mean(LD[factor_group == comb[1]]) - mean(LD[factor_group == comb[2]]))
    wfinal <- w.unit * effect.size
    coeff <- abs(wfinal)
    gm <- abs(mm[match(comb[1], rownames(mm)),,drop=FALSE]- mm[match(comb[2], rownames(mm)),,drop=FALSE])
    tmpres <- (gm+coeff) *0.5
    compareres[[p]] <- tmpres
  }
  LDA_score <- apply(do.call(rbind, compareres), 2, function(x) x[which(x == max(x))][1]) # + 1
  LDA_score_scaled <- ((LDA_score - min(LDA_score)) / (max(LDA_score) - min(LDA_score))) * (10^6 - 1)
  LDA_score_final <- log10(LDA_score_scaled)
  LDA_score_final[LDA_score_final == "-Inf"] <- NA
  return(LDA_score_final)
}








