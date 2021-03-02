################################################################################
#  _______  _____  _____ ___ ____  _   _ __  __ _____ ____  
# | ____\ \/ / _ \|  ___|_ _/ ___|| | | |  \/  | ____|  _ \ 
# |  _|  \  / | | | |_   | |\___ \| |_| | |\/| |  _| | | | |
# | |___ /  \ |_| |  _|  | | ___) |  _  | |  | | |___| |_| |
# |_____/_/\_\___/|_|   |___|____/|_| |_|_|  |_|_____|____/ 
#
################################################################################
#
# Utility functions to process sequences using dada2
#
# Arthur Escalas January 2019
# arthur.escalas@gmail.com
################################################################################


### 01_select_samples.R ########################################################

# retrieve the path of fastq files corresponding to a set of samples ----
# need a vector containing the name of the runs and a vector of sample names

get_path_to_fastq <- function(nm_runs, nm_sples, path_to_seq_data) {
  
  lapply(as.character(nm_runs), function(run) {
    dir_run <- path_to_seq_data[[run]]
    path_to_seq_run <- list.files(dir_run, full.names = TRUE)
    idx <- lapply(unique(nm_sples) %>% 
                    as.character(), function(x) {
                      grep(x, path_to_seq_run)
                    }) %>% unlist()
    paths_fastq <- path_to_seq_run[idx]
    return(paths_fastq)
  }) %>% unlist()
}

# copy the fastq files, save them in a directory and unzip them ----

copy_and_save_fastq <- function(file_paths, dir_output) {
  
  # remove the files with 'extendedFrags" in the name
  mask <- grep(pattern = "extendedFrags", file_paths)
  if (length(mask) != 0) {
    file_paths <- file_paths[- mask]
  }
  
  # Remove the files from the sample folder
  for (f in list.files(dir_output, full.names = TRUE)) {
    file.remove(f)
  }
  
  # grab the files of interest in the source data folder and copy them in the project data folder
  file.copy(file_paths, dir_output, overwrite = TRUE)
  
  # unzip the fastq files in the data folder of the project
  lapply(list.files(dir_output), function(file) {
    gunzip(paste0(dir_output, file))
  })
}

### 02_process_sequences.R #####################################################

# make the plot of the rrror rate, save it and print it

plot_error_rates <- function(err_rate, r1r2, dir_save) {
  
  plot_err <- plotErrors(err_rate, nominalQ = TRUE)
  ggsave(plot = plot_err, device = "png", width = 15, height = 15, 
         filename = paste0(dir_save, "plot_error_rate_", r1r2, ".png"), 
         scale = 1, units = "cm")
  print(plot_err)
}


# transform_and_save_ASV_sequences

transform_and_save_ASV_sequences <- function(seqtab, dir_save, nm_save) {
  
  # Clean the sequence table
  seqtab_trans <- as.data.frame(t(seqtab)) %>% 
    rownames_to_column(var = "sequence") %>% 
    rowid_to_column(var = "OTUNumber") %>% 
    mutate(OTUNumber = sprintf("otu_%04d", OTUNumber)) %>% 
    mutate(sequence = str_replace_all(sequence, "(-|\\.)", ""))
  
  # Extract the ASVs sequences from the sequence table
  seq_out <- Biostrings::DNAStringSet(seqtab_trans$sequence)
  names(seq_out) <- seqtab_trans$OTUNumber
  
  # Save the sequences
  Biostrings::writeXStringSet(seq_out, 
                              str_c(dir_save, nm_save), 
                              compress = FALSE, width = 20000)
  
}


### 04_data_cleaning.R #####################################################


create_dir <- function(dir_parent, date, nm_obj) {
  
  path <- paste0(dir_parent, date, "/")
  assign(nm_obj, path, envir = .GlobalEnv)
  dir.create(path)
}

## @knitr compute_taxa_prevalence

# Compute prevalence of each feature and add taxonomy and total read counts

get_prevalence_tab <- function(ps) {
  
  prev <- apply(otu_table(ps),
                ifelse(taxa_are_rows(ps), yes = 1, no = 2),
                function(x){ sum(x > 0) })
  prevdf <- data.frame(Prevalence = prev,
                       TotalAbundance = taxa_sums(ps),
                       ps %>% tax_table() %>% as("matrix") %>% 
                         data.frame())
  return(prevdf)
}

compute_prevalence_stats <- function(prevdf, taxon_level = "Phylum") {
  
  out <- lapply(split(prevdf, prevdf[,taxon_level]), function(x) {
    data.frame(num_asv  = nrow(x),
               avg_prev = round(mean(x$Prevalence), 1),
               sd_prev  = round(sd(x$Prevalence), 1),
               min_prev = min(x$Prevalence),
               max_prev = max(x$Prevalence),
               sum_prev = sum(x$Prevalence))
  }) %>% bind_rows()
  out[[taxon_level]] <- prevdf[,taxon_level] %>% levels()
  out <- out %>% arrange(avg_prev) %>% dplyr::select(taxon_level, everything())
  return(out)
}

get_effect_of_prev_filtering <- function(ps, prevdf) {
  
  num_sples   <- nsamples(ps)
  out <- lapply(seq(0, 0.2, 0.01), function(x) {
    num_occ <- round(x * num_sples)
    tmp <- prevdf[prevdf$Prevalence >= num_occ,] %>% dim()
    return(data.frame(prevalence_thd = x*100, num_occ_thd = num_occ,
                      num_remaining_asv = tmp[1], 
                      prop_remaining_asv = round(tmp[1] / nrow(prevdf) * 100, 2)))
  }) %>% bind_rows()
  return(out)
}


filter_asv_low_prevalence <- function(ps, thd_prev, prevdf) {
  
  num_sples  <- nsamples(ps)
  
  # Define prevalence threshold according to proportion of samples
  prevalenceThreshold <- round(thd_prev * num_sples)
  
  # list taxa to keep and remove them from the summary prevalence table
  taxa_to_keep <- rownames(prevdf)[prevdf$Prevalence >= prevalenceThreshold]
  
  # prune taxa from the phyloseq object
  out <- prune_taxa(taxa_to_keep, ps)
  return(out)
}


plot_asv_and_otu_numbers <- function(ps) {
  
  asv <- data.frame(sorted = 1:ntaxa(ps), nreads = sort(taxa_sums(ps), TRUE))
  sple <- data.frame(sorted = 1:nsamples(ps), nreads = sort(sample_sums(ps), TRUE))
  
  par(mfrow = c(2,2), mar = c(3,3,1,1), oma = c(1,1,1,1), mgp = c(2,0.5,0))
  
  plot(asv, cex= 0.2, pch = 21, xlab = "ordered ASVs", ylab = "number of reads")
  plot(sple, cex= 0.2, pch = 21, xlab = "ordered samples", ylab = "number of reads")
  plot(asv$sorted, log(asv$nreads), cex= 0.2, pch = 21, xlab = "ordered ASVs", 
       ylab = "log(number of reads)")
  plot(sple$sorted, log(sple$nreads), cex= 0.2, pch = 21, xlab = "ordered samples", 
       ylab = "log(number of reads)")
}


show_num_seq_per_sample <- function(ps) {
  
  seq_num <- colSums(ps@otu_table@.Data)
  
  hist(seq_num %>% sort(), 20, main = "Number of reads per sample", 
       xlab = "Number of sequences", ylim = c(0,50))
  text(max(seq_num) - 0.5 * (max(seq_num) - min(seq_num)), 50,
       labels = paste0(names(summary(seq_num)), collapse = "\t"))
  text(max(seq_num) - 0.5 * (max(seq_num) - min(seq_num)), 45,
       labels = paste0(round(summary(seq_num)), collapse = "\t"))
}



check_remaining_samples_after_rarefaction <- function(ps) {
  
  # Check the effect of various rarefaction thrshold on the number of remaining samples
  
  tmp <- ps@sam_data@.Data %>% as.data.frame()
  seq_num <- readcount(ps)
  
  tab_rar_thd_effect <- lapply(seq(1000, 20000, 1000), function(x) {
    res <- tmp[seq_num >= x,] %>% nrow()
    return(data.frame(raref_thd = x, num_sple = res, 
                      prop_sple = round(res/nrow(tmp) * 100, 0)))
  }) %>% bind_rows()
  
  knitr::kable(tab_rar_thd_effect, 
               caption = "Number of remaning samples at various rarefaction thresholds") %>%
    kable_styling(bootstrap_options = "striped", full_width = F)
  
}


save_sequence_table <- function(ps, taxonomy_table, dir_output) {
  
  nms_otu <- row.names(ps@otu_table)
  otu_seq <- taxonomy_table %>%
    distinct(OTUNumber, sequence) %>% # OTUNumber == name of the MOTU
    mutate(sequence = toupper(sequence)) %>% 
    filter(OTUNumber %in% nms_otu)
  
  fa <- character(2 * nrow(otu_seq))
  fa[c(TRUE, FALSE)] = sprintf(">%s", otu_seq$OTUNumber)
  fa[c(FALSE, TRUE)] = as.character(otu_seq$sequence)
  
  dir.create(dir_output)
  writeLines(fa, paste(dir_output, "OTU_sequences.fasta", sep=""))
  
}

save_list_of_sequence_table_from_phyloseq <- function(ps, taxonomy_table, dir_output) {
  
  nms_otu <- row.names(ps@otu_table)
  num_otu <- length(nms_otu)
  num_chunks <- round(num_otu / 2000)
  otu_seq <- taxonomy_table %>%
    distinct(OTUNumber, sequence) %>% # OTUNumber == name of the MOTU
    mutate(sequence = toupper(sequence)) %>% 
    filter(OTUNumber %in% nms_otu)
  
  seq(1, num_chunks * 2000, length.out = num_chunks)  
  counter <- 1
  
  for (i in 1:num_chunks) {
    
    if (i != num_chunks) {
      which_seq <- counter:(i * 2000)
    } else {
      which_seq <- counter:num_otu
    }
    counter <- max(which_seq) + 1
    
    tmp_tab <- otu_seq[which_seq, ]
    
    fa <- character(2 * nrow(tmp_tab))
    fa[c(TRUE, FALSE)] = sprintf(">%s", tmp_tab$OTUNumber)
    fa[c(FALSE, TRUE)] = as.character(tmp_tab$sequence)
    
    dir.create(dir_output)
    writeLines(fa, paste(dir_output, "OTU_sequences_chunk_", i, ".fasta", sep=""))
  }
}


save_list_of_sequence_table <- function(sequence_table, dir_output) {
  
  nms_otu <- sequence_table[, "OTUNumber"] %>% unlist() %>% as.vector()
  num_otu <- length(nms_otu)
  num_chunks <- round(num_otu / 2000)
  otu_seq <- sequence_table %>%
    distinct(OTUNumber, sequence) %>% # OTUNumber == name of the MOTU
    mutate(sequence = toupper(sequence)) %>% 
    filter(OTUNumber %in% nms_otu)
  
  seq(1, num_chunks * 2000, length.out = num_chunks)  
  counter <- 1
  
  for (i in 1:num_chunks) {
    
    if (i != num_chunks) {
      which_seq <- counter:(i * 2000)
    } else {
      which_seq <- counter:num_otu
    }
    counter <- max(which_seq) + 1
    
    tmp_tab <- otu_seq[which_seq, ]
    
    fa <- character(2 * nrow(tmp_tab))
    fa[c(TRUE, FALSE)] = sprintf(">%s", tmp_tab$OTUNumber)
    fa[c(FALSE, TRUE)] = as.character(tmp_tab$sequence)
    
    dir.create(dir_output)
    writeLines(fa, paste(dir_output, "OTU_sequences_chunk_", i, ".fasta", sep=""))
  }
}








