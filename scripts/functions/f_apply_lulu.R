apply_lulu <- function(taxonomy_table, sample_names, path_lulu, match_lulu = 84, co_occurence_lulu = 0.95, match_list = NULL) {
 
  # Clean the taxonomy_table
  taxonomy_table <- taxonomy_table %>% dplyr::select( - grep("_boot", names(taxonomy_table)))
  
  # Create necessary files
  otu_seq <- taxonomy_table %>%
    distinct(OTUNumber, sequence) %>% # OTUNumber == name of the MOTU
    mutate(sequence = toupper(sequence))
 
  if (is.null(match_list)) {
    # Convert to fasta
    fa <- character(2 * nrow(otu_seq))
    fa[c(TRUE, FALSE)] = sprintf(">%s", otu_seq$OTUNumber)
    fa[c(FALSE, TRUE)] = as.character(otu_seq$sequence)
    writeLines(fa, paste(path_lulu, "OTU_sequences.fasta", sep=""))
   
    # Launch blastn with bash
    system(paste("makeblastdb -in ", path_lulu, "OTU_sequences.fasta -parse_seqids -dbtype nucl", sep=""))
   
    # Matchlist with blast
    # We can change -task to 1) "blastn" or 2) "megablast"
    system(paste("blastn -db ", path_lulu, 
                 "OTU_sequences.fasta -outfmt '6 qseqid sseqid pident' -out  ", 
                 path_lulu, "match_list_blastn.txt -task blastn -qcov_hsp_perc 60 -perc_identity 60 -query ", 
                 path_lulu, "OTU_sequences.fasta", sep=""))
    # Import matchlist
    match_list <- read.table(paste(path_lulu, "match_list_blastn.txt", sep=""), 
                            header=FALSE,as.is=TRUE, stringsAsFactors=FALSE)
  }

  # ------------------------------ #
  # Analysis
 
  # Format for LULU function
  otutab <- taxonomy_table %>% dplyr::select(OTUNumber, sample_names) %>% 
    as.data.frame()
  otutab <- column_to_rownames(otutab, var = "OTUNumber")

  # Put the otu name in row name

  # Run LULU
  res_lulu <- lulu(otutab, match_list, minimum_match = match_lulu, 
                   minimum_relative_cooccurence = co_occurence_lulu) # 94% is 3 mismatchs
 
  # Clean the logs
  system("rm lulu.log*")
 
  # Output
  return(res_lulu)
}
