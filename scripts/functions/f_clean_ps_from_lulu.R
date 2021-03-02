clean_ps_from_lulu <- function(ps, res_lulu) {
  otu <- otu_table(res_lulu$curated_table, taxa_are_rows = TRUE)
  tax <- prune_taxa(ps, taxa = res_lulu$curated_otus) %>% tax_table()
  
  phyloseq(otu, tax, ps@sam_data)
}

