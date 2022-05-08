library(tidyverse)
library(seqminer)

# ENLOC results on EUR individuals
process_enloc_table <- function(tissue) {
  file_names <- list.files(path = "data/", pattern = tissue)
  file_names_path <- paste0("data",file_names)
  enloc_tables <- lapply(file_names_path, read.table, header = T, sep = "\t", stringsAsFactors = F, fill = TRUE, col.names = paste0("V",seq_len(6)))
  file_names <- gsub(paste0("__",tissue,".enloc.rst.gz"),"",file_names)
  names(enloc_tables) <- file_names
  enloc_tables <- enloc_tables[sapply(enloc_tables, function(x) dim(x)[1]) > 0]
  enloc_tables <- lapply(enloc_tables, function(df) mutate_at(df, .vars = c("V5"), as.character))
  enloc_tables_df <- bind_rows(enloc_tables, .id = "GWAS")
  colnames(enloc_tables_df) <- c("GWAS","gwas_locus", "molecular_qtl_trait", "locus_gwas_pip", "locus_rcp", "lead_coloc_SNP", "lead_snp_scp")
  enloc_tables_df <- enloc_tables_df[enloc_tables_df$locus_rcp > 0.5,]
  intron_gene_map <- read.table(paste0("data/intron_gene_map_",tissue,".txt.gz"), header=TRUE)
  enloc_tables_df <- merge(enloc_tables_df, intron_gene_map, by.x = "molecular_qtl_trait", by.y="intron_id")
  return(enloc_tables_df)
}


#Read in ASTS
asts_quant_ont_flair_df <- read_tsv("analysis/ase/flair/all_asts_quant_ont_flair_df.txt")
asts_quant_ont_flair_df$tissue <- gsub("-", "", asts_quant_ont_flair_df$tissue, fixed = TRUE)
asts_quant_ont_flair_df$tissue <- gsub(")", "", asts_quant_ont_flair_df$tissue, fixed = TRUE)
asts_quant_ont_flair_df$tissue <- gsub("(", "", asts_quant_ont_flair_df$tissue, fixed = TRUE)
asts_quant_ont_flair_df$tissue <- str_replace(gsub("\\s+", "_", str_trim(asts_quant_ont_flair_df$tissue)), "sk", "Sk")

for (i in unique(asts_quant_ont_flair_df$tissue)){
  print(i)
  enloc_df <- process_enloc_table(i)
  print(length(unique(enloc_df$gene_id)))
  asts_df <- dplyr::filter(asts_quant_ont_flair_df, tissue == i)
  print(length(unique(asts_df$Gene)))
  asts_df_sig <- dplyr::filter(asts_df, FDR <= 0.05)
  print(length(unique(asts_df_sig$Gene)))
  enloc_asts_overlap <- dplyr::filter(asts_df, Gene %in% enloc_df$gene_id)
  print(length(unique(enloc_asts_overlap$Gene)))
  enloc_asts_overlap_sig <- dplyr::filter(enloc_asts_overlap, FDR <= 0.05)
  print(length(unique(enloc_asts_overlap_sig$Gene)))
}