library('tidyverse')
library('qvalue')
library('magrittr')
library('gridExtra')
library('annotables')

#Generate extended figure 10B-D

########### LOAD DATA
#### Select samples to process
samps <- read.table("data/metadata.txt", header=TRUE, fill = TRUE, sep = "\t")
samps <- samps[!grepl("K562",samps$tissue),]
samps <- samps[!grepl("direct",samps$protocol),]
samps <- samps[!grepl("experimental",samps$name),]
samps <- samps[!grepl("CVD",samps$name),]
#no vcf
samps <- samps[!samps$path %in% c("FAK44759","FAK41473","FAK49131","FAK46829","FAK49024"),]
#too few features
samps <- samps[!samps$path %in% c("FAK91589"),]
#remove replicates 
samps <- samps[!grepl("rep",samps$name),]
samps <- droplevels(samps)

## load names dict
gff_class <- read_tsv("analysis/flair/new_noSJ_promoters.isoforms_annotation.txt")

########### ANALYSIS
### Analyse ASE events
file_names <- unique(paste0("analysis/ase/flair/",samps$path,"_asts_quant.tsv"))
ase_quant_ont_flair <- lapply(file_names, process_ase_quant_flair, 20, gff_class)
names(ase_quant_ont_flair) <- samps$name

ase_quant_ont_flair[["GTEX_WY7C_heartatrialappendage"]] <- process_ase_quant_flair("analysis/ase/flair/GTEX_WY7C_heartatrialappendage_asts_quant.tsv", 20, gff_class)
ase_quant_ont_flair[["GTEX_13QJ3_muscle"]] <- process_ase_quant_flair("analysis/ase/flair/GTEX_13QJ3_muscle_asts_quant.tsv", 20, gff_class)
ase_quant_ont_flair[["GTEX_14BMU_lung"]] <- process_ase_quant_flair("analysis/ase/flair/GTEX_14BMU_lung_asts_quant.tsv", 20, gff_class)
ase_quant_ont_flair[["GTEX_Q2AG_braincerebellarhemisphere"]] <- process_ase_quant_flair("analysis/ase/flair/GTEX_Q2AG_braincerebellarhemisphere_asts_quant.tsv", 20, gff_class)
ase_quant_ont_flair[["GTEX_Y5LM_liver"]] <- process_ase_quant_flair("analysis/ase/flair/GTEX_Y5LM_liver_asts_quant.tsv", 20, gff_class)
ase_quant_ont_flair[["GTEX_1GN1W_heartatrialappendage"]] <- process_ase_quant_flair("analysis/ase/flair/GTEX_1GN1W_heartatrialappendage_asts_quant.tsv", 20, gff_class)
ase_quant_ont_flair[["GTEX_SM3QNG2_fibs"]] <- process_ase_quant_flair("analysis/ase/flair/GTEX_RWS6_fibs_asts_quant.tsv", 20, gff_class)
ase_quant_ont_flair[["QV44_control"]] <- process_ase_quant_flair("analysis/ase/flair/QV44_control_asts_quant.tsv", 20, gff_class)
ase_quant_ont_flair[["GTEX_ZT9X_muscle"]] <- process_ase_quant_flair("analysis/ase/flair/GTEX_ZT9X_muscle_asts_quant.tsv", 20, gff_class)

all <- list()
sig <- list()
for (i in 1:length(ase_quant_ont_flair)){
  all[[i]] <- length(unique(ase_quant_ont_flair[[i]]$Gene))
}

ase_quant_ont_flair_sig <- lapply(ase_quant_ont_flair, function(x) filter(x, qvalue <= 0.05))
for (i in 1:length(ase_quant_ont_flair_sig)){
  sig[[i]] <- length(unique(ase_quant_ont_flair_sig[[i]]$Gene))
}

all <- unlist(all)
sig <- unlist(sig)

df <- data.frame(all = all,
                 sig = sig,
                 samps)
df$ratio <- df$sig/df$all

ggplot(df, aes(x=all, y=sig, color=tissue)) +
  geom_point(size=4) +
  theme_classic(base_size=14) +
  geom_abline(slope = median(df$ratio), colour = "grey30") +
  ylab("Number of significant genes") + xlab("Number of genes tested") +
  scale_color_manual(values = c("#FFC0CB","#8D5B96","#7776B1","#9773BA","#B873BA","#C893C9",
                                "#FF69B4","#D4A910","#C4625D","#BC3C28","#815375","#0072B5",
                                "#1F854E","#E18726")) +
  theme(panel.grid.major = element_line(colour = "grey90")) +
  ggtitle("ASE")
ggsave("number_sig_ase_by_tissue_lorals_ase_quant_flair_filter.pdf", plot = last_plot(), useDingbats=FALSE,
       scale = 1, width = 9, height = 5, units = "in", dpi = 300)

ase_quant_ont_flair_df <- bind_rows(ase_quant_ont_flair, .id = "names")
ase_quant_ont_flair_df$symbol <- grch38$symbol[match(gsub("\\..*","",ase_quant_ont_flair_df$Gene), grch38$ensgene)]
ase_quant_ont_flair_df$tissue <- samps$tissue[match(ase_quant_ont_flair_df$names,samps$name)]
ase_quant_ont_flair_df$Donor <- samps$sample_name[match(ase_quant_ont_flair_df$names,samps$name)]
write.table(ase_quant_ont_flair_df, "analysis/ase/flair/all_ase_quant_ont_flair_df.txt", col.names = TRUE, 
            row.names = FALSE, quote = FALSE, sep="\t")

##PER TISSUE ANALYSIS
# number of donors in which tissue was tested
ase_quant_ont_flair_df_tiss <- ase_quant_ont_flair_df %>%
  mutate(sig = ifelse(qvalue <= 0.05, "FDR <= 0.05", "FDR > 0.05")) %>%
  select(-names, -variantID, -contig, -position, -refAllele, -altAllele, -refCount, -altCount, -Total_counts,
         -refRatio, -afc, -pvalue, -qvalue) %>%
  unique() %>%
  filter(!tissue %in% c("Pancreas","Breast - Mammary Tissue","Adipose - Subcutaneous",
                        "Brain - Anterior cingulate cortex (BA24)","Brain - Caudate (basal ganglia)")) %>%
  group_by(tissue, Gene, symbol) %>%
  summarise(counts = n()) %>%
  ungroup() %>%
  group_by(counts) %>%
  summarise(counts_sum = n()) %>%
  ungroup() %>%
  mutate(method = "ASE")

ase_quant_ont_flair_df_TS <- ase_quant_ont_flair_df %>%
  filter(!tissue %in% c("Pancreas","Breast - Mammary Tissue","Adipose - Subcutaneous",
                        "Brain - Anterior cingulate cortex (BA24)","Brain - Caudate (basal ganglia)")) %>%
  mutate(sig = ifelse(qvalue <= 0.05, "FDR <= 0.05", "FDR > 0.05")) %>%
  select(-names, -variantID, -contig, -position, -refAllele, -altAllele, -refCount, -altCount, -Total_counts,
         -refRatio, -afc, -pvalue, -qvalue) %>%
  unique() %>%
  group_by(tissue, Gene, symbol) %>%
  filter(n() > 1) %>%
  mutate(observer = toString(unique(sig))) %>%
  filter(observer != "FDR > 0.05") %>%
  group_by(tissue, Gene, symbol) %>%
  mutate(counts_all = n()) %>%
  group_by(tissue, Gene, sig, symbol) %>%
  mutate(counts_cat = n()) %>%
  filter(sig == "FDR <= 0.05") %>%
  filter(counts_cat > 1) %>%
  select(tissue, Gene, symbol, counts_cat, counts_all) %>%
  ungroup()

write.table(ase_quant_ont_flair_df_TS, "analysis/ase/flair/constistent_ase_genes.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")

### Analyse ASTS events quant
###
asts_quant_ont_flair <- lapply(file_names, process_asts_quant_flair, 10, 36, gff_class)
names(asts_quant_ont_flair) <- samps$name

asts_quant_ont_flair[["GTEX_WY7C_heartatrialappendage"]] <- process_asts_quant_flair("ont/flair_filter/GTEX_WY7C_heartatrialappendage_asts_quant.tsv", 10, 36, gff_class)
asts_quant_ont_flair[["GTEX_13QJ3_muscle"]] <- process_asts_quant_flair("ont/flair_filter/GTEX_13QJ3_muscle_asts_quant.tsv", 10, 36, gff_class)
asts_quant_ont_flair[["GTEX_14BMU_lung"]] <- process_asts_quant_flair("ont/flair_filter/GTEX_14BMU_lung_asts_quant.tsv", 10, 36, gff_class)
asts_quant_ont_flair[["GTEX_Q2AG_braincerebellarhemisphere"]] <- process_asts_quant_flair("ont/flair_filter/GTEX_Q2AG_braincerebellarhemisphere_asts_quant.tsv", 10, 36, gff_class)
asts_quant_ont_flair[["GTEX_Y5LM_liver"]] <- process_asts_quant_flair("ont/flair_filter/GTEX_Y5LM_liver_asts_quant.tsv", 10, 36, gff_class)
asts_quant_ont_flair[["GTEX_1GN1W_heartatrialappendage"]] <- process_asts_quant_flair("ont/flair_filter/GTEX_1GN1W_heartatrialappendage_asts_quant.tsv", 10, 36, gff_class)
asts_quant_ont_flair[["GTEX_SM3QNG2_fibs"]] <- process_asts_quant_flair("ont/flair_filter/GTEX_RWS6_fibs_asts_quant.tsv", 10, 36, gff_class)
asts_quant_ont_flair[["QV44_control"]] <- process_asts_quant_flair("ont/flair_filter/QV44_control_asts_quant.tsv", 10, 36, gff_class)
asts_quant_ont_flair[["GTEX_ZT9X_muscle"]] <- process_asts_quant_flair("ont/flair_filter/GTEX_ZT9X_muscle_asts_quant.tsv", 10, 36, gff_class)

all <- list()
sig <- list()
for (i in 1:length(asts_quant_ont_flair)){
  all[[i]] <- length(unique(asts_quant_ont_flair[[i]]$Gene))
}

asts_quant_ont_flair_sig <- lapply(asts_quant_ont_flair, function(x) filter(x, FDR <= 0.05))
for (i in 1:length(asts_quant_ont_flair_sig)){
  sig[[i]] <- length(unique(asts_quant_ont_flair_sig[[i]]$Gene))
}

all <- unlist(all)
sig <- unlist(sig)

df <- data.frame(all = all,
                 sig = sig,
                 samps)
df$ratio <- df$sig/df$all

ggplot(df, aes(x=all, y=sig, color=tissue)) +
  geom_point(size=4) +
  theme_classic(base_size=14) +
  geom_abline(slope = median(df$ratio), colour = "grey30") +
  ylab("Number of significant genes") + xlab("Number of genes tested") +
  scale_color_manual(values = c("#FFC0CB","#8D5B96","#7776B1","#9773BA","#B873BA","#C893C9",
                                "#FF69B4","#D4A910","#C4625D","#BC3C28","#815375","#0072B5",
                                "#1F854E","#E18726")) +
  theme(panel.grid.major = element_line(colour = "grey90")) +
  ggtitle("ASTS")
ggsave("number_sig_asts_by_tissue_lorals_ase_quant_flair_filter.pdf", plot = last_plot(), useDingbats=FALSE,
       scale = 1, width = 9, height = 5, units = "in", dpi = 300)

asts_quant_ont_flair_df <- bind_rows(asts_quant_ont_flair, .id = "names")
asts_quant_ont_flair_df$symbol <- grch38$symbol[match(gsub("\\..*","",asts_quant_ont_flair_df$Gene), grch38$ensgene)]
asts_quant_ont_flair_df$tissue <- samps$tissue[match(asts_quant_ont_flair_df$names,samps$name)]
asts_quant_ont_flair_df$Donor <- samps$sample_name[match(asts_quant_ont_flair_df$names,samps$name)]
write.table(asts_quant_ont_flair_df, "analysis/ase/flair/all_asts_quant_ont_flair_df.txt", col.names = TRUE, 
            row.names = FALSE, quote = FALSE, sep="\t")

##PER TISSUE ANALYSIS
asts_quant_ont_flair_df_tiss <- asts_quant_ont_flair_df %>%
  mutate(sig = ifelse(FDR <= 0.05, "FDR <= 0.05", "FDR > 0.05")) %>%
  select(-names, -variantID, -contig, -position, -refAllele, -altAllele, -refCount, -altCount, -Transcript,
         -sum, -cohen, -pvalue, -FDR) %>%
  unique() %>%
  filter(!tissue %in% c("Pancreas","Breast - Mammary Tissue","Adipose - Subcutaneous",
                        "Brain - Anterior cingulate cortex (BA24)","Brain - Caudate (basal ganglia)")) %>%
  group_by(tissue, Gene, symbol) %>%
  summarise(counts = n()) %>%
  ungroup() %>%
  group_by(counts) %>%
  summarise(counts_sum = n()) %>%
  ungroup() %>%
  mutate(method = "ASTS")

by_tissue <- rbind(ase_quant_ont_flair_df_tiss, asts_quant_ont_flair_df_tiss)
ggplot(by_tissue, aes(x=counts, y=counts_sum, color=method)) +
  geom_point(size=3, alpha=0.8) +
  theme_classic(base_size=14) +
  xlab("# donors") +
  ylab("# genes") +
  scale_color_manual(values=c("#a8479a","#99a745")) 
ggsave("number_consistent_allelic_flair_filter.pdf", plot = last_plot(), useDingbats=FALSE,
       scale = 1, width = 6, height = 6, units = "in", dpi = 300)

asts_quant_ont_flair_df_TS <- asts_quant_ont_flair_df %>%
  mutate(sig = ifelse(FDR <= 0.05, "FDR <= 0.05", "FDR > 0.05")) %>%
  group_by(tissue, Gene, symbol) %>%
  mutate(no_transcripts = length(unique(Transcript))) %>%
  select(-names, -variantID, -contig, -position, -refAllele, -altAllele, -refCount, -altCount, -Transcript,
         -sum, -cohen, -pvalue, -FDR) %>%
  unique() %>%
  filter(!tissue %in% c("Pancreas","Breast - Mammary Tissue","Adipose - Subcutaneous",
                        "Brain - Anterior cingulate cortex (BA24)","Brain - Caudate (basal ganglia)")) %>%
  filter(n() > 1) %>%
  mutate(observer = toString(unique(sig))) %>%
  filter(observer != "FDR > 0.05") %>%
  group_by(tissue, Gene, symbol) %>%
  mutate(counts_all = n()) %>%
  group_by(tissue, Gene, sig, symbol) %>%
  mutate(counts_cat = n()) %>%
  filter(sig == "FDR <= 0.05") %>%
  filter(counts_cat > 1) %>%
  group_by(tissue, Gene, symbol, counts_cat, counts_all) %>%
  summarise(no_transcripts = max(no_transcripts)) %>%
  ungroup()
write.table(asts_quant_ont_flair_df_TS, "analysis/ase/flair/constistent_asts_genes.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")

asts_quant_ont_flair_df_gene <- asts_quant_ont_flair_df %>%
  group_by(names, Gene) %>%
  mutate(Ref_count = sum(refCount), Alt_count = sum(altCount), number_of_isoforms = length(Transcript)) %>%
  ungroup() %>%
  dplyr::select(names, Gene, number_of_isoforms, variantID, Ref_count, Alt_count, sum, cohen,  pvalue, FDR) %>%
  unique() %>%
  mutate(sig = ifelse(FDR <= 0.05, "SIG","NON-SIG"))

toplot_trans <- asts_quant_ont_flair_df_gene %>%
  mutate(number_of_isoforms = ifelse(number_of_isoforms >= 7, 7, number_of_isoforms)) %>%
  group_by(number_of_isoforms) %>%
  summarise(allelic_n = n()) %>%
  ungroup() 

ggplot(toplot_trans, aes(x=number_of_isoforms, y=allelic_n)) +
  geom_bar(stat="identity", position = "dodge") +
  theme_classic(base_size = 14) +
  ylab("# genes") + xlab("# isoforms per gene") +
  scale_fill_manual(values=c("#a3a1a2","#d66227")) 
ggsave("number_isoforms_asts.pdf", plot = last_plot(), useDingbats=FALSE,
       scale = 1, width = 4, height = 6, units = "in", dpi = 300)
