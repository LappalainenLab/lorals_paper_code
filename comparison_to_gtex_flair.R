#### FLAIR ANALYSIS IN GTEX
library(tidyverse)
library(qvalue)
library(seqminer)

### Analyse ASE QUANT FLAIR events
my_ase_quant_ont_flair_df <- read_tsv("analysis/ase/flair/all_ase_quant_ont_flair_df.txt")
ase_quant_ont_flair <- split( my_ase_quant_ont_flair_df , f = my_ase_quant_ont_flair_df$names )

### ASE QUANT in EQTLs
all_geno <- list()
for (i in 1:length(ase_quant_ont_flair)){
  print(paste0("Processing sample ",ase_quant_ont_flair[[i]]$names[1]))
  name <- gsub("-", "", ase_quant_ont_flair[[i]]$tissue[1], fixed = TRUE)
  name <- gsub(")", "", name, fixed = TRUE)
  name <- gsub("(", "", name, fixed = TRUE)
  name <- str_replace(gsub("\\s+", "_", str_trim(name)), "sk", "Sk")
  eqtls_gtex <- read.table(paste0("data",
                                  name,".v8.egenes.txt.gz"), header=TRUE)
  eqtls_gtex <- eqtls_gtex[gsub("\\..*","",eqtls_gtex$gene_id) %in% gsub("\\..*","",ase_quant_ont_flair[[i]]$Gene),]
  eqtls_gtex <- eqtls_gtex[eqtls_gtex$qval <= 0.05,]
  afc_gtex <- read.table(paste0("data",name,".aFC.txt.gz"), header=TRUE)
  afc_gtex <- afc_gtex[abs(afc_gtex$log2_aFC_lower) >= 0.2,]
  eqtls_gtex <- eqtls_gtex[eqtls_gtex$gene_id %in% afc_gtex$pid,]
  print(paste0("Found ", dim(eqtls_gtex)[1]," genes"))
  common_df <- ase_quant_ont_flair[[i]][gsub("\\..*","",ase_quant_ont_flair[[i]]$Gene) %in% gsub("\\..*","",eqtls_gtex$gene_id),]
  eqtls_gtex$ranges <- paste0(eqtls_gtex$chr,":",eqtls_gtex$variant_pos - 1,"-",eqtls_gtex$variant_pos)
  donor <- gsub("_", "-", common_df$Donor[1], fixed = TRUE)
  fileName = paste0("analysis/vcf/",donor,"_all_variants.vcf.gz")
  my_ranges <- paste0( eqtls_gtex$ranges, collapse=",")
  snp <- tabix.read.table(fileName, my_ranges)
  snp$names <- ase_quant_ont_flair[[i]]$names[1]
  snp$GENE_ID <- eqtls_gtex$gene_id[match(gsub("_b38","",snp$ID), paste0(eqtls_gtex$chr,"_",eqtls_gtex$variant_pos,"_",eqtls_gtex$ref,"_",eqtls_gtex$alt))]
  all_geno[[i]] <- snp
}

colnames_geno <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","GENOTYPE","names","GENE_ID") 
all_geno <- lapply(all_geno, setNames, colnames_geno)
all_geno_df <- bind_rows(all_geno, .id = NULL)
all_geno_df %<>%
  dplyr::select(-CHROM, -POS, -REF, -ALT, -QUAL, -FILTER, -INFO, -FORMAT) %>%
  mutate(GT = ifelse(GENOTYPE %in% c("0|0","1|1"), "HOMOZYGOUS","HETEROZYGOUS"))

all_geno_df_ase_quant_flair <- merge(all_geno_df, my_ase_quant_ont_flair_df, by.x=c("names","GENE_ID"), by.y=c("names","Gene"), keep="all")
all_geno_df_ase_quant_flair$is_sig <- ifelse(all_geno_df_ase_quant_flair$pvalue <= 0.05, "Yes", "No")
all_geno_df_ase_quant_flair %<>%
  group_by(GT, is_sig) %>%
  mutate(proportion1 = n()) %>%
  ungroup() %>%
  group_by(GT) %>%
  mutate(proportion2 = n()) %>%
  ungroup() %>%
  mutate(proportion = proportion1/proportion2) %>%
  dplyr::select(GT, is_sig, proportion, proportion1) %>%
  unique() %>%
  arrange(GT, is_sig)
  
eqtl_in_ase <-
  matrix(all_geno_df_ase_quant_flair$proportion1,
         nrow = 2,
         dimnames = list(Sig = c("No", "Yes"),
                         Geno = c("HET", "HOM")))
eqtl_in_ase <- eqtl_in_ase[ nrow(eqtl_in_ase):1, ]
f1 <- fisher.test(eqtl_in_ase, alternative = "two.sided")
all_geno_df_ase_quant_flair <- all_geno_df_ase_quant_flair[all_geno_df_ase_quant_flair$is_sig=="Yes",]
p1 <- ggplot(all_geno_df_ase_quant_flair, aes(x=is_sig, fill=GT, y=proportion*100)) + 
  geom_bar(position='dodge', stat="identity") +
  theme_classic(base_size=14) +
  xlab("") + ylab("% of sig. ASE quant events") +
  scale_fill_manual(values=c("#e69f00","#56b4e9")) +
  ggtitle("eQTLs") +
  ylim(0,50) +
  annotate(geom="text", x=1, y=50, label=paste0("p = ", f1$p.value)) +
  annotate(geom="text", x=1, y=40, label=paste0("OR = ", round(f1$estimate,3))) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "bottom")
ggsave("all_eqtls_in_ase_quant_flair.pdf", plot = p1,
       scale = 1, width = 2.5, height = 3.75, units = "in", dpi = 300)

### ANALYSE ASE QUANT IN SQTLs
all_geno <- list()
for (i in 1:length(ase_quant_ont_flair)){
  print(paste0("Processing sample ",ase_quant_ont_flair[[i]]$names[1]))
  name <- gsub("-", "", ase_quant_ont_flair[[i]]$tissue[1], fixed = TRUE)
  name <- gsub(")", "", name, fixed = TRUE)
  name <- gsub("(", "", name, fixed = TRUE)
  name <- str_replace(gsub("\\s+", "_", str_trim(name)), "sk", "Sk")
  sqtls_gtex <- read.table(paste0("data/",name,".v8.sgenes.txt.gz"), header=TRUE)
  sqtls_gtex <- sqtls_gtex[sqtls_gtex$group_id %in% ase_quant_ont_flair[[i]]$Gene,]
  sqtls_gtex <- sqtls_gtex[sqtls_gtex$qval <= 0.05,]
  print(paste0("Found ", dim(sqtls_gtex)[1]," genes"))
  common_df <- ase_quant_ont_flair[[i]][ase_quant_ont_flair[[i]]$Gene %in% sqtls_gtex$group_id,]
  sqtls_gtex$ranges <- paste0(sqtls_gtex$chr,":",sqtls_gtex$variant_pos - 1,"-",sqtls_gtex$variant_pos)
  donor <- gsub("_", "-", common_df$Donor[1], fixed = TRUE)
  fileName = paste0("analysis/vcf/",donor,"_all_variants.vcf.gz")
  my_ranges <- paste0( sqtls_gtex$ranges, collapse=",")
  snp <- tabix.read.table(fileName, my_ranges)
  snp$names <- ase_quant_ont_flair[[i]]$names[1]
  snp$GENE_ID <- sqtls_gtex$group_id[match(snp$ID, sqtls_gtex$variant_id)]
  all_geno[[i]] <- snp
}

colnames_geno <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","GENOTYPE","names","GENE_ID") 
all_geno <- lapply(all_geno, setNames, colnames_geno)
all_geno_df <- bind_rows(all_geno, .id = NULL)
all_geno_df %<>%
  dplyr::select(-CHROM, -POS, -REF, -ALT, -QUAL, -FILTER, -INFO, -FORMAT) %>%
  mutate(GT = ifelse(GENOTYPE %in% c("0|0","1|1"), "HOMOZYGOUS","HETEROZYGOUS"))

all_geno_df_ase_quant_s_flair <- merge(all_geno_df, my_ase_quant_ont_flair_df, by.x=c("names","GENE_ID"), by.y=c("names","Gene"), keep="all")

all_geno_df_ase_quant_s_flair$is_sig <- ifelse(all_geno_df_ase_quant_s_flair$pvalue <= 0.05, "Yes", "No")
all_geno_df_ase_quant_s_flair %<>%
  group_by(GT, is_sig) %>%
  mutate(proportion1 = n()) %>%
  ungroup() %>%
  group_by(GT) %>%
  mutate(proportion2 = n()) %>%
  ungroup() %>%
  mutate(proportion = proportion1/proportion2) %>%
  dplyr::select(GT, is_sig, proportion, proportion1) %>%
  unique() %>%
  arrange(GT, is_sig)

sqtl_in_ase <-
  matrix(all_geno_df_ase_quant_s_flair$proportion1,
         nrow = 2,
         dimnames = list(Sig = c("No", "Yes"),
                         Geno = c("HET", "HOM")))
sqtl_in_ase <- sqtl_in_ase[ nrow(sqtl_in_ase):1, ]
f2 <- fisher.test(sqtl_in_ase, alternative = "two.sided")
all_geno_df_ase_quant_s_flair <- all_geno_df_ase_quant_s_flair[all_geno_df_ase_quant_s_flair$is_sig=="Yes",]
p2 <- ggplot(all_geno_df_ase_quant_s_flair, aes(x=is_sig, fill=GT, y=proportion*100)) + 
  geom_bar(position='dodge', stat="identity") +
  theme_classic(base_size=14) +
  xlab("") + ylab("% of sig. ASE quant events") +
  scale_fill_manual(values=c("#e69f00","#56b4e9")) +
  ggtitle("sQTLs") +
  ylim(0,50) +
  annotate(geom="text", x=1, y=50, label=paste0("p = ", f2$p.value)) +
  annotate(geom="text", x=1, y=40, label=paste0("OR = ", round(f2$estimate,3))) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "bottom")
ggsave("all_sqtls_in_ase_quant_flair.pdf", plot = p2,
       scale = 1, width = 2.5, height = 3.75, units = "in", dpi = 300)

### Analyse ASTS quant flair events in sQTLs
my_asts_quant_ont_flair_df <- read_tsv("analysis/ase/flair/all_asts_quant_ont_flair_df.txt")
asts_quant_ont_flair <- split( my_asts_quant_ont_flair_df , f = my_asts_quant_ont_flair_df$names )

all_geno <- list()
for (i in 1:length(asts_quant_ont_flair)){
  print(paste0("Processing sample ",asts_quant_ont_flair[[i]]$names[1]))
  name <- gsub("-", "", asts_quant_ont_flair[[i]]$tissue[1], fixed = TRUE)
  name <- gsub(")", "", name, fixed = TRUE)
  name <- gsub("(", "", name, fixed = TRUE)
  name <- str_replace(gsub("\\s+", "_", str_trim(name)), "sk", "Sk")
  sqtls_gtex <- read.table(paste0("data/",name,".v8.sgenes.txt.gz"), header=TRUE)
  sqtls_gtex <- sqtls_gtex[sqtls_gtex$group_id %in% asts_quant_ont_flair[[i]]$Gene,]
  sqtls_gtex <- sqtls_gtex[sqtls_gtex$qval <= 0.05,]
  print(paste0("Found ", dim(sqtls_gtex)[1]," genes"))
  common_df <- asts_quant_ont_flair[[i]][asts_quant_ont_flair[[i]]$Gene %in% sqtls_gtex$group_id,]
  sqtls_gtex$ranges <- paste0(sqtls_gtex$chr,":",sqtls_gtex$variant_pos - 1,"-",sqtls_gtex$variant_pos)
  donor <- gsub("_", "-", common_df$Donor[1], fixed = TRUE)
  fileName = paste0("analysis/vcf/",donor,"_all_variants.vcf.gz")
  my_ranges <- paste0( sqtls_gtex$ranges, collapse=",")
  snp <- tabix.read.table(fileName, my_ranges)
  snp$names <- asts_quant_ont_flair[[i]]$names[1]
  snp$Gene <- sqtls_gtex$group_id[match(snp$ID, sqtls_gtex$variant_id)]
  all_geno[[i]] <- snp
}

colnames_geno <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","GENOTYPE","names","Gene") 
all_geno <- lapply(all_geno, setNames, colnames_geno)
for (i in 1:length(all_geno)){
  all_geno[[i]]$REF <- ifelse(all_geno[[i]]$REF==TRUE, "T", all_geno[[i]]$REF)
  all_geno[[i]]$ALT <- ifelse(all_geno[[i]]$ALT==TRUE, "T", all_geno[[i]]$ALT)
}
all_geno_df <- bind_rows(all_geno, .id = NULL)
all_geno_df %<>%
  dplyr::select(-CHROM, -POS, -REF, -ALT, -QUAL, -FILTER, -INFO, -FORMAT) %>%
  mutate(GT = ifelse(GENOTYPE %in% c("0|0","1|1"), "HOMOZYGOUS","HETEROZYGOUS"))

all_geno_df_asts_flair <- merge(all_geno_df, my_asts_quant_ont_flair_df, by=c("names","Gene"), keep="all")
all_geno_df_asts_flair$is_sig <- ifelse(all_geno_df_asts_flair$pvalue <= 0.05, "Yes", "No")
all_geno_df_asts_flair <- na.omit(all_geno_df_asts_flair)
all_geno_df_asts_flair %<>%
  group_by(GT, is_sig) %>%
  mutate(proportion1 = n()) %>%
  ungroup() %>%
  group_by(GT) %>%
  mutate(proportion2 = n()) %>%
  ungroup() %>%
  mutate(proportion = proportion1/proportion2) %>%
  dplyr::select(GT, is_sig, proportion, proportion1) %>%
  unique() %>%
  arrange(GT, is_sig)

sqtl_in_asts <-
  matrix(all_geno_df_asts_flair$proportion1,
         nrow = 2,
         dimnames = list(Sig = c("No", "Yes"),
                         Geno = c("HET", "HOM")))
sqtl_in_asts <- sqtl_in_asts[ nrow(sqtl_in_asts):1, ]
f3 <- fisher.test(sqtl_in_asts, alternative = "two.sided")
all_geno_df_asts_flair <- all_geno_df_asts_flair[all_geno_df_asts_flair$is_sig=="Yes",]
p3 <- ggplot(all_geno_df_asts_flair, aes(x=is_sig, fill=GT, y=proportion*100)) + 
  geom_bar(position='dodge', stat="identity") +
  theme_classic(base_size=14) +
  xlab("") + ylab("% of sig. ASTS quant events") +
  scale_fill_manual(values=c("#e69f00","#56b4e9")) +
  ggtitle("sQTLs") +
  ylim(0,50) +
  annotate(geom="text", x=1, y=50, label=paste0("p = ", f3$p.value)) +
  annotate(geom="text", x=1, y=40, label=paste0("OR = ", round(f3$estimate,3))) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "bottom")
ggsave("all_sqtls_in_asts_quant_flair.pdf", plot = p3,
       scale = 1, width = 2.5, height = 3.75, units = "in", dpi = 300)

### Analyse ASTS events in eQTLs
all_geno <- list()
for (i in 1:length(asts_quant_ont_flair)){
  print(paste0("Processing sample ",asts_quant_ont_flair[[i]]$names[1]))
  name <- gsub("-", "", asts_quant_ont_flair[[i]]$tissue[1], fixed = TRUE)
  name <- gsub(")", "", name, fixed = TRUE)
  name <- gsub("(", "", name, fixed = TRUE)
  name <- str_replace(gsub("\\s+", "_", str_trim(name)), "sk", "Sk")
  eqtls_gtex <- read.table(paste0("data/", name,".v8.egenes.txt.gz"), header=TRUE)
  eqtls_gtex <- eqtls_gtex[eqtls_gtex$gene_id %in% asts_quant_ont_flair[[i]]$Gene,]
  eqtls_gtex <- eqtls_gtex[eqtls_gtex$qval <= 0.05,]
  afc_gtex <- read.table(paste0("data/", name,".aFC.txt.gz"), header=TRUE)
  afc_gtex <- afc_gtex[abs(afc_gtex$log2_aFC_lower) >= 0.2,]
  eqtls_gtex <- eqtls_gtex[eqtls_gtex$gene_id %in% afc_gtex$pid,]
  print(paste0("Found ", dim(eqtls_gtex)[1]," genes"))
  if (dim(eqtls_gtex)[1] == 0){
    next
  }
  common_df <- asts_quant_ont_flair[[i]][asts_quant_ont_flair[[i]]$Gene %in% eqtls_gtex$gene_id,]
  eqtls_gtex$ranges <- paste0(eqtls_gtex$chr,":",eqtls_gtex$variant_pos - 1,"-",eqtls_gtex$variant_pos)
  donor <- gsub("_", "-", common_df$Donor[1], fixed = TRUE)
  fileName = paste0("analysis/vcf/",donor,"_all_variants.vcf.gz")
  my_ranges <- paste0( eqtls_gtex$ranges, collapse=",")
  snp <- tabix.read.table(fileName, my_ranges)
  snp$names <- asts_quant_ont_flair[[i]]$names[1]
  snp$Gene <- eqtls_gtex$gene_id[match(gsub("_b38","",snp$ID), paste0(eqtls_gtex$chr,"_",eqtls_gtex$variant_pos,"_",eqtls_gtex$ref,"_",eqtls_gtex$alt))]
  all_geno[[i]] <- snp
}

all_geno <- all_geno[lapply(all_geno,length)>0]
colnames_geno <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","GENOTYPE","names","Gene") 
all_geno <- lapply(all_geno, setNames, colnames_geno)
for (i in 1:length(all_geno)){
  all_geno[[i]]$REF <- ifelse(all_geno[[i]]$REF==TRUE, "T", all_geno[[i]]$REF)
  all_geno[[i]]$ALT <- ifelse(all_geno[[i]]$ALT==TRUE, "T", all_geno[[i]]$ALT)
}
all_geno_df <- bind_rows(all_geno, .id = NULL)
all_geno_df %<>%
  dplyr::select(-CHROM, -POS, -REF, -ALT, -QUAL, -FILTER, -INFO, -FORMAT) %>%
  mutate(GT = ifelse(GENOTYPE %in% c("0|0","1|1"), "HOMOZYGOUS","HETEROZYGOUS"))

all_geno_df_asts_e_flair <- merge(all_geno_df, my_asts_quant_ont_flair_df, by=c("names","Gene"), keep="all")
all_geno_df_asts_e_flair$is_sig <- ifelse(all_geno_df_asts_e_flair$pvalue <= 0.05, "Yes", "No")
all_geno_df_asts_e_flair <- na.omit(all_geno_df_asts_e_flair)
all_geno_df_asts_e_flair %<>%
  group_by(GT, is_sig) %>%
  mutate(proportion1 = n()) %>%
  ungroup() %>%
  group_by(GT) %>%
  mutate(proportion2 = n()) %>%
  ungroup() %>%
  mutate(proportion = proportion1/proportion2) %>%
  dplyr::select(GT, is_sig, proportion, proportion1) %>%
  unique() %>%
  arrange(GT, is_sig)

eqtl_in_asts <-
  matrix(all_geno_df_asts_e_flair$proportion1,
         nrow = 2,
         dimnames = list(Sig = c("No", "Yes"),
                         Geno = c("HET", "HOM")))
eqtl_in_asts <- eqtl_in_asts[ nrow(eqtl_in_asts):1, ]
f4 <- fisher.test(eqtl_in_asts, alternative = "two.sided")
all_geno_df_asts_e_flair <- all_geno_df_asts_e_flair[all_geno_df_asts_e_flair$is_sig=="Yes",]
p4 <- ggplot(all_geno_df_asts_e_flair, aes(x=is_sig, fill=GT, y=proportion*100)) + 
  geom_bar(position='dodge', stat="identity") +
  theme_classic(base_size=14) +
  xlab("") + ylab("% of sig. ASTS quant events") +
  scale_fill_manual(values=c("#e69f00","#56b4e9")) +
  ggtitle("eQTLs") +
  ylim(0,50) +
  annotate(geom="text", x=1, y=50, label=paste0("p = ", f4$p.value)) +
  annotate(geom="text", x=1, y=40, label=paste0("OR = ", round(f4$estimate,3))) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "bottom")
ggsave("all_eqtls_in_asts_quant_flair.pdf", plot = p4,
       scale = 1, width = 2.5, height = 3.75, units = "in", dpi = 300)

all_geno_df_asts_e_flair$ONT <- "ASTS"
all_geno_df_asts_flair$ONT <- "ASTS"
all_geno_df_ase_quant_flair$ONT <- "ASE"
all_geno_df_ase_quant_s_flair$ONT <- "ASE"
all_geno_df_asts_e_flair$GTEx <- "eQTLs"
all_geno_df_asts_flair$GTEx <- "sQTLs"
all_geno_df_ase_quant_flair$GTEx <- "eQTLs"
all_geno_df_ase_quant_s_flair$GTEx <- "sQTLs"
all_together <- rbind(all_geno_df_asts_e_flair, all_geno_df_asts_flair, all_geno_df_ase_quant_flair, all_geno_df_ase_quant_s_flair)
ggplot(all_together, aes(x=GTEx, fill=GT, y=proportion*100)) + 
  geom_bar(position='dodge', stat="identity") +
  theme_classic(base_size=14) +
  facet_grid(.~ONT) +
  xlab("") + ylab("% of sig. events") +
  scale_fill_manual(values=c("#e69f00","#56b4e9")) +
  ylim(0,40) +
  theme(legend.position = "bottom")
ggsave("gtex_qtls_in_allelic_flair.pdf", plot = last_plot(),
       scale = 1, width = 5, height = 3.75, units = "in", dpi = 300)
