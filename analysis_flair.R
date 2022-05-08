library(tidyverse)
library(magrittr)
library(pheatmap)
library(RColorBrewer)
library(UpSetR)

#Generate Extended Figure 6A,B,C,

## set parameters for pheatmap
breaksList = seq(0, 1, by = 0.01)

my_colour = list(
  tissue = c(Adipose_Subcutaneous = "#ffc0cb",
             Brain_Anterior_cingulate_cortex_BA24 = "#8d5b96",
             Brain_Caudate_basal_ganglia = "#7776b1",
             Brain_Cerebellar_Hemisphere = "#9773ba",
             Brain_Frontal_Cortex_BA9 = "#b873ba",
             Brain_Putamen_basal_ganglia = "#c893c9",
             Breast_Mammary_Tissue = "#ff69b4",
             Cells_Cultured_fibroblasts = "#d4a910",
             Heart_Atrial_Appendage = "#c4625d",
             Heart_Left_Ventricle = "#bc3c28",
             K562 = "#B09638",
             Liver = "#815375",
             Lung = "#0072b5",
             Muscle_Skeletal = "#1f854e",
             Pancreas = "#e18726"))

## Read in tracking info to determine novel transcripts
gff_class <- read_tsv("analysis/flair/new_noSJ_promoters.isoforms_annotation.txt")
gff_class_novel <- gff_class[gff_class$Class_code != "=",]
gff_class_annot <- gff_class[gff_class$Class_code == "=",]

## Read in sample info
samps <- read.table("data/metadata.txt",
                    header=TRUE, fill = TRUE, sep = "\t")
samps <- samps[samps$name != "GTEX_T5JC_brainfrontalcortex",]
samps <- samps[!grepl("CVD",samps$name),]
samps <- samps[!grepl("rep",samps$name),]
samps <- samps[!grepl("experimental",samps$name),]
samps <- samps[!grepl("direct",samps$name),]
samps <- samps[!grepl("K562",samps$name),]
samps$group <- gsub(" - ","_",samps$tissue)
samps$group <- gsub(" \\(","_",samps$group)
samps$group <- gsub("\\)","",samps$group)
samps$group <- gsub(" ","_",samps$group)

ont_nano_tpm_df <- read.table("analysis/count_tables/trans_flair.tpm.txt", header = TRUE)
rownames(ont_nano_tpm_df) <- ont_nano_tpm_df$transcript_id
ont_nano_tpm_df <- ont_nano_tpm_df[,-1]
ont_nano_tpm_df <- ont_nano_tpm_df[,colnames(ont_nano_tpm_df) %in% samps$name]
ont_nano_tpm_df <- ont_nano_tpm_df[rowSums(ont_nano_tpm_df) > 0, ]
ont_nano_tpm_df <- na.omit(ont_nano_tpm_df)

### REPLACE COLUMNS WITH THE MERGED COUNTS
merged_file_names <- c("GTEX_WY7C_heartatrialappendage","GTEX_13QJ3_muscle","GTEX_14BMU_lung","GTEX_ZT9X_muscle",
                       "GTEX_Q2AG_braincerebellarhemisphere","GTEX_Y5LM_liver","GTEX_1GN1W_heartatrialappendage",
                       "GTEX_SM3QNG2_fibs","QV44_control")
merged_files <- file.path("analysis/count_tables/",paste0(merged_file_names,".tpm.txt"))
names(merged_files) <- merged_file_names

merged_ont_nano_tpm = lapply(merged_files, read_tsv, col_names=TRUE)
merged_ont_nano_tpm_df <- bind_rows(merged_ont_nano_tpm, .id = "names")
merged_ont_nano_tpm_df %<>%
  spread(names, counts)
merged_ont_nano_tpm_df[is.na(merged_ont_nano_tpm_df)] <- 0

matches <- unique (grep(paste(merged_file_names,collapse="|"), 
                        samps$name, value=TRUE))

ont_nano_tpm_uniq = ont_nano_tpm_df
for (i in matches) {
  ont_nano_tpm_uniq %<>%
    dplyr::select(-i)
}

ont_nano_tpm_uniq$transcript_id <- rownames(ont_nano_tpm_uniq)

ont_nano_tpm_uniq <- plyr::join_all(list(ont_nano_tpm_uniq,merged_ont_nano_tpm_df), by = "transcript_id", type = 'full')
rownames(ont_nano_tpm_uniq) <- ont_nano_tpm_uniq$transcript_id 
ont_nano_tpm_uniq$transcript_id <- NULL

ont_nano_tpm_uniq[is.na(ont_nano_tpm_uniq)] <- 0
ont_nano_tpm_uniq <- ont_nano_tpm_uniq[rowSums(ont_nano_tpm_uniq) > 0, ]
ont_nano_tpm_uniq <- na.omit(ont_nano_tpm_uniq)

write.table(ont_nano_tpm_uniq, "analysis/count_tables/trans_flair_reps.tpm.txt", col.names = TRUE,
            row.names = TRUE, quote = FALSE, sep = "\t")

# PCA of annotated transcripts
ont_nano_tpm_uniq_filter <- ont_nano_tpm_uniq[ rowSums(ont_nano_tpm_uniq >= 5) >= 3, ]

ont_nano_tpm_uniq_filter_annot <- ont_nano_tpm_uniq_filter[rownames(ont_nano_tpm_uniq_filter) %in% gff_class_annot$fish,]
highly_variable_lcpm <- log10(ont_nano_tpm_uniq_filter_annot+1)
genes.pca <- prcomp(t(highly_variable_lcpm), center = TRUE, scale. = TRUE)
percentVar <- round(100 * apply(genes.pca$x, 2, var)/sum(apply(genes.pca$x, 2, var)),1)
pcs <- as.data.frame(genes.pca$x[, c(1:6)])
pcs$name <- rownames(pcs)
rownames(pcs) <- NULL

pcs <- merge(pcs, samps, by = "name")
ggplot(pcs, aes(PC1, PC2)) +
  geom_point(aes(fill=tissue),colour="grey20",pch=21, size=5) +
  theme_classic(base_size =14) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  guides(fill = guide_legend(override.aes = list(size=6))) +
  scale_fill_manual(values = c("#ffc0cb","#8d5b96","#7776b1","#9773ba","#b873ba","#c893c9",
                               "#ff69b4","#d4a910","#c4625d","#bc3c28","#815375",
                               "#0072b5", "#1f854e","#e18726"))
ggsave("transcripts_PCA_all_min5in3_annot_flair.pdf", plot = last_plot(),
       width = 7.8, height = 5.7, units = "in", useDingbats=FALSE)

ont_tpm_df_filter_cor1 <- cor(ont_nano_tpm_uniq_filter_annot, method = "spearman")
sampleDists <- dist(ont_tpm_df_filter_cor1, method = "euclidean")
sampleDistMatrix <- as.matrix(ont_tpm_df_filter_cor1)
my_sample_col <- as.data.frame(samps[samps$name %in% colnames(ont_tpm_df_filter_cor1),c("group")])
row.names(my_sample_col) <- samps[samps$name %in% colnames(ont_tpm_df_filter_cor1),]$name
colnames(my_sample_col) <- c("tissue")

pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists, show_rownames = FALSE,
         clustering_distance_cols=sampleDists, show_colnames = FALSE,
         annotation_colors = my_colour,
         col=colorRampPalette(brewer.pal(9, "BrBG"))(length(breaksList)),
         annotation_col = my_sample_col, drop_levels = TRUE,
         breaks=breaksList)

# PCA of novel transcripts 
ont_nano_tpm_uniq_filter_novel <- ont_nano_tpm_uniq_filter[!rownames(ont_nano_tpm_uniq_filter) %in% gff_class_annot$fish,]

highly_variable_lcpm <- log10(ont_nano_tpm_uniq_filter_novel+1)
genes.pca <- prcomp(t(highly_variable_lcpm), center = TRUE, scale. = TRUE)
percentVar <- round(100 * apply(genes.pca$x, 2, var)/sum(apply(genes.pca$x, 2, var)),1)
pcs <- as.data.frame(genes.pca$x[, c(1:6)])
pcs$name <- rownames(pcs)
rownames(pcs) <- NULL

pcs <- merge(pcs, samps, by = "name")
ggplot(pcs, aes(PC1, PC2)) +
  geom_point(aes(fill=tissue),colour="grey20",pch=21, size=5) +
  theme_classic(base_size =14) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  guides(fill = guide_legend(override.aes = list(size=6))) +
  scale_fill_manual(values = c("#ffc0cb","#8d5b96","#7776b1","#9773ba","#b873ba","#c893c9",
                               "#ff69b4","#d4a910","#c4625d","#bc3c28","#815375",
                               "#0072b5", "#1f854e","#e18726"))
ggsave("transcripts_PCA_all_min5in3_novel_flair.pdf", plot = last_plot(),
       width = 7.8, height = 5.7, units = "in", useDingbats=FALSE)

ont_tpm_df_filter_cor2 <- cor(ont_nano_tpm_uniq_filter_novel, method = "spearman")
sampleDists <- dist(ont_tpm_df_filter_cor2, method = "euclidean")
sampleDistMatrix <- as.matrix(ont_tpm_df_filter_cor2)
my_sample_col <- as.data.frame(samps[samps$name %in% colnames(ont_tpm_df_filter_cor2),c("group")])
row.names(my_sample_col) <- samps[samps$name %in% colnames(ont_tpm_df_filter_cor2),]$name
colnames(my_sample_col) <- c("tissue")

pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists, show_rownames = FALSE,
         clustering_distance_cols=sampleDists, show_colnames = FALSE,
         annotation_colors = my_colour,
         col=colorRampPalette(brewer.pal(9, "BrBG"))(length(breaksList)),
         annotation_col = my_sample_col, drop_levels = TRUE,
         breaks=breaksList)

#Analyse transcripts at different thresholds
transcript_metadata_flair <- as.data.frame(rownames(ont_nano_tpm_uniq))
transcript_metadata_flair$transcript_id <- gff_class$transcript_id[match(transcript_metadata_flair$`rownames(ont_nano_tpm_uniq)`, gff_class$fish)]
transcript_metadata_flair$gene_id <- gff_class$gene_id[match(transcript_metadata_flair$`rownames(ont_nano_tpm_uniq)`, gff_class$fish)]

ont_nano_tpm_uniq_a <- ont_nano_tpm_uniq[,-1]
ont_nano_tpm_uniq_a$transcript_id <- transcript_metadata_flair$transcript_id
ont_nano_tpm_uniq_a$gene_id <- transcript_metadata_flair$gene_id

ont_nano_tpm_uniq_m <- reshape2::melt(ont_nano_tpm_uniq_a)
colnames(ont_nano_tpm_uniq_m) <- c("transcript_id","gene_id","sample_id","TPM")
ont_nano_tpm_uniq_m$Tissue <- samps$group[match(ont_nano_tpm_uniq_m$sample_id, samps$name)]
ont_nano_tpm_uniq_m <- ont_nano_tpm_uniq_m[ont_nano_tpm_uniq_m$Tissue %in% c("Muscle_Skeletal","Lung","Liver","Heart_Left_Ventricle","Heart_Atrial_Appendage","Brain_Cerebellar_Hemisphere",
                                                                             "Brain_Frontal_Cortex_BA9","Brain_Putamen_basal_ganglia","Cells_Cultured_fibroblasts"),]

gff_class <- read_tsv("analysis/flair/new_noSJ_promoters.isoforms_annotation.txt")
gff_class_novel <- gff_class[gff_class$Class_code != "=",]
gff_class_annot <- gff_class[gff_class$Class_code == "=",]

transcripts_per_tissue <- function(df, TPM_thresh, isNovel) {
  transcripts_per_tissue_thresh <- df %>%
    filter(TPM > TPM_thresh) %>%
    filter(transcript_id %in% gff_class[gff_class$Annotation == isNovel,]$transcript_id) %>%
    group_by(transcript_id, Tissue) %>%
    mutate(freq = n()) %>% 
    ungroup() %>% 
    filter(freq > 1) %>%
    select(-freq) %>%
    select(-sample_id, -TPM, -gene_id) %>%
    distinct() %>%
    group_by(Tissue) %>%
    group_split(.keep = FALSE)
  names(transcripts_per_tissue_thresh) <- sort(unique(df$Tissue))
  listInput <-  lapply(transcripts_per_tissue_thresh, function(x) as.vector(x[['transcript_id']]))
  bin_table <- fromList(listInput)
  bin_table$sum <- rowSums(bin_table)
  bin_table$TPM_filter <- TPM_thresh
  bin_table$Transcript <- isNovel
  return(bin_table)
}

bin_table_novel_01 <- transcripts_per_tissue(ont_nano_tpm_uniq_m, 0.1, "Novel")
bin_table_novel_1 <- transcripts_per_tissue(ont_nano_tpm_uniq_m, 1, "Novel")
bin_table_novel_10 <- transcripts_per_tissue(ont_nano_tpm_uniq_m, 10, "Novel")

bin_table_annot_01 <- transcripts_per_tissue(ont_nano_tpm_uniq_m, 0.1, "Annotated")
bin_table_annot_1 <- transcripts_per_tissue(ont_nano_tpm_uniq_m, 1, "Annotated")
bin_table_annot_10 <- transcripts_per_tissue(ont_nano_tpm_uniq_m, 10, "Annotated")

bin_table <- rbind(bin_table_annot_01, bin_table_novel_01,
                   bin_table_annot_1, bin_table_novel_1,
                   bin_table_annot_10, bin_table_novel_10)
bin_table$sum <- as.factor(bin_table$sum)
bin_table$TPM <- as.factor(bin_table$TPM)

bin_table %<>%
  group_by(Transcript, TPM_filter, sum) %>%
  tally() %>%
  ungroup() %>%
  group_by(Transcript, TPM_filter) %>%
  mutate(Prop = n/sum(n), Sum = sum(n))

ggplot(bin_table, aes(x=Transcript, y=Prop, fill=Transcript, alpha=TPM_filter)) +
  geom_bar(stat="identity", position="dodge") +
  facet_grid(~sum, space = "free") +
  theme_classic(base_size=14) +
  xlab("Number of tissues") +
  ylab("Proportion of transcripts") +
  scale_alpha_discrete(c(0.2, 0.6, 1)) +
  scale_fill_manual(values = c("#54ba5d","#c0504d")) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave("at_least_two_tpm_filter_perc.pdf", plot = last_plot(),
       width = 7, height = 4, units = "in", useDingbats=FALSE)
