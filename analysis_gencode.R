library(tidyverse)
library(magrittr)
library(pheatmap)
library(RColorBrewer)

#Generates extended Figure 1D,E and Figure 1A,B

#Pheatmap parameters
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

samps <- read.table("data/metadata.txt",
                    header=TRUE, fill = TRUE, sep = "\t")
samps <- samps[samps$path != "FAK91589",]
samps <- samps[!grepl("CVD",samps$name),]
samps$group <- gsub(" - ","_",samps$tissue)
samps$group <- gsub(" \\(","_",samps$group)
samps$group <- gsub("\\)","",samps$group)
samps$group <- gsub(" ","_",samps$group)
samps$sample_id <- gsub("-",".",samps$sample_id)

ont_nano_tpm_df <- read.table("analysis/counts_tables/trans_gencode.tpm.txt", header = TRUE)
rownames(ont_nano_tpm_df) <- ont_nano_tpm_df$sample
ont_nano_tpm_df <- ont_nano_tpm_df[,-1]

ont_nano_tpm_df <- ont_nano_tpm_df[,colnames(ont_nano_tpm_df) %in% samps$name]
ont_nano_tpm_df <- ont_nano_tpm_df[rowSums(ont_nano_tpm_df) > 0, ]
ont_nano_tpm_df <- na.omit(ont_nano_tpm_df)

### PCA of Transcripts
ont_nano_tpm_df_filter <- ont_nano_tpm_df[ rowSums(ont_nano_tpm_df >= 5) >= 3, ]
highly_variable_lcpm <- log10(ont_nano_tpm_df_filter+1)
genes.pca <- prcomp(t(highly_variable_lcpm), center = TRUE, scale. = TRUE)
percentVar <- round(100 * apply(genes.pca$x, 2, var)/sum(apply(genes.pca$x, 2, var)))
pcs <- as.data.frame(genes.pca$x[, c(1:6)])
pcs$name <- rownames(pcs)
rownames(pcs) <- NULL

pcs <- merge(pcs, samps, by = "name")
ggplot(pcs, aes(PC1, PC2)) +
  geom_point(aes(fill=tissue),colour="grey20",pch=21, size=5) +
  theme_classic(base_size =14) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + guides(fill = guide_legend(override.aes = list(size=6))) +
  scale_fill_manual(values = c("#ffc0cb","#8d5b96","#7776b1","#9773ba","#b873ba","#c893c9",
                               "#ff69b4","#d4a910","#c4625d","#bc3c28","#B09638","#815375",
                               "#0072b5", "#1f854e","#e18726"))
ggsave("transcripts_PCA_all_min5in3_ALL_gencode.pdf", plot = last_plot(),
       width = 7.8, height = 5.7, units = "in", useDingbats=FALSE)

### REPLICATE ANALYSIS
ont_nano_tpm_df_reps <- ont_nano_tpm_df[,colnames(ont_nano_tpm_df) %in% c("QV44_control","QV44_control_rep",
                                                                          "S4Z8_experimental","S4Z8_experimental_rep",
                                                                          "S95S_experimental","S95S_experimental_rep",
                                                                          "GTEX_1GN1W_heartatrialappendage","GTEX_1GN1W_heartatrialappendage_rep",
                                                                          "GTEX_13QJ3_muscle","GTEX_13QJ3_muscle_rep",
                                                                          "GTEX_SM3QNG2_fibs","GTEX_SM3QNG2_fibs_rep",
                                                                          "GTEX_Y5LM_liver","GTEX_Y5LM_liver_rep","GTEX_Y5LM_liver_rep2",
                                                                          "GTEX_14BMU_lung","GTEX_14BMU_lung_rep","GTEX_14BMU_lung_rep2",
                                                                          "GTEX_ZT9X_muscle","GTEX_ZT9X_muscle_rep","GTEX_ZT9X_muscle_rep2",
                                                                          "GTEX_WY7C_heartatrialappendage","GTEX_WY7C_heartatrialappendage_rep","GTEX_WY7C_heartatrialappendage_rep2",
                                                                          "GTEX_Q2AG_braincerebellarhemisphere","GTEX_Q2AG_braincerebellarhemisphere_rep","GTEX_Q2AG_braincerebellarhemisphere_rep2")]
ont_nano_tpm_df_filter <- ont_nano_tpm_df_reps[ rowSums(ont_nano_tpm_df_reps >= 1) >= 5, ]
ont_nano_counts_df_filter_cor <- cor(ont_nano_tpm_df_filter, method = "spearman")
sampleDists_reps <- dist(ont_nano_counts_df_filter_cor, method = "euclidean")
sampleDistMatrix_reps <- as.matrix(ont_nano_counts_df_filter_cor)
my_sample_col_reps <- as.data.frame(samps[samps$name %in% colnames(ont_nano_counts_df_filter_cor),c("name","tissue")])
row.names(my_sample_col_reps) <- my_sample_col_reps$name
my_sample_col_reps <- my_sample_col_reps[,-1, drop=F]

my_heatmap <- pheatmap(sampleDistMatrix_reps,
                       clustering_distance_rows=sampleDists_reps, show_rownames = FALSE,
                       clustering_distance_cols=sampleDists_reps, show_colnames = FALSE,
                       col=colorRampPalette(brewer.pal(9, "BrBG"))(length(breaksList)),
                       annotation_col = my_sample_col_reps, drop_levels = TRUE,
                       breaks=breaksList)
save_pheatmap_pdf(my_heatmap, "spearman_heatmap_transcripts_min5in3_SEL.pdf")

### REPLACE COLUMNS WITH THE MERGED COUNTS
merged_file_names <- c("GTEX_WY7C_heartatrialappendage","GTEX_13QJ3_muscle","GTEX_14BMU_lung","GTEX_ZT9X_muscle",
                       "GTEX_Q2AG_braincerebellarhemisphere","GTEX_Y5LM_liver","GTEX_1GN1W_heartatrialappendage",
                       "GTEX_SM3QNG2_fibs","QV44_control")
merged_files <- file.path("analysis/count_tables/",paste0(merged_file_names,"_gencode_tpm.txt"))
names(merged_files) <- merged_file_names

merged_ont_nano_tpm = lapply(merged_files, read_tsv, col_names=TRUE)
merged_ont_nano_tpm_df <- bind_rows(merged_ont_nano_tpm, .id = "names")

## merge EVERYTHING and turn into dataframe
merged_ont_nano_tpm_df %<>%
  spread(names, counts)

merged_ont_nano_tpm_df[is.na(merged_ont_nano_tpm_df)] <- 0
colnames(merged_ont_nano_tpm_df) <- samps$name[match(colnames(merged_ont_nano_tpm_df), samps$name)]
colnames(merged_ont_nano_tpm_df)[1] <- "rn"

toMatch <- c("rep", "K562", "direct", "experimental","GTEX_WY7C_heartatrialappendage","GTEX_13QJ3_muscle",
             "GTEX_14BMU_lung","GTEX_ZT9X_lung", "GTEX_Q2AG_braincerebellarhemisphere","QV44_control",
             "GTEX_Y5LM_liver","GTEX_1GN1W_heartatrialappendage", "GTEX_SM3QNG2_fibs")

matches <- unique (grep(paste(toMatch,collapse="|"), 
                        samps$name, value=TRUE))

ont_nano_tpm_uniq = ont_nano_tpm_df
for (i in matches) {
  ont_nano_tpm_uniq %<>%
    dplyr::select(-i)
}

ont_nano_tpm_uniq$rn <- rownames(ont_nano_tpm_uniq)

ont_nano_tpm_uniq <- plyr::join_all(list(ont_nano_tpm_uniq,merged_ont_nano_tpm_df), by = "rn", type = 'full')
rownames(ont_nano_tpm_uniq) <- ont_nano_tpm_uniq$rn 
ont_nano_tpm_uniq$rn <- NULL

ont_nano_tpm_uniq[is.na(ont_nano_tpm_uniq)] <- 0
ont_nano_tpm_uniq <- ont_nano_tpm_uniq[rowSums(ont_nano_tpm_uniq) > 0, ]
ont_nano_tpm_uniq <- na.omit(ont_nano_tpm_uniq)

write.table("analysis/count_tables/trans_gencode_reps.counts.txt", col.names = TRUE, row.names = TRUE,
            quote = FALSE, sep="\t")

# PCA of Transcripts for selected samples
ont_nano_tpm_uniq_filter <- ont_nano_tpm_uniq[ rowSums(ont_nano_tpm_uniq >= 5) >= 3, ]

highly_variable_lcpm <- log10(ont_nano_tpm_uniq_filter+1)
genes.pca <- prcomp(t(highly_variable_lcpm), center = TRUE, scale. = TRUE)
percentVar <- round(100 * apply(genes.pca$x, 2, var)/sum(apply(genes.pca$x, 2, var)),2)
pcs <- as.data.frame(genes.pca$x[, c(1:6)])
pcs$name <- rownames(pcs)
rownames(pcs) <- NULL

pcs <- merge(pcs, samps, by = "name")
ggplot(pcs, aes(PC1, PC2)) +
  geom_point(aes(fill=tissue),colour="grey20",pch=21, size=5) +
  theme_classic(base_size =14) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + guides(fill = guide_legend(override.aes = list(size=6))) +
  scale_fill_manual(values = c("#ffc0cb","#8d5b96","#7776b1","#9773ba","#b873ba","#c893c9",
                               "#ff69b4","#d4a910","#c4625d","#bc3c28","#815375",
                               "#0072b5", "#1f854e","#e18726"))
ggsave("transcripts_PCA_selected_min5in3_ALL_gencode.pdf", plot = p1,
       width = 7.8, height = 5.7, units = "in", useDingbats=FALSE)

samps <- samps[samps$name %in% colnames(ont_nano_tpm_uniq_filter),]
ont_nano_tpm_uniq_filter_cor <- cor(ont_nano_tpm_uniq_filter, method = "spearman")
sampleDists <- dist(ont_nano_tpm_uniq_filter_cor, method = "euclidean")
sampleDistMatrix <- as.matrix(ont_nano_tpm_uniq_filter_cor)
my_sample_col <- as.data.frame(samps[samps$name %in% colnames(ont_nano_tpm_uniq_filter_cor),c("group")])

row.names(my_sample_col) <- samps$name
colnames(my_sample_col) <- c("tissue")

breaksList = seq(0, 1, by = 0.01)

my_heatmap <- pheatmap(sampleDistMatrix,
                       clustering_distance_rows=sampleDists, show_rownames = FALSE,
                       clustering_distance_cols=sampleDists, show_colnames = FALSE,
                       annotation_colors = my_colour,
                       col=colorRampPalette(brewer.pal(9, "BrBG"))(length(breaksList)),
                       annotation_col = my_sample_col, drop_levels = TRUE,
                       breaks=breaksList)
save_pheatmap_pdf(my_heatmap, "spearman_heatmap_transcripts_min5in3_SEL.pdf")
