library(GenomicFeatures)
library(tidyverse)
library(magrittr)
library(gridExtra)
library(data.table)
library(annotables)

#Generates Figure 1C and Extended Figure 2

# Load stranscript to gene naming
txdb = makeTxDbFromGFF('data/gencode.v26.annotation.gtf')
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
tx2gene %<>%
  dplyr::group_by(GENEID) %>%
  dplyr::mutate(n_tx = n())
txlengths <- GenomicFeatures::transcriptLengths(txdb)

# Load samples text file and subset to files also sequenced in the GTEx
samps <- read.table("data/metadata.txt",
                    header=TRUE, fill = TRUE, sep = "\t")
samps <- samps[!grepl("experimental", samps$old_name),]
samps <- samps[!grepl("direct", samps$old_name),]
samps <- samps[!grepl("CVD", samps$old_name),]
samps <- samps[!grepl("K562", samps$old_name),]
#too few features
samps <- samps[!samps$path %in% c("FAK91589"),]
#remove replicates 
samps <- samps[!grepl("rep",samps$sample_id),]
# Do some renaming to match GTEx names
samps$ID <- paste0(samps$sample_name,"_",samps$tissue)
samps$ID <- gsub(" - ","_",samps$ID)
samps$ID <- gsub(" \\(","_",samps$ID)
samps$ID <- gsub("\\)","",samps$ID)
samps$ID <- gsub(" ","_",samps$ID)
samps$ID <- gsub("-","_",samps$ID)

# Load GTEx naming file and subset to files also sequenced in the ONT
gtex_char <- read_tsv("data/GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.txt", col_names=TRUE)
gtex_char %<>% 
  filter(SMTSD %in% samps$tissue) %>%
  separate(SAMPID, into=c('GTEX','donor','other'), sep="-", remove=FALSE) %>%
  mutate(ID = paste0(GTEX,"_",donor,"_",SMTSD)) %>%
  mutate(ID = gsub(" - ","_",ID)) %>%
  mutate(ID = gsub(" \\(","_",ID)) %>%
  mutate(ID = gsub("\\)","",ID)) %>%
  mutate(ID = gsub(" ","_",ID)) %>%
  filter(ID %in% samps$ID) %>%
  select(SAMPID, SMRIN, SMTSD, SMNABTCHT, ID)

##### Load ONT abundance files
ont_nano_counts <- read.table("analysis/count_tables/trans_gencode_reps.counts.txt", header=TRUE)
ont_nano_counts$ids <- rownames(ont_nano_counts)
ont_nano_counts <- ont_nano_counts %>%
  pivot_longer(-ids, names_to = "sample_id", values_to = "counts")

ont_nano_counts$variable <- samps$ID[match(ont_nano_counts$sample_id, samps$old_name)]
ont_nano_counts$gene_id <- tx2gene$GENEID[match(ont_nano_counts$ids, tx2gene$TXNAME)]

ont_nano_counts_tpm <- ont_nano_counts %>%
  dplyr::group_by(sample_id, variable) %>%
  dplyr::mutate(tpm = counts / sum(counts) * 1e6) %>%
  dplyr::mutate(log_tpm = log10(tpm+1)) %>%
  dplyr::ungroup()

## SUM UP TO GENES
ont_nano_gene_counts <- ont_nano_counts %>%
  dplyr::group_by(gene_id, sample_id, variable) %>%
  dplyr::summarise(counts = sum(counts)) %>%
  dplyr::ungroup(gene_id) %>%
  dplyr::mutate(tpm = counts / sum(counts) * 1e6, log_tpm = log10(tpm+1)) %>%
  dplyr::ungroup()

##### READ IN ILLUMINA (can be downloaded from the GTEx portal)
#Genes
#genes_illumina <- read_tsv("data/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_gene_tpm.gct.gz",
#                             skip = 2, col_names = TRUE)
#genes_illumina1 <- genes_illumina[,c(1,2)]
#genes_illumina2 <- genes_illumina[,gtex_char$SAMPID]
#colnames(genes_illumina2) <- gtex_char$ID
#genes_illumina2_log <- log10(genes_illumina2 + 1)
#genes_illumina2_log <- cbind(genes_illumina1,genes_illumina2_log)
#write.table(genes_illumina2_log, "analysis/count_tables/illumina_genes_gencode.log10tpm.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")
genes_illumina2_log <- read.table("analysis/count_tables/illumina_genes_gencode.log10tpm.txt", header = TRUE, sep = "\t")
genes_illumina_melt <- reshape2::melt(genes_illumina2_log)

#Transcripts
#transcripts_illumina = fread("data/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz",
#                             select = c("transcript_id", "gene_id", gtex_char$SAMPID), skip = 2)
#transcripts_illumina <- as.data.frame(transcripts_illumina)
#transcripts_illumina1 <- transcripts_illumina[,c(1,2)]
#transcripts_illumina2 <- transcripts_illumina[,gtex_char$SAMPID]
#colnames(transcripts_illumina2) <- gtex_char$ID
#transcripts_illumina2_log <- log10(transcripts_illumina2 + 1)
#transcripts_illumina2_log <- cbind(transcripts_illumina1,transcripts_illumina2_log)
#write.table(transcripts_illumina2_log, "analysis/count_tables/illumina_transcripts_gencode.log10tpm.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")
transcripts_illumina2_log <- read.table("analysis/count_tables/illumina_transcripts_gencode.log10tpm.txt", header = TRUE, sep = "\t")
transcripts_illumina_melt <- reshape2::melt(transcripts_illumina2_log)

# COMPARE GENE COUNTS BETWEEN ONT AND ILLUMINA
tab_gene <- merge(genes_illumina_melt,ont_nano_gene_counts,by=c("gene_id","variable"))
threshold = threshold
tab_gene_sub = tab_gene[tab_gene$value > threshold & tab_gene$log_tpm > threshold, ]
tab_gene_sub <- tab_gene_sub[,c(1,4,5,8)]
tab_gene_sub <- droplevels(tab_gene_sub)
eq_gene <- plyr::ddply(tab_gene_sub,"sample_id", lm_eqn_table)
eq_gene$Category <- "Gene"

g1 <- ggplot(tab_gene_sub[tab_gene_sub$sample_id=="GTEX_1LVA9_muscle",], aes(x=value, y=log_tpm)) +
  stat_binhex(bins = 80) +
  #geom_point(alpha=0.5, color="grey20") +
  geom_abline(slope = 1, color="#1e854f", size=2, alpha=0.5) +
  #geom_smooth(method = "lm", se=FALSE, color="#1e854f", formula = y ~ x) +
  theme_classic(base_size=14) +
  scale_fill_gradient(low = "grey80", high = "black") +
  xlab("Illumina LOG10 gene TPM") +
  ylab("ONT LOG10 gene TPM") +
  ylim(0,4.8) + xlim(0,4.8)

tab_gene_sub_list <- split( tab_gene_sub , f = tab_gene_sub$sample_id )
names(tab_gene_sub_list) <- unique(tab_gene_sub$sample_id)

for (i in 1:length(tab_gene_sub_list)){
  df = tab_gene_sub_list[[i]]
  row.names(df) <- df$gene_id
  df <- df[,-c(1,3)]
  tab_gene_sub_list[[i]] <- df
}
genes.lm <- lapply(tab_gene_sub_list, FUN = fitLM)

extreme_genes <- list()
for (i in 1:length(genes.lm)){
  myresid = resid(genes.lm[[i]])
  names(myresid) <- rownames(tab_gene_sub_list[[i]])
  myresid <- myresid[abs(myresid) > 1]
  extreme_genes[[names(tab_gene_sub_list)[[i]]]] <- as.data.frame(myresid)
}

extreme_genes_df <- do.call(rbind.data.frame, extreme_genes)
extreme_genes_df$names <- rownames(extreme_genes_df)
extreme_genes_df %<>%
  separate(names, c("sample_id", "gene_id"), sep = "\\.ENSG", extra = "merge", remove = FALSE) %>%
  mutate(gene_id = paste0("ENSG",gene_id))

### ANALYZE LOG10 EXPRESSION
tab_gene$names <- paste0(tab_gene$sample_id, ".", tab_gene$gene_id)

extreme_genes_df2 <- merge(extreme_genes_df, tab_gene, by=c("sample_id","gene_id","names"))
extreme_genes_df2 <- extreme_genes_df2[,c(1:3,7,10)]
extreme_genes_df_expr <- reshape2::melt(extreme_genes_df2, measure_ids=c("value","log_tpm"))

p1 <- ggplot(extreme_genes_df_expr, aes(x=value, color=variable)) +
  geom_density(size=2) +
  theme_classic(base_size=14) +
  xlab("LOG10 gene expression") +
  theme(legend.position = "bottom") +
  scale_color_manual(values=c("#ffb441","#0b5971"))

### ANALYZE LENGTH
extreme_genes_df2$Category <- "Extreme"
nonextreme_genes_df2 <- tab_gene[!tab_gene$names %in% extreme_genes_df2$names,]
nonextreme_genes_df2 <- nonextreme_genes_df2[,c(5,1,9,4,8)]
nonextreme_genes_df2$Category <- "Non-extreme"

genes_df2 <- rbind(extreme_genes_df2, nonextreme_genes_df2)
grch38$length <- grch38$end-grch38$start
genes_df2$length <- grch38$length[match(gsub("\\..*","",genes_df2$gene_id), grch38$ensgene)]

t1 <- wilcox.test(genes_df2[genes_df2$Category=="Extreme",]$length, genes_df2[genes_df2$Category=="Non-extreme",]$length)
p2 <- ggplot(genes_df2, aes(x=log2(length), color=Category)) + 
  geom_density(size=2, alpha=0.5) +
  theme_classic(base_size=14) +
  xlab("LOG2 Length") +
  theme(legend.position = "bottom") +
  scale_color_manual(values=c("#c28834","#45a19a")) +
  annotate(geom = "text", x=10, y=0.2, label=paste0("p-value = ", t1$p.value))

genes_df2$n_tx <- tx2gene$n_tx[match(genes_df2$gene_id, tx2gene$GENEID)]
t2 <- wilcox.test(genes_df2[genes_df2$Category=="Extreme",]$n_tx, genes_df2[genes_df2$Category=="Non-extreme",]$n_tx)
genes_df2_sum <- genes_df2 %>%
  dplyr::mutate(n_tx = ifelse(n_tx > 30, 30, n_tx)) %>%
  dplyr::group_by(n_tx, Category) %>%
  dplyr::summarize(Freq = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(Category) %>%
  dplyr::mutate(Perc=Freq/sum(Freq)) %>%
  dplyr::ungroup()

p3 <- ggplot() + 
  geom_bar(data=genes_df2_sum[genes_df2_sum$Category=="Extreme",], aes(x=n_tx, y=Perc),
           fill="#c28834", alpha=0.5, stat = "identity", position="dodge") +
  geom_bar(data=genes_df2_sum[genes_df2_sum$Category!="Extreme",], aes(x=n_tx, y=Perc), 
           fill="#45a19a", alpha=0.5, stat = "identity", position="dodge") +
  theme_classic(base_size=14) +
  xlab("Number of transcripts") +
  theme(legend.position = "bottom")

df_genes <- as.data.frame(table(extreme_genes_df2$gene_id)[table(extreme_genes_df2$gene_id) > 10])
df_genes$symbol <- grch38$symbol[match(gsub("\\..*","",df_genes$Var1), grch38$ensgene)]

# COMPARE TRANSCRIPT COUNTS BETWEEN ONT AND ILLUMINA
tab_transcript <- merge(transcripts_illumina_melt,ont_nano_counts_tpm,
                        by.x=c("transcript_id","variable","gene_id"), by.y=c("ids","variable","gene_id"))

tab_transcript_sub <- tab_transcript[tab_transcript$value > 0.0000001 & tab_transcript$log_tpm > 0.0000001,]
tab_transcript_sub <- tab_transcript_sub[,c(1,4,5,8)]
tab_transcript_sub <- droplevels(tab_transcript_sub)
eq_transcript <- plyr::ddply(tab_transcript_sub,"sample_id", lm_eqn_table)
eq_transcript$Category <- "Transcript"

g2 <- ggplot(tab_transcript_sub[tab_transcript_sub$sample_id=="GTEX_1LVA9_muscle",], aes(x=value, y=log_tpm)) +
  stat_binhex(bins = 80) +
  #geom_point(alpha=0.5, color="grey20") +
  geom_abline(slope = 1, color="#1e854f", size=2, alpha=0.5) +
  #geom_smooth(method = "lm", se=FALSE, color="#1e854f", formula = y ~ x) +
  theme_classic(base_size=14) +
  scale_fill_gradient(low = "grey80", high = "black") +
  xlab("Illumina LOG10 transcript TPM") +
  ylab("ONT LOG10 transcript TPM") +
  ylim(0,4.8) + xlim(0,4.8)

tab_transcript_sub_list <- split( tab_transcript_sub , f = tab_transcript_sub$sample_id )
names(tab_transcript_sub_list) <- unique(tab_transcript_sub$sample_id)

for (i in 1:length(tab_transcript_sub_list)){
  df = tab_transcript_sub_list[[i]]
  row.names(df) <- df$transcript_id
  df <- df[,-c(1,3)]
  tab_transcript_sub_list[[i]] <- df
}

transcript.lm <- lapply(tab_transcript_sub_list, FUN = fitLM)

extreme_transcripts <- list()
for (i in 1:length(transcript.lm)){
  myresid = resid(transcript.lm[[i]])
  names(myresid) <- rownames(tab_transcript_sub_list[[i]])
  myresid <- myresid[abs(myresid) > 1]
  extreme_transcripts[[names(tab_transcript_sub_list)[[i]]]] <- as.data.frame(myresid)
}

extreme_transcripts_df <- do.call(rbind.data.frame, extreme_transcripts)
extreme_transcripts_df$names <- rownames(extreme_transcripts_df)
extreme_transcripts_df %<>%
  separate(names, c("sample_id", "transcript_id"), sep = "\\.ENST", extra = "merge", remove = FALSE) %>%
  mutate(transcript_id = paste0("ENST",transcript_id))
  
### ANALYZE LOG10 EXPRESSION
tab_transcript$names <- paste0(tab_transcript$sample_id, ".", tab_transcript$transcript_id)
extreme_transcripts_df2 <- merge(extreme_transcripts_df, tab_transcript, by=c("sample_id","transcript_id","names"))
extreme_transcripts_df2 <- extreme_transcripts_df2[,c(1:3,7,10)]
extreme_transcripts_df_expr <- reshape2::melt(extreme_transcripts_df2, measure_ids=c("value","log_tpm"))

p4 <- ggplot(extreme_transcripts_df_expr, aes(x=value, color=variable)) +
  #geom_histogram(position="dodge", bins = 20) +
  geom_density(size=2) +
  theme_classic(base_size=14) +
  xlab("LOG10 transcript expression") +
  theme(legend.position = "bottom") +
  scale_color_manual(values=c("#ffb441","#0b5971"))

### ANALYZE LENGTH
extreme_transcripts_df2$Category <- "Extreme"
nonextreme_transcripts_df2 <- tab_transcript[!tab_transcript$names %in% extreme_genes_df2$names,]
nonextreme_transcripts_df2 <- nonextreme_transcripts_df2[,c(5,1,9,4,8)]
nonextreme_transcripts_df2$Category <- "Non-extreme"

transcripts_df2 <- rbind(extreme_transcripts_df2, nonextreme_transcripts_df2)
transcripts_df2$length <- txlengths$tx_len[match(transcripts_df2$transcript_id, txlengths$tx_name)]
transcripts_df2$nexon <- txlengths$nexon[match(transcripts_df2$transcript_id, txlengths$tx_name)]

t3 <- wilcox.test(log2(transcripts_df2[transcripts_df2$Category=="Extreme",]$length), log2(transcripts_df2[transcripts_df2$Category=="Non-extreme",]$length))
p5 <- ggplot(transcripts_df2, aes(x=log2(length), color=Category)) + 
  geom_density(size=2) +
  theme_classic(base_size=14) +
  xlab("LOG2 Length") +
  theme(legend.position = "bottom") +
  scale_color_manual(values=c("#c28834","#45a19a")) +
  annotate(geom = "text", x=10, y=0.2, label=paste0("p-value = ", t3$p.value))

t4 <- wilcox.test(transcripts_df2[transcripts_df2$Category=="Extreme",]$nexon, transcripts_df2[transcripts_df2$Category=="Non-extreme",]$nexon)
transcripts_df2$nexon <- ifelse(transcripts_df2$nexon > 40, 40, transcripts_df2$nexon)
transcripts_df2_sum <- transcripts_df2 %>%
  dplyr::mutate(nexon = ifelse(nexon > 30, 30, nexon)) %>%
  dplyr::group_by(nexon, Category) %>%
  dplyr::summarize(Freq = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(Category) %>%
  dplyr::mutate(Perc=Freq/sum(Freq)) %>%
  dplyr::ungroup()
p6 <- ggplot() + 
  geom_bar(data=transcripts_df2_sum[transcripts_df2_sum$Category=="Extreme",], aes(x=nexon, y=Perc),
           fill="#c28834", alpha=0.5, stat = "identity", position="dodge") +
  geom_bar(data=transcripts_df2_sum[transcripts_df2_sum$Category!="Extreme",], aes(x=nexon, y=Perc), 
           fill="#45a19a", alpha=0.5, stat = "identity", position="dodge") +
  theme_classic(base_size=14) +
  xlab("Number of exons")

pdf(file = "correlation_analysis_expression.pdf", width =6, height = 10)
multiplot(p1, p4, cols=1)
dev.off()
pdf(file = "correlation_analysis_overview.pdf", width =10, height = 10)
multiplot(p2, p5, p3, p6, cols=2)
dev.off()
pdf(file = "expression_correlation_example.pdf", width =10, height = 6)
multiplot(g1,g2, cols=2)
dev.off()

### PLOT CORRELATIONS
eq_table <- rbind(eq_gene, eq_transcript)
eq_table <- na.omit(eq_table)
eq_table$V1 <- as.numeric(as.character(eq_table$V1))

ggplot(eq_table, aes(x=Category, y=V1)) +
  geom_jitter(size=4, colour="grey20",pch=21, fill="#c4c4c4") +
  theme_classic(base_size=16) +
  theme(panel.grid.major.y = element_line(color="grey70")) +
  ylab(expression(paste(rho," between Illumina and ONT"))) + xlab("") +
  ylim(0,1) + guides(fill = guide_legend(override.aes = list(size=6)))
ggsave("cor_rho_ont_illumina_flair_style_updated.pdf", plot = last_plot(),
       width = 4, height = 4, units = "in", useDingbats=FALSE)

