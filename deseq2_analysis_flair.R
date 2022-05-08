library(tidyverse)
library(DESeq2)
library(GenomicFeatures)
library(annotables)

gff_class <- read_tsv("analysis/flair/new_noSJ_promoters.isoforms_annotation.txt")

### DESEq2 on Illumina
# Load illumina reads
counts_ill <- read.table("analysis/ptbp1_kd/illumina_kd_count_table.txt", header=TRUE)
rownames(counts_ill) <- counts_ill$Geneid
counts_ill <- counts_ill[,-c(1:6)]

counts_ill <- counts_ill [ rowSums(counts_ill[,c(1,3,5,7,9)] > 10) >= 2 | rowSums(counts_ill[,c(2,4,6,8,10)] > 10) >= 2, ]

sample_info <- read.table("analysis/ptbp1_kd/metadata_illumina.txt", header=TRUE, sep="\t")
rownames(sample_info) <- sample_info$Name
sample_info <- sample_info[,-1]
rownames(sample_info) <- gsub("-",".",rownames(sample_info))

dds <- DESeqDataSetFromMatrix(countData = counts_ill,
                              colData = sample_info,
                              design = ~ Sample + Treatment)
dds <- DESeq(dds)

vsd <- vst(dds, blind=FALSE)
res <- DESeq2::results(dds, alpha=0.05)
res_sig <- na.omit(res)
res_sig <- res_sig[res_sig$padj <= 0.05,]
res_sig <- res_sig[order(res_sig$padj),] 

### DESEQ2 on ONT
# Load ONT reads
# Prepare metadata table
samps <- read_tsv(file = "data/metadata.txt", col_names = TRUE)

#Transform replicates info to make it compatible with DESeq2
samps <- samps %>%
  filter(date_of_sequencing == "11218") %>%
  dplyr::select(-RNA_extraction_method, -WGS, -data_center, -amount_loaded_ng, -apprx_run_time, -flush_buffer,
                -date_of_sequencing, -tissue, -protocol) %>%
  separate(old_name, into = c("sample_name","condition","rep"), sep = "_", remove = FALSE) %>%
  mutate(new_name = paste0(sample_name,"_",condition))

samps$rep <- ifelse(is.na(samps$rep), 1, 2)
samps <- as.data.frame(samps)
rownames(samps) <- samps$old_name

counts_ont <- read.table("analysis/trans_flair.counts.txt",
                         header=TRUE, sep="\t")
counts_ont <- counts_ont[,grepl("control", colnames(counts_ont)) | grepl("experimental", colnames(counts_ont)) | grepl("transcript_id", colnames(counts_ont))]

genes_counts_ont <- counts_ont
genes_counts_ont$gene_id <- gff_class$gene_id[match(genes_counts_ont$transcript_id, gff_class$fish)]

genes_counts_ont <- aggregate(.~gene_id, genes_counts_ont[,-1], sum)
rownames(genes_counts_ont) <- genes_counts_ont$gene_id
genes_counts_ont <- genes_counts_ont[,-1]

genes_counts_ont <- genes_counts_ont [ rowSums(genes_counts_ont[,grepl("control",colnames(genes_counts_ont))] > 10) >= 2 |
                             rowSums(genes_counts_ont[,grepl("experimental",colnames(genes_counts_ont))] > 10) >= 2, ]
genes_counts_ont = genes_counts_ont[ , rownames(samps) ]

dds_ont <- DESeqDataSetFromMatrix(countData = genes_counts_ont,
                                  colData = samps,
                                  design = ~ sample_name + condition)
dds_ont <- DESeq(dds_ont)

#Collapse replicates
dds_ont <- collapseReplicates(dds_ont, dds_ont$new_name, dds_ont$rep)

vsd_ont <- vst(dds_ont, blind=FALSE)
res_ont <- DESeq2::results(dds_ont, alpha=0.05)

res_ont$symbol <- grch38$symbol[match(gsub("\\..*","",rownames(res_ont)),grch38$ensgene)]

df <- as.data.frame(res_ont)

df$illumina = ifelse(rownames(df) %in% rownames(res_sig), "illumina","ont only")

df$threshold <- ifelse(abs(df$log2FoldChange) >= 1 & df$padj <= 0.05, "p-value and FC",
                       ifelse(abs(df$log2FoldChange) >= 1 & df$padj > 0.05, "FC", 
                              ifelse(abs(df$log2FoldChange) < 1 & df$padj <= 0.05, "p-value","not significant")))

ggplot(data=df, aes(x=log2FoldChange, y=-log10(padj), label=symbol, shape=illumina)) +
  geom_point(alpha=0.8, size=2, aes(color = threshold)) +
  scale_color_manual(values=c("#666666","#c1c1c1","#666666","#d95f02")) +
  theme_bw(base_size=16) +
  xlim(-5.5,5.5) +
  theme(legend.title=element_blank(),panel.grid.minor = element_blank(),legend.position = "bottom") +
  xlab(expression(paste("LOG"[2]," fold-change"))) + ylab(expression(paste("LOG"[10]," p-value")))
ggsave("volcano_plot.pdf", height = 5, width = 4)

d1 <- plotCounts(dds_ont, gene="ENSG00000117569.18", intgroup=c("condition", "sample_name"), returnData=TRUE)
d1$symbol <- "PTBP2"
d2 <- plotCounts(dds_ont, gene="ENSG00000011304.19", intgroup=c("condition", "sample_name"), returnData=TRUE)
d2$symbol <- "PTBP1"
d <- rbind(d1,d2)

ggplot(d, aes(x=condition, y=count, color=sample_name)) + 
  geom_point(position=position_jitter(w=0.1,h=0), size=5, alpha=0.5) +
  theme_bw(base_size=16) +
  facet_grid(.~symbol) +
  ylab("Normalized read count") + xlab("Condition") +
  scale_y_log10() +
  scale_color_brewer(palette="Dark2")
ggsave("ptbp1_kd_counts.pdf", height = 5, width = 4)

## DIFFERENTIAL TRANSCRIPT EXPRESSION
counts_ont <- as.data.frame(counts_ont)
rownames(counts_ont) <- counts_ont$transcript_id
counts_ont <- counts_ont[,-1]

counts_ont <- counts_ont [ rowSums(counts_ont[,grepl("control",colnames(counts_ont))] > 10) >= 2 |
                             rowSums(counts_ont[,grepl("experimental",colnames(counts_ont))] > 10) >= 2, ]
counts_ont = counts_ont[ , rownames(samps) ]

dds_ont <- DESeqDataSetFromMatrix(countData = counts_ont,
                                  colData = samps,
                                  design = ~ sample_name + condition)
dds_ont <- DESeq(dds_ont)

#Collapse replicates
dds_ont <- collapseReplicates(dds_ont, dds_ont$new_name, dds_ont$rep)

vsd_ont <- vst(dds_ont, blind=FALSE)
res_ont <- DESeq2::results(dds_ont, alpha=0.05)
res_ont <- as.data.frame(res_ont)
res_ont$transcript_id <- rownames(res_ont)

### ANALYSE ASTS EVENTS by AS CATEGORY
events <- read.table("analysis/diffsplice/new_noSJ_promoters_all_events.tsv", header=TRUE)
events$diffsplice_events <- gsub("AF5","A5",events$diffsplice_events)
events$diffsplice_events <- gsub("AL3","A3",events$diffsplice_events)
events$fish <- gff_class$fish[match(events$transcript_id, gff_class$transcript_id)]
events <- events[events$fish %in% rownames(res_ont),]
tracking_char_events <- merge(events, res_ont, by.x="fish",by.y="transcript_id")
tracking_char_events <- tracking_char_events[tracking_char_events$Type=="inclusion",]

tracking_char_events_sum_sig <- tracking_char_events %>%
  filter(padj <= 0.05) %>%
  mutate(direction=ifelse(log2FoldChange < 0,"negative","positive")) %>%
  dplyr::select(transcript_id, event_id, direction, gene_id, diffsplice_events) %>%
  unique() %>%
  group_by(diffsplice_events, direction) %>%
  summarise(No_event=n()) %>%
  ungroup() %>%
  group_by(direction) %>%
  mutate(Proportion = (No_event/sum(No_event))*100) %>%
  ungroup()

ggplot(tracking_char_events_sum_sig, aes(x=diffsplice_events, fill=direction, group=direction, y=Proportion)) +
  geom_bar(stat="identity", position='dodge') +
  theme_classic(base_size=14) +
  xlab("") +
  ylab("Proportion of AltTS events") +
  theme(legend.position = "bottom") +
  scale_fill_manual(values=c("#c1997b","#25857d")) +
  coord_flip()
ggsave("altts_ptbp1_kd.pdf", height = 5, width = 4)

pvalues_prop <- tracking_char_events_sum_sig %>%
  group_by(direction) %>%
  mutate(No_all = sum(No_event)) %>%
  ungroup() %>%
  dplyr::select(-direction, -Proportion) %>%
  dplyr::nest_by(diffsplice_events) %>%
  dplyr::mutate(value = list(prop.test(data$No_event, data$No_all)$p.value)) %>%
  ungroup() %>%
  dplyr::select(diffsplice_events, value) %>%
  unnest(value) %>%
  mutate(qvalue = p.adjust(value, method = "fdr"))
