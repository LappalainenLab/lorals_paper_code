library(tidyverse)
library(GenomicFeatures)
library(magrittr)

#Process output of calculate_aligned_read_length_mit.sh
#Generates extended Figure 3B,C,D and Figure 1E

samps <- read.table("data/metadata.txt",
                    header=TRUE, fill = TRUE, sep = "\t")
samps <- samps[!grepl("CVD", samps$name),]
#samps <- samps[!grepl("experimental", samps$name),]
#too few features
samps <- samps[!samps$path %in% c("FAK91589"),]
samps$ID <- paste0(samps$sample_name,"_",samps$tissue)
samps$ID2 <- gsub(" - ","_",samps$ID)
samps$ID2 <- gsub(" \\(","_",samps$ID2)
samps$ID2 <- gsub("\\)","",samps$ID2)
samps$ID2 <- gsub(" ","_",samps$ID2)
samps <- droplevels(samps)

# Load GTEx naming file and subset to files also sequenced in the ONT
# (can be downloaded from here: https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt)
gtex_char <- read_tsv("data/GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.txt", col_names=TRUE)
gtex_char %<>% 
  filter(SMTSD %in% samps$tissue) %>%
  separate(SAMPID, into=c('GTEX','donor','other'), sep="-", remove=FALSE) %>%
  mutate(ID = paste0(GTEX,"-",donor,"_",SMTSD)) %>%
  filter(ID %in% samps$ID)
gtex_char$ID2 <- gsub(" - ","_",gtex_char$ID)
gtex_char$ID2 <- gsub(" \\(","_",gtex_char$ID2)
gtex_char$ID2 <- gsub("\\)","",gtex_char$ID2)
gtex_char$ID2 <- gsub(" ","_",gtex_char$ID2)

#Read in output of calculate_aligned_read_length_mit.sh
file_names <- unique(paste0("analysis/three_prime_bias/",samps$path,"_mt_coding_lengths.txt"))
my_lengths_ont <- lapply(file_names, read_delim, col_names=c("read_no","read_length","transcript_id"), delim=" ")
names(my_lengths_ont) <- samps$name

my_lengths_ont_df <- bind_rows(my_lengths_ont, .id = "sample_id")

#read in GTF file
gtf <- rtracklayer::import('data/gencode.v26.annotation.gtf')
gtf_df=as.data.frame(gtf)

gtf_df %<>%
  filter(type == "transcript") %>%
  filter(transcript_type == "protein_coding") %>%
  dplyr::select(-gene_type, -level, -havana_gene, -transcript_name, -transcript_support_level, -tag,
                -havana_transcript, -exon_number, -exon_id, -ont, -protein_id, -ccdsid, -score, -phase,
                -source, -transcript_type, -type)

gene_table <- gtf_df %>%
  filter(seqnames == "chrM")

my_lengths_ont_mt <- my_lengths_ont_df %>%
  filter(transcript_id %in% gene_table$transcript_id) %>%
  left_join(gene_table) %>%
  dplyr::slice(rep(1:n(), read_no)) %>% 
  dplyr::select(-read_no) %>%
  mutate(gene_name = gsub("MT-","",gene_name))

my_lengths_ont_mt$Tissue <- samps$tissue[match(my_lengths_ont_mt$sample_id, samps$name)]

#Combine chimeric transcripts manually
my_lengths_ont_mt$width = ifelse(my_lengths_ont_mt$gene_name %in% c("ATP6","ATP8"), 886, 
                        ifelse(my_lengths_ont_mt$gene_name %in% c("ND4L","ND4"), 1673, my_lengths_ont_mt$width))
my_lengths_ont_mt$gene_name = ifelse(my_lengths_ont_mt$gene_name %in% c("ATP6","ATP8"), "ATP8/ATP6", 
                        ifelse(my_lengths_ont_mt$gene_name %in% c("ND4L","ND4"), "ND4L/ND4", my_lengths_ont_mt$gene_name))

median_lengths <-  my_lengths_ont_mt %>%
  dplyr::group_by(sample_id, gene_name, width, Tissue) %>%
  summarise(median_read_length = median(read_length), sd_read_length = sd(read_length), total_reads = n()) %>%
  ungroup() %>%
  mutate(perc_med = median_read_length/width) %>%
  mutate(perc_sd = sd_read_length/width) %>%
  mutate(sample_id = factor(sample_id))

median_lengths$sample_id2 <- factor(x = as.vector(median_lengths$sample_id), levels = unique(median_lengths[order(median_lengths$Tissue),]$sample_id))
            
#Get names from a unique sample                        
idl = median_lengths %>% filter(sample_id == "QV44_control")
ggplot(median_lengths, aes(x=as.factor(width), y=perc_med, group=sample_id2)) +
  geom_jitter(aes(color=Tissue), size=1, alpha=0.8, position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin=perc_med-perc_sd, ymax=perc_med+perc_sd, color=Tissue), 
                width=0.1, alpha=0.8,
                position=position_dodge(width=0.5)) +
  theme_classic(base_size=14) +
  xlab("Expected length") +
  ylab("Median observed length/expected length") +
  scale_color_manual(values = c("#ffc0cb","#8d5b96","#7776b1","#9773ba","#b873ba","#c893c9",
                                "#ff69b4","#d4a910","#c4625d","#bc3c28","#B09638","#815375",
                                "#0072b5", "#1f854e","#e18726")) +
  geom_text(data = idl, aes(label=gene_name),hjust=0, vjust=0) +
  geom_hline(aes(yintercept = 1), linetype = "dashed")
ggsave("all_samples_mt_lengths.pdf", width = 8, height = 4)

median_lengths_temp <- median_lengths[median_lengths$Tissue=="Cells - Cultured fibroblasts",]
median_lengths_temp$ID <- ifelse(grepl("control",median_lengths_temp$sample_id), "Center2 - PCR",
                                 ifelse(grepl("direct",median_lengths_temp$sample_id), "Center1 - direct","Center1 - PCR"))
median_lengths_temp$sample_id2 <- factor(x = as.vector(median_lengths_temp$sample_id), levels = unique(median_lengths_temp[order(median_lengths_temp$ID),]$sample_id))
ggplot(median_lengths_temp, aes(x=as.factor(width), y=perc_med, group=sample_id2)) +
  geom_jitter(aes(color=ID), size=1, alpha=0.8, position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin=perc_med-perc_sd, ymax=perc_med+perc_sd, color=ID), 
                width=0.1, alpha=0.8,
                position=position_dodge(width=0.5)) +
  theme_classic(base_size=14) +
  xlab("Expected length") +
  ylab("Median observed length/expected length") +
  geom_text(data = idl, aes(label=gene_name),hjust=0, vjust=0) +
  geom_hline(aes(yintercept = 1), linetype = "dashed")
ggsave("fibs_samples_mt_lengths.pdf", width = 8, height = 4)

median_lengths_sum <- median_lengths %>%
  dplyr::mutate(gene_name = as.factor(gene_name)) %>%
  dplyr::group_by(gene_name, width) %>%
  dplyr::summarise(median_read_length2 = median(median_read_length),
            sd_read_length2 = sd(median_read_length),
            total_reads2 = sum(total_reads),
            sample_no = length(unique(sample_id))) %>%
  ungroup() %>%
  mutate(perc_med = median_read_length2/width) %>%
  mutate(perc_sd = sd_read_length2/width)

ggplot(median_lengths_sum, aes(x=width, y=perc_med, label=gene_name)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_jitter(size=3,width = 0.25, color="#0e5a72") +
  geom_errorbar(aes(ymin=perc_med-perc_sd, ymax=perc_med+perc_sd), width=.2,
                position=position_dodge(0.05), color="#0e5a72") +
  theme_classic(base_size=14) +
  xlab("Expected length") +
  xlim(200,2000) +
  ylab("Median fraction of full-length") +
  geom_text(aes(label=gene_name),hjust=1, vjust=-1) +
  ylim(0.2,1.8)
ggsave("all_samples_median_mt_lengths.pdf", width = 4, height = 4)

median_lengths_sum <- median_lengths %>%
  dplyr::group_by(sample_id, Tissue) %>%
  dplyr::summarise(median_read_length2 = median(perc_med),
            sd_read_length2 = sd(perc_med)) %>%
  ungroup() 

median_lengths_sum$ID2 <- samps$ID2[match(median_lengths_sum$sample_id, samps$name)]
median_lengths_sum$SMNABTCHT <- gtex_char$SMNABTCHT[match(median_lengths_sum$ID2,gtex_char$ID2)]
median_lengths_sum$SMNABTCHT[is.na(median_lengths_sum$SMNABTCHT)] <- "RNA isolation_Trizol Manual (Cell Pellet)"

ggplot(median_lengths_sum, aes(x=SMNABTCHT, y=median_read_length2)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_violin(draw_quantiles = c(0.5)) +
  geom_jitter(size=3,width = 0.25) +
  theme_classic(base_size=14) +
  xlab("Method") +
  ylab("Median fraction of full-length") +
  coord_flip()
ggsave("all_samples_median_mt_lengths_method.pdf", width = 8, height = 8)
write.table(median_lengths, "analysis/three_prime_bias/median_lengths_sum.txt", col.names = TRUE, row.names = FALSE, sep="\t", quote = FALSE)
