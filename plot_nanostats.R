library(tidyverse)
library(magrittr)

#Read in result of running submit_nanoplot.sh

#Generates extended figure 1a,b,c

# Load samples metadata and subset to files also sequenced in the GTEx
samps <- read.table("data/metadata.txt",
                    header=TRUE, fill = TRUE, sep = "\t")
samps <- samps[!grepl("CVD", samps$sample_id),]
#too few features
samps <- samps[!samps$path %in% c("FAK91589"),]

all_files <- list.files(path = "analysis/nanoplot_fastq", pattern = "_NanoStats.txt", full.names = TRUE)
nanostats_fastq <- lapply(all_files, read_delim, col_names = FALSE, skip = 1, delim = ":", n_max = 8,trim_ws = TRUE)
names(nanostats_fastq) <-  gsub("_fastq_NanoStats.txt","",basename(all_files))
nanostats_fastq <- bind_rows(nanostats_fastq, .id = "Sample")

all_files <- list.files(path = "analysis/nanoplot_bam", pattern = "_NanoStats.txt", full.names = TRUE)
nanostats_bam <- lapply(all_files, read_delim, col_names = FALSE, skip = 1, delim = ":", n_max = 8,trim_ws = TRUE)
names(nanostats_bam) <- gsub("_bam_NanoStats.txt","",basename(all_files))
nanostats_bam <- bind_rows(nanostats_bam, .id = "Sample")

nanostats_combined <- merge(nanostats_fastq, nanostats_bam, by=c("Sample","X1"), suffix = c("_fastq","_bam"))
nanostats_combined <- merge(nanostats_combined, samps, by.x="Sample", by.y="path")

no_reads <- nanostats_combined[nanostats_combined$X1 == "Number of reads",]
ggplot(no_reads, aes(y=X2_bam, x=X2_fastq)) +
  geom_point(size=5, aes(fill=tissue),colour="grey20",pch=21) +
  theme_classic(base_size=14) +
  theme(legend.position = "right") +
  scale_fill_manual(values = c("#ffc0cb","#8d5b96","#7776b1","#9773ba","#b873ba","#c893c9",
                               "#ff69b4","#d4a910","#c4625d","#bc3c28","#B09638","#815375",
                               "#0072b5", "#1f854e","#e18726")) +
  ylab("Number of aligned reads") + xlab("Number of raw reads") +
  scale_y_continuous(labels = fancy_scientific) +
  scale_x_continuous(labels = fancy_scientific) +
  geom_abline(intercept = 0, color="grey70")
ggsave("number_of_reads.pdf", plot = last_plot(),
       width = 8.5, height = 5.1, units = "in", useDingbats=FALSE)

read_length <- nanostats_combined[nanostats_combined$X1 == "Median read length",]
ggplot(read_length, aes(y=X2_bam, x=X2_fastq)) +
  geom_point(size=5, aes(fill=tissue),colour="grey20",pch=21) +
  theme_classic(base_size=14) +
  theme(legend.position = "right") +
  scale_fill_manual(values = c("#ffc0cb","#8d5b96","#7776b1","#9773ba","#b873ba","#c893c9",
                               "#ff69b4","#d4a910","#c4625d","#bc3c28","#B09638","#815375",
                               "#0072b5", "#1f854e","#e18726")) +
  ylab("Median length of aligned reads") + xlab("Median length of raw reads") +
  geom_abline(intercept = 0, color="grey70")
ggsave("median_read_length.pdf", plot = last_plot(),
       width = 8.5, height = 5.1, units = "in", useDingbats=FALSE)

## Look into direct VS PCR
general_nanobam_FAK46618 <- read.table("analysis/nanoplot_bam/FAK46618_bam_NanoPlot-data.tsv.gz", header=TRUE)
general_nanobam_FAK46618$Sample <- "FAK46618"
general_nanobam_FAK46687 <- read.table("analysis/nanoplot_bam/FAK46687_bam_NanoPlot-data.tsv.gz", header=TRUE)
general_nanobam_FAK46687$Sample <- "FAK46687"
general_nanobam_FAK46867 <- read.table("analysis/nanoplot_bam/FAK46867_bam_NanoPlot-data.tsv.gz", header=TRUE)
general_nanobam_FAK46867$Sample <- "FAK46867"
general_nanobam_FAK49083 <- read.table("analysis/nanoplot_bam/FAK49083_bam_NanoPlot-data.tsv.gz", header=TRUE)
general_nanobam_FAK49083$Sample <- "FAK49083"
general_nanobam_FAK49207 <- read.table("analysis/nanoplot_bam/FAK49207_bam_NanoPlot-data.tsv.gz", header=TRUE)
general_nanobam_FAK49207$Sample <- "FAK49207"

df <- rbind(general_nanobam_FAK46618, general_nanobam_FAK46687, general_nanobam_FAK46867, general_nanobam_FAK49083, general_nanobam_FAK49207)

df1 <- merge(df, samps, by.y = "path", by.x="Sample")
ggplot(df1, aes(x=log2(aligned_lengths), color = sample_name, linetype=protocol)) +
  geom_density() +
  scale_color_manual(values = c("#c28834","#45a19a")) +
  xlab("LOG2 Aligned read length") +
  theme_classic(base_size=14) +
  theme(legend.position = "bottom")
ggsave("aligned_length_dist_fibs.pdf", plot = last_plot(),
       width = 6, height = 4, units = "in", useDingbats=FALSE)

summary <- data_summary(df, "aligned_lengths", "Sample")
summary1 <- merge(summary, samps, by.y = "path", by.x="Sample")
ggplot(summary1, aes(x=sample_id, y=aligned_lengths, group=protocol)) +
  geom_bar(position="dodge", stat="identity", fill="#d4aa2a", width = 1.5) +
  facet_wrap(~protocol) +
  geom_errorbar(aes(ymin=aligned_lengths-sd, ymax=aligned_lengths+sd), width=.2,
                position=position_dodge(.9)) +
  ylab("Aligned read length") + xlab("") +
  theme_classic(base_size=14) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave("median_read_length_fibs.pdf", plot = last_plot(),
       width = 6, height = 4, units = "in", useDingbats=FALSE)

tt1 <- t.test(summary1[summary1$protocol=="cDNA-direct",]$aligned_lengths,
       summary1[summary1$protocol=="cDNA-PCR",]$aligned_lengths, paired = FALSE)
print(tt$p.value)

no_reads_test <- no_reads[no_reads$Sample %in% summary1$Sample,]
ggplot(no_reads_test, aes(x=sample_id, y=X2_bam, group=protocol)) +
  geom_bar(position="dodge", stat="identity", fill="#d4aa2a", width = 1.5) +
  facet_wrap(~protocol) +
  ylab("Aligned read number") + xlab("") +
  theme_classic(base_size=14) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave("read_count_fibs.pdf", plot = last_plot(),
       width = 6, height = 4, units = "in", useDingbats=FALSE)

tt2 <- t.test(no_reads_test[no_reads_test$protocol=="cDNA-direct",]$X2_bam,
       no_reads_test[no_reads_test$protocol=="cDNA-PCR",]$X2_bam, paired = FALSE)
print(tt2$p.value)
