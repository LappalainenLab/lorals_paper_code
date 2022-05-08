library(tidyverse)
library(gridExtra)
library(qvalue)

#Generate extended Figure 9A-D, 11A,B,C,F

#Meta/background files
# Downloaded this file from the GTEx portal
tissue_anot_file = read_tsv('data/tissue_abbreviations.txt')

samps <- read.table("data/metadata.txt",
                    header=TRUE, fill = TRUE, sep = "\t")

samps <- samps[!grepl("experimental", samps$name),]
samps <- samps[!grepl("direct", samps$name),]
samps <- samps[!grepl("CVD", samps$name),]
samps <- samps[!grepl("K562", samps$name),]
samps$sample_name = gsub("_","-",samps$sample_name)
samps <- samps[!samps$sample_name %in% c("GTEX-1HBPH","GTEX-1IDJU","GTEX-11H98","GTEX-1H3NZ","GTEX-1H11D"),]
samps <- samps[samps$name != "GTEX_T5JC_brainfrontalcortex",]
samps$tissue_abbrv <- tissue_anot_file$tissue_abbrv[match(samps$tissue, tissue_anot_file$tissue_name)]
samps <- samps[!grepl("rep", samps$name),]

### LOAD ILLUMINA (downloaded from TERRA)
file_names <- unique(paste0("data/illumina_ase/",
                            samps$sample_name,".v8.ase_table.",samps$tissue_abbrv,"_only.txt"))
dataFiles_illumina <- lapply(file_names, read_tsv,
                             col_names=c("CHR","POS","VARIANT_ID","REF_ALLELE","ALT_ALLELE","SAMPLE_ID",
                                         "SUBJECT_ID","TISSUE_ID","REF_COUNT","ALT_COUNT","TOTAL_COUNT",
                                         "REF_RATIO","OTHER_ALLELE_COUNT","NULL_RATIO","BINOM_P",
                                         "BINOM_P_ADJUSTED","MAMBA_POST_SINGLETIS","MAMBA_POST_MULTITIS",
                                         "GENOTYPE","VARIANT_ANNOTATION","GENE_ID","LOW_MAPABILITY",
                                         "MAPPING_BIAS_SIM","GENOTYPE_WARNING"))
names(dataFiles_illumina) <- samps$name

dataFiles_illumina_filtered <- lapply(dataFiles_illumina, completeFun, "GENE_ID")
dataFiles_illumina_filtered <- lapply(dataFiles_illumina_filtered, function(x) filter(x, TOTAL_COUNT >= 20))
dataFiles_illumina_filtered <- lapply(dataFiles_illumina_filtered, function(x) filter(x, CHR != "chrX"))
dataFiles_illumina_filtered <- lapply(dataFiles_illumina_filtered, function(x) filter(x, CHR != "chrY"))
dataFiles_illumina_filtered <- lapply(dataFiles_illumina_filtered, function(x) filter(x, CHR != "chrM"))

### ONT LOAD
file_names <- unique(paste0("analysis/ase/gencode/snp_counts.a",
                            samps$path,".txt"))
dataFiles_ont <- lapply(file_names, read_tsv)
names(dataFiles_ont) <- samps$name

dataFiles_ont[["GTEX_WY7C_heartatrialappendage"]] <- read_tsv("ont/gencode/snp_counts.aGTEX_WY7C_heartatrialappendage.txt")
dataFiles_ont[["GTEX_13QJ3_muscle"]] <- read_tsv("ont/gencode/snp_counts.aGTEX_13QJ3_muscle.txt")
dataFiles_ont[["GTEX_14BMU_lung"]] <- read_tsv("ont/gencode/snp_counts.aGTEX_14BMU_lung.txt")
dataFiles_ont[["GTEX_Q2AG_braincerebellarhemisphere"]] <- read_tsv("ont/gencode/snp_counts.aGTEX_Q2AG_braincerebellarhemisphere.txt")
dataFiles_ont[["GTEX_ZT9X_muscle"]] <- read_tsv("ont/gencode/snp_counts.aGTEX_ZT9X_muscle.txt")
dataFiles_ont[["GTEX_Y5LM_liver"]] <- read_tsv("ont/gencode/snp_counts.aGTEX_Y5LM_liver.txt")
dataFiles_ont[["GTEX_1GN1W_heartatrialappendage"]] <- read_tsv("ont/gencode/snp_counts.aGTEX_1GN1W_heartatrialappendage.txt")
dataFiles_ont[["GTEX_SM3QNG2_fibs"]] <- read_tsv("ont/gencode/snp_counts.aGTEX_RWS6_fibs.txt")
dataFiles_ont[["QV44_control"]] <- read_tsv("ont/gencode/snp_counts.aQV44_control.txt")

dataFiles_ont_filtered <- lapply(dataFiles_ont, function(x) filter(x, GENE_ID != "."))
dataFiles_ont_filtered <- lapply(dataFiles_ont_filtered, function(x) filter(x, totalCount >= 20))
dataFiles_ont_filtered <- lapply(dataFiles_ont_filtered, function(x) filter(x, contig != "chrX"))
dataFiles_ont_filtered <- lapply(dataFiles_ont_filtered, function(x) filter(x, contig != "chrY"))
dataFiles_ont_filtered <- lapply(dataFiles_ont_filtered, function(x) filter(x, contig != "chrM"))

# Split each line so that there is a single gene
s = list()
for (i in names(dataFiles_ont_filtered)){
  s[[i]] <- strsplit(dataFiles_ont_filtered[[i]]$GENE_ID, split = ",")
  dataFiles_ont_filtered[[i]] <- data.frame(variantID = rep(dataFiles_ont_filtered[[i]]$variantID, sapply(s[[i]], length)),
                                            refCount = rep(dataFiles_ont_filtered[[i]]$refCount, sapply(s[[i]], length)),
                                            altCount = rep(dataFiles_ont_filtered[[i]]$altCount, sapply(s[[i]], length)),
                                            totalCount = rep(dataFiles_ont_filtered[[i]]$totalCount, sapply(s[[i]], length)),
                                            refIndelCount = rep(dataFiles_ont_filtered[[i]]$refIndelCount, sapply(s[[i]], length)),
                                            altIndelCount = rep(dataFiles_ont_filtered[[i]]$altIndelCount, sapply(s[[i]], length)),
                                            otherBases = rep(dataFiles_ont_filtered[[i]]$otherBases, sapply(s[[i]], length)),
                                            GENOTYPE = rep(dataFiles_ont_filtered[[i]]$GENOTYPE, sapply(s[[i]], length)),
                                            NULL_RATIO = rep(dataFiles_ont_filtered[[i]]$NULL_RATIO, sapply(s[[i]], length)),
                                            MULTI_MAPPING = rep(dataFiles_ont_filtered[[i]]$MULTI_MAPPING, sapply(s[[i]], length)),
                                            BLACKLIST = rep(dataFiles_ont_filtered[[i]]$BLACKLIST, sapply(s[[i]], length)),
                                            GENOTYPE_WARNING = rep(dataFiles_ont_filtered[[i]]$GENOTYPE_WARNING, sapply(s[[i]], length)),
                                            HIGH_INDEL_WARNING = rep(dataFiles_ont_filtered[[i]]$HIGH_INDEL_WARNING, sapply(s[[i]], length)),
                                            OTHER_ALLELE_WARNING = rep(dataFiles_ont_filtered[[i]]$OTHER_ALLELE_WARNING, sapply(s[[i]], length)),
                                            BINOM_P = rep(dataFiles_ont_filtered[[i]]$BINOM_P, sapply(s[[i]], length)),
                                            BINOM_P_ADJUSTED = rep(dataFiles_ont_filtered[[i]]$BINOM_P_ADJUSTED, sapply(s[[i]], length)),
                                            GENE_ID = unlist(s[[i]]))
  dataFiles_ont_filtered[[i]]$GENE_ID <- gsub("\\..*","",dataFiles_ont_filtered[[i]]$GENE_ID)
}

dataFiles_ont_filtered_df <- bind_rows(dataFiles_ont_filtered, .id = "names")
pa <- ggplot(dataFiles_ont_filtered_df, aes(x=(refIndelCount+altIndelCount)/(totalCount+refIndelCount+altIndelCount+otherBases))) +
  geom_histogram(fill="#0f475e", color="grey30", bins=40) +
  theme_classic(base_size=14) +
  geom_vline(xintercept = 0.2, linetype="dashed") +
  xlab("Ratio of reads containing INDELs")

pb <- ggplot(dataFiles_ont_filtered_df, aes(x=otherBases/(totalCount+refIndelCount+altIndelCount+otherBases))) +
  geom_histogram(fill="#0f475e", color="grey30", bins=40) +
  theme_classic(base_size=14) +
  geom_vline(xintercept = 0.2, linetype="dashed") +
  xlab("Ratio of reads containing other allele to REF and ALT")

multiplot(pa,pb, cols=2)

dataFiles_ont_filtered_df_temp <- dataFiles_ont_filtered_df[order(dataFiles_ont_filtered_df$GENE_ID, (dataFiles_ont_filtered_df$totalCount)),]
dataFiles_ont_filtered_df_temp <- dataFiles_ont_filtered_df_temp[!duplicated(dataFiles_ont_filtered_df_temp[c("GENE_ID","names")]), ]  
sum_dataFiles_ont <- dataFiles_ont_filtered_df_temp %>%
  select(names, GENE_ID,MULTI_MAPPING,BLACKLIST,GENOTYPE_WARNING,HIGH_INDEL_WARNING,OTHER_ALLELE_WARNING) %>%
  pivot_longer(cols = c("MULTI_MAPPING","BLACKLIST","GENOTYPE_WARNING","HIGH_INDEL_WARNING","OTHER_ALLELE_WARNING")) %>%
  group_by(names, name, value) %>%
  summarise(freq = n()) %>%
  ungroup() %>%
  group_by(names, name) %>%
  mutate(ratio = freq/sum(freq)) %>%
  filter(value == 1)
ggplot(sum_dataFiles_ont, aes(x=name, y=ratio)) +
  geom_boxplot() +
  xlab("") + ylab("Proportion of genes which are flagged") +
  coord_flip() +
  theme_classic()
ggsave("proportion_flagged_sites.pdf",width=6.3, height = 3.3)

dataFiles_illumina_filtered_df <- bind_rows(dataFiles_illumina_filtered, .id = "names")
dataFiles_illumina_filtered_df$VARIANT_ID <- gsub("_b38","",dataFiles_illumina_filtered_df$VARIANT_ID)

dataFiles_df <- merge(dataFiles_ont_filtered_df, dataFiles_illumina_filtered_df, by.x=c("names", "GENE_ID","variantID"), by.y=c("names", "GENE_ID","VARIANT_ID"), all=TRUE)

#Keep variant with the highest count in both ILLUMINA AND ONT
dataFiles_df <- dataFiles_df[order(dataFiles_df$GENE_ID, (dataFiles_df$totalCount+dataFiles_df$TOTAL_COUNT)),]
dataFiles_df <- dataFiles_df[!duplicated(dataFiles_df[c("GENE_ID","names")]), ]  

### RUN ONE FLAG AT A TIME
## ILLUMINA FLAGS
temp <- dataFiles_df[dataFiles_df$MAPPING_BIAS_SIM == 1,]
p1 <- ggplot(temp) +
  geom_histogram(aes(x=(refCount/totalCount)), fill="#0f475e", color="grey30", alpha=0.5, bins=40) +
  geom_histogram(aes(x=REF_RATIO), fill="#fda633", color="grey30", alpha=0.5, bins=40) +
  theme_classic(base_size=14) +
  xlab("Ref ratios") +
  ggtitle("Mapping bias")

p1a <- ggplot(temp) +
  geom_histogram(aes(x=log2(totalCount)), fill="#0f475e", color="grey30", alpha=0.5, bins=40) +
  geom_histogram(aes(x=log2(TOTAL_COUNT)), fill="#fda633", color="grey30", alpha=0.5, bins=40) +
  theme_classic(base_size=14) +
  xlab("LOG2 Total Counts") +
  ggtitle("Mapping bias")

temp <- dataFiles_df[dataFiles_df$LOW_MAPABILITY == 1,]
p2 <- ggplot(temp) +
  geom_histogram(aes(x=(refCount/totalCount)), fill="#0f475e", color="grey30", alpha=0.5, bins=40) +
  geom_histogram(aes(x=REF_RATIO), fill="#fda633", color="grey30", alpha=0.5, bins=40) +
  theme_classic(base_size=14) +
  xlab("Ref ratios") +
  ggtitle("Low mappability")

p2a <- ggplot(temp) +
  geom_histogram(aes(x=log2(totalCount)), fill="#0f475e", color="grey30", alpha=0.5, bins=40) +
  geom_histogram(aes(x=log2(TOTAL_COUNT)), fill="#fda633", color="grey30", alpha=0.5, bins=40) +
  theme_classic(base_size=14) +
  xlab("LOG2 Total Counts") +
  ggtitle("Low mappability")

temp <- dataFiles_df[dataFiles_df$GENOTYPE_WARNING.y == 1,]
p3 <- ggplot(temp) +
  geom_histogram(aes(x=(refCount/totalCount)), fill="#0f475e", color="grey30", alpha=0.5, bins=40) +
  geom_histogram(aes(x=REF_RATIO), fill="#fda633", color="grey30", alpha=0.5, bins=40) +
  theme_classic(base_size=14) +
  xlab("Ref ratios") +
  ggtitle("Genotype warning")

p3a <- ggplot(temp) +
  geom_histogram(aes(x=log2(totalCount)), fill="#0f475e", color="grey30", alpha=0.5, bins=40) +
  geom_histogram(aes(x=log2(TOTAL_COUNT)), fill="#fda633", color="grey30", alpha=0.5, bins=40) +
  theme_classic(base_size=14) +
  xlab("LOG2 Total Counts") +
  ggtitle("Genotype warning")

temp <- dataFiles_df[dataFiles_df$OTHER_ALLELE_WARNING == 1,]
p4 <- ggplot(temp) +
  geom_histogram(aes(x=(refCount/totalCount)), fill="#0f475e", color="grey30", alpha=0.5, bins=40) +
  geom_histogram(aes(x=REF_RATIO), fill="#fda633", color="grey30", alpha=0.5, bins=40) +
  theme_classic(base_size=14) +
  xlab("Ref ratios") +
  ggtitle("Other allele warning")

p4a <- ggplot(temp) +
  geom_histogram(aes(x=log2(totalCount)), fill="#0f475e", color="grey30", alpha=0.5, bins=40) +
  geom_histogram(aes(x=log2(TOTAL_COUNT)), fill="#fda633", color="grey30", alpha=0.5, bins=40) +
  theme_classic(base_size=14) +
  xlab("LOG2 Total Counts") +
  ggtitle("Other allele warning")

temp <- dataFiles_df[dataFiles_df$BLACKLIST == 1,]
p5 <- ggplot(temp) +
  geom_histogram(aes(x=(refCount/totalCount)), fill="#0f475e", color="grey30", alpha=0.5, bins=40) +
  geom_histogram(aes(x=REF_RATIO), fill="#fda633", color="grey30", alpha=0.5, bins=40) +
  theme_classic(base_size=14) +
  xlab("Ref ratios") +
  ggtitle("Blacklist")

p5a <- ggplot(temp) +
  geom_histogram(aes(x=log2(totalCount)), fill="#0f475e", color="grey30", alpha=0.5, bins=40) +
  geom_histogram(aes(x=log2(TOTAL_COUNT)), fill="#fda633", color="grey30", alpha=0.5, bins=40) +
  theme_classic(base_size=14) +
  xlab("LOG2 Total Counts") +
  ggtitle("Blacklist")

temp <- dataFiles_df[dataFiles_df$MULTI_MAPPING == 1,]
p6 <- ggplot(temp) +
  geom_histogram(aes(x=(refCount/totalCount)), fill="#0f475e", color="grey30", alpha=0.5, bins=40) +
  geom_histogram(aes(x=REF_RATIO), fill="#fda633", color="grey30", alpha=0.5, bins=40) +
  theme_classic(base_size=14) +
  xlab("Ref ratios") +
  ggtitle("Multi-mapping")

p6a <- ggplot(temp) +
  geom_histogram(aes(x=log2(totalCount)), fill="#0f475e", color="grey30", alpha=0.5, bins=40) +
  geom_histogram(aes(x=log2(TOTAL_COUNT)), fill="#fda633", color="grey30", alpha=0.5, bins=40) +
  theme_classic(base_size=14) +
  xlab("LOG2 Total Counts") +
  ggtitle("Multi-mapping")

temp <- dataFiles_df[dataFiles_df$HIGH_INDEL_WARNING == 1,]
p7 <- ggplot(temp) +
  geom_histogram(aes(x=(refCount/totalCount)), fill="#0f475e", color="grey30", alpha=0.5, bins=40) +
  geom_histogram(aes(x=REF_RATIO), fill="#fda633", color="grey30", alpha=0.5, bins=40) +
  theme_classic(base_size=14) +
  xlab("Ref ratios") +
  ggtitle("High indel warning")

p7a <- ggplot(temp) +
  geom_histogram(aes(x=log2(totalCount)), fill="#0f475e", color="grey30", alpha=0.5, bins=40) +
  geom_histogram(aes(x=log2(TOTAL_COUNT)), fill="#fda633", color="grey30", alpha=0.5, bins=40) +
  theme_classic(base_size=14) +
  xlab("LOG2 Total Counts") +
  ggtitle("High indel warning")



dataFiles_illumina_filtered <- lapply(dataFiles_illumina_filtered, function(x) filter(x, GENOTYPE_WARNING == "0"))
dataFiles_illumina_filtered <- lapply(dataFiles_illumina_filtered, function(x) filter(x, MAPPING_BIAS_SIM == "0"))
dataFiles_illumina_filtered <- lapply(dataFiles_illumina_filtered, function(x) filter(x, LOW_MAPABILITY == "0"))

tt <- dataFiles_ont_filtered_df %>%
  group_by(names) %>%
  summarise(count = n()) %>%
  ungroup()

dataFiles_ont_filtered <- lapply(dataFiles_ont_filtered, function(x) filter(x, OTHER_ALLELE_WARNING == "0"))
dataFiles_ont_filtered <- lapply(dataFiles_ont_filtered, function(x) filter(x, MULTI_MAPPING == "0"))
dataFiles_ont_filtered <- lapply(dataFiles_ont_filtered, function(x) filter(x, BLACKLIST == "0"))
dataFiles_ont_filtered <- lapply(dataFiles_ont_filtered, function(x) filter(x, HIGH_INDEL_WARNING == "0"))
dataFiles_ont_filtered <- lapply(dataFiles_ont_filtered, function(x) filter(x, GENOTYPE_WARNING == "0"))

for (i in names(dataFiles_ont_filtered)){
  dataFiles_ont_filtered[[i]]$BINOM_P_ADJUSTED <- qvalue(dataFiles_ont_filtered[[i]]$BINOM_P)$qvalue
}

dataFiles_ont_filtered_df <- bind_rows(dataFiles_ont_filtered, .id = "names")

tt2 <- dataFiles_ont_filtered_df %>%
  group_by(names) %>%
  summarise(count2 = n()) %>%
  ungroup()

dataFiles_illumina_filtered_df <- bind_rows(dataFiles_illumina_filtered, .id = "names")
dataFiles_illumina_filtered_df$VARIANT_ID <- gsub("_b38","",dataFiles_illumina_filtered_df$VARIANT_ID)

dataFiles_df <- merge(dataFiles_ont_filtered_df, dataFiles_illumina_filtered_df, by.x=c("names", "GENE_ID","variantID"), by.y=c("names", "GENE_ID","VARIANT_ID"), all=TRUE)
dataFiles_df$BINOM_P_ADJUSTED.y[is.na(dataFiles_df$BINOM_P_ADJUSTED.y)] <- "ABS"
dataFiles_df$BINOM_P_ADJUSTED.x[is.na(dataFiles_df$BINOM_P_ADJUSTED.x)] <- "ABS"

dataFiles_df <- dataFiles_df[order(dataFiles_df$GENE_ID, (dataFiles_df$totalCount+dataFiles_df$TOTAL_COUNT)),]
dataFiles_df <- dataFiles_df[!duplicated(dataFiles_df[c("GENE_ID","names")]), ]  

p8 <- ggplot(dataFiles_df) +
  geom_histogram(aes(x=(refCount/totalCount)), fill="#0f475e", color="grey30", alpha=0.5, bins = 40) +
  geom_histogram(aes(x=REF_RATIO), fill="#fda633", color="grey30", alpha=0.5, bins = 40) +
  theme_classic(base_size=14) +
  xlab("Ref ratios") +
  ggtitle("All sites post-filtering")

p8a <- ggplot(dataFiles_df) +
  geom_histogram(aes(x=log2(totalCount)), fill="#0f475e", color="grey30", alpha=0.5, bins=40) +
  geom_histogram(aes(x=log2(TOTAL_COUNT)), fill="#fda633", color="grey30", alpha=0.5, bins=40) +
  theme_classic(base_size=14) +
  xlab("LOG2 Total Counts") +
  ggtitle("All sites post-filtering")

p9 <- ggplot(dataFiles_df) +
  geom_histogram(aes(x=(refCount/totalCount), y=..count../sum(..count..)), fill="#0f475e", color="grey30", alpha=0.5, bins = 40) +
  geom_histogram(aes(x=REF_RATIO, y=..count../sum(..count..)), fill="#fda633", color="grey30", alpha=0.5, bins = 40) +
  theme_classic(base_size=14) +
  xlab("Ref ratios") +
  ggtitle("All sites post-filtering")

p9a <- ggplot(dataFiles_df) +
  geom_histogram(aes(x=log2(totalCount), y=..count../sum(..count..)), fill="#0f475e", color="grey30", alpha=0.5, bins=40) +
  geom_histogram(aes(x=log2(TOTAL_COUNT), y=..count../sum(..count..)), fill="#fda633", color="grey30", alpha=0.5, bins=40) +
  theme_classic(base_size=14) +
  xlab("LOG2 Total Counts") +
  ggtitle("All sites post-filtering")

grid.arrange(p1, p2, p3,
             p5, p6, p4, p7, 
             p8, p9, nrow = 3)
grid.arrange(p1a, p2a, p3a,
             p5a, p6a, p4a, p7a, 
             p8a, p9a, nrow=3)

dataFiles_df$sig_category <- ifelse(dataFiles_df$BINOM_P_ADJUSTED.x <= 0.05 & dataFiles_df$BINOM_P_ADJUSTED.y <= 0.05, "Both",
                                    ifelse(dataFiles_df$BINOM_P_ADJUSTED.x <= 0.05 & dataFiles_df$BINOM_P_ADJUSTED.y == "ABS", "ONT only NA",
                                           ifelse(dataFiles_df$BINOM_P_ADJUSTED.x <= 0.05 & dataFiles_df$BINOM_P_ADJUSTED.y > 0.05, "ONT only sig",
                                                  ifelse(dataFiles_df$BINOM_P_ADJUSTED.y <= 0.05 & dataFiles_df$BINOM_P_ADJUSTED.x == "ABS", "Illumina only NA",
                                                         ifelse(dataFiles_df$BINOM_P_ADJUSTED.y <= 0.05 & dataFiles_df$BINOM_P_ADJUSTED.x > 0.05, "Illumina only sig", "Other"
                                                         )))))
both <- dataFiles_df[dataFiles_df$sig_category == "Both",]
ggplot(both) +
  geom_point(size=3, aes(x=log2(altCount/refCount), y=log2(ALT_COUNT/REF_COUNT)), alpha=0.6) +
  xlab("ONT aFC") + ylab("ILLUMINA aFC") +
  xlim(-7,7) + ylim(-7,7) +
  theme_classic(base_size=14) +
  geom_vline(xintercept = 0, color="grey80") +
  geom_hline(yintercept = 0, color="grey80") +
  geom_text(aes(x=-5,y=-5, label=dim(both[log2(both$altCount/both$refCount) < 0 & log2(both$ALT_COUNT/both$REF_COUNT) < 0,])[1])) +
  geom_text(aes(x=5,y=-5, label=dim(both[log2(both$altCount/both$refCount) > 0 & log2(both$ALT_COUNT/both$REF_COUNT) < 0,])[1])) +
  geom_text(aes(x=5,y=5, label=dim(both[log2(both$altCount/both$refCount) > 0 & log2(both$ALT_COUNT/both$REF_COUNT) > 0,])[1])) +
  geom_text(aes(x=-5,y=5, label=dim(both[log2(both$altCount/both$refCount) < 0 & log2(both$ALT_COUNT/both$REF_COUNT) > 0,])[1]))

both$Direction <- ifelse(log2(both$altCount/both$refCount) * log2(both$ALT_COUNT/both$REF_COUNT) > 0, "Same","Opposite")
both$Tissue <- ifelse(grepl("brain",both$names), "Brain",
                      ifelse(grepl("heart",both$names), "Heart",
                             ifelse(grepl("lung",both$names), "Lung",
                                    ifelse(grepl("muscle",both$names), "Muscle",
                                           ifelse(grepl("liver",both$names), "Liver",
                                                  ifelse(grepl("fib|control",both$names), "Fibroblasts","Adipose"))))))

summ_tissue <- both %>%
  group_by(names, Direction, Tissue) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(names) %>%
  mutate(prop = count/sum(count)) %>%
  ungroup() %>%
  filter(Direction == "Opposite")

ggplot(summ_tissue, aes(x=reorder(names,count), y=count, fill=Tissue)) +
  geom_bar(position = "dodge", stat="identity") +
  theme_classic(base_size=14) +
  scale_fill_manual(values = c("#FF69B4","#9773BA","#D4A910","#BC3C28",
                               "#815375","#0072B5","#1F854E","#E18726")) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  xlab("Samples") + ylab("# of discordant events")

temp_1 <- dataFiles_df[,c("sig_category","TOTAL_COUNT")]
temp_1$Method <- "Illumina"
colnames(temp_1)[2] <- "TOTAL_COUNT"
#temp_1 <- temp_1[temp_1$sig_category != "ONT only NA",]
temp_2 <- dataFiles_df[,c("sig_category","totalCount")]
temp_2$Method <- "ONT"
colnames(temp_2)[2] <- "TOTAL_COUNT"
temp <- rbind(temp_1, temp_2)
temp <- temp[temp$sig_category != "Other",]
temp$Group <- ifelse(temp$sig_category == "Illumina only NA", "Missing",
                     ifelse(temp$sig_category == "ONT only NA", "Missing",
                            ifelse(temp$sig_category == "Both", "Zhared","Significant")))
temp$Meth <- ifelse(temp$sig_category == "Illumina only NA", "Illumina",
                     ifelse(temp$sig_category == "Illumina only sig", "Illumina",
                            ifelse(temp$sig_category == "Both", "Zhared","ONT")))

ggplot(temp, aes(x=log2(TOTAL_COUNT), fill=Method)) +
  geom_histogram() +
  theme_classic(base_size=14) +
  scale_fill_manual(values = c("#fda633","#0f475e")) +
  xlab("LOG2 Total counts") + 
  facet_wrap(Group~Meth, scales = "free") +
  theme(panel.grid.major.x = element_line(colour = "grey80"))
ggsave("depth_significance_ase_lorals.pdf", plot = last_plot(),
       scale = 1, width = 10, height = 4, units = "in", dpi = 300)

#### CALCULATE PI1
table1 <- dataFiles_illumina_filtered_df[order(dataFiles_illumina_filtered_df$GENE_ID, dataFiles_illumina_filtered_df$TOTAL_COUNT),]
table1 <- table1[!duplicated(table1[c("GENE_ID","names")]), ]  

table2 <- dataFiles_ont_filtered_df[order(dataFiles_ont_filtered_df$GENE_ID, dataFiles_ont_filtered_df$totalCount),]
table2 <- table2[!duplicated(table2[c("GENE_ID","names")]), ]  

table1_hits = dplyr::filter(table1, BINOM_P_ADJUSTED < 0.05)
table2_hits = dplyr::semi_join(table2, table1_hits, by = c("GENE_ID","names"))
pi1 = 1 - qvalue::qvalue(table2_hits$BINOM_P, lambda = seq(0, 0.8, 0.001))$pi0
podj <- qvalue::qvalue(table2_hits$BINOM_P, lambda = seq(0, 0.8, 0.001))
hist(podj)
#    0.4150494
table2_hits = dplyr::filter(table2, BINOM_P_ADJUSTED < 0.05)
table1_hits = dplyr::semi_join(table1, table2_hits, by = c("GENE_ID","names"))
pi1 = 1 - qvalue::qvalue(table1_hits$BINOM_P, lambda = seq(0, 0.8, 0.001))$pi0
podj <- qvalue::qvalue(table1_hits$BINOM_P, lambda = seq(0, 0.8, 0.001))
hist(podj)
#  0.2326331
