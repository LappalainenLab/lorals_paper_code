library(dplyr)
library(seqminer)
library(readr)
library(tidyr)
library(purrr)
library(ggplot2)
library(magrittr)

vep_severity <- read.table("data/vep_severity.txt")

gencode <- read_tsv("analysis/vep/gtex.most_severe_gencode.txt", col_names = TRUE)
flair <- read_tsv("analysis/vep/gtex.most_severe_flair.txt", col_names = TRUE)

vep_anno <- merge(gencode[,c(1,2,7)], flair[,c(1,2,7)], suffixes = c("_GENCODE","_FLAIR"), by=c("Uploaded_variation","Location"))

tab_long <- tab %>%
  pivot_longer(cols = c("Consequence_GENCODE","Consequence_FLAIR")) %>%
  group_by(name, value) %>%
  summarise(Freq = n()) %>%
  filter(!value %in% c("intron_variant","intergenic_variant")) %>%
  mutate(value =  gsub("_variant","",value))

tab_long$Severity <- vep_severity$V2[match(tab_long$value,gsub("_variant","",vep_severity$V1))]

ggplot(tab_long, aes(x=reorder(value, Severity), y=Freq, fill=name)) +
  geom_bar(stat="identity", position="dodge") +
  theme_classic(base_size=14) +
  theme(legend.position = "bottom") +
  coord_flip() +
  ylab("Number of variants") +
  xlab("VEP Consequence") +
  scale_fill_manual(values = c("#23847c", "#8c5322")) 
ggsave("barplot_num_reassign.pdf", width = 6.7, height = 7.5)

#Look at difference between reassigned and non
reassign <- tab %>%
  filter(Consequence_FLAIR != Consequence_GENCODE) %>%
  separate(Uploaded_variation, into = c("chr","pos","ref","alt","build"), sep = "_") %>%
  mutate(pos = as.numeric(pos)) %>%
  mutate(ranges = paste0(gsub("chr","",chr),":",pos-1,"-",pos))

my_ranges <- paste0( reassign$ranges, collapse="," )

snp_cadd <- tabix.read.table("/gpfs/commons/groups/lappalainen_lab/data/cadd/v1.5/whole_genome_SNVs.tsv.gz", my_ranges)
snp_cadd$Chrom <- gsub("^","chr",snp_cadd$Chrom)

reassign <- merge(reassign, snp_cadd, by.x=c("chr","pos","ref","alt"), by.y=c("Chrom","Pos","Ref","Alt"))
reassign <- reassign %>%
  distinct()

reassign <- droplevels(reassign)
reassign$Severity <- vep_severity$V2[match(reassign$Consequence_FLAIR, vep_severity$V1)]

reassign_numbers <- reassign %>%
  arrange(Consequence_GENCODE) %>%
  group_by(Consequence_GENCODE) %>%
  summarise(Count = n())

set.seed(90210)
same_assign <- tab %>%
  filter(Consequence_FLAIR == Consequence_GENCODE) %>%
  filter(Consequence_GENCODE %in% reassign_numbers$Consequence_GENCODE) %>%
  arrange(Consequence_GENCODE) %>%
  group_by(Consequence_GENCODE) %>% 
  nest() %>%            
  ungroup() %>% 
  mutate(n = reassign_numbers$Count) %>% 
  mutate(samp = map2(data, n, sample_n)) %>% 
  select(-data) %>%
  unnest(samp) %>%
  separate(Uploaded_variation, into = c("chr","pos","ref","alt","build"), sep = "_") %>%
  mutate(pos = as.numeric(pos)) %>%
  mutate(ranges = paste0(gsub("chr","",chr),":",as.integer(pos-1),"-",as.integer(pos)))

my_ranges <- paste0( same_assign$ranges, collapse=",")

snp_cadd_n <- tabix.read.table("/gpfs/commons/groups/lappalainen_lab/data/cadd/v1.5/whole_genome_SNVs.tsv.gz", my_ranges)
snp_cadd_n$Chrom <- gsub("^","chr",snp_cadd_n$Chrom)

same_assign1 <- merge(same_assign, snp_cadd_n, by.x=c("chr","pos","ref","alt"), by.y=c("Chrom","Pos","Ref","Alt"))
same_assign1 <- same_assign1 %>%
  distinct() %>%
  mutate(Consequence_FLAIR = "same")

same_assign1 <- droplevels(same_assign1)
same_assign1$Severity <- vep_severity$V2[match(same_assign1$Consequence_GENCODE, vep_severity$V1)]

ready_plot <- rbind(reassign, same_assign1[,c(1:4,7,5,8:12)])

ready_plot_temp <- ready_plot %>%
  filter(!Consequence_FLAIR %in% c("splice_donor_variant","splice_region_variant","upstream_gene_variant",
                                   "intron_variant","downstream_gene_variant")) %>%
  filter(!Consequence_GENCODE %in% c("coding_sequence_variant","intergenic_variant","upstream_gene_variant",
                                     "downstream_gene_variant")) %>%
  mutate(Consequence_FLAIR =  gsub("_variant","",Consequence_FLAIR)) %>%
  mutate(Consequence_GENCODE =  gsub("_variant","",Consequence_GENCODE))

mypal <- c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02",
           "grey","#027cd9","#9e1b42","#aeb370","#666666")
ggplot(ready_plot_temp, aes(x=reorder(Consequence_FLAIR, Severity), y=PHRED, fill=reorder(Consequence_FLAIR, Severity))) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic(base_size=14) +
  facet_wrap(.~Consequence_GENCODE, scale="free") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "bottom") +
  xlab("Reassigned VEP Consequence") +
  ylab("CADD Score") +
  scale_fill_manual(values = mypal, name = "VEP Consequence")
ggsave("reassigned_cadd_boxplot.pdf", width = 8, height = 7.5)

library(rstatix)
mypvalues <- ready_plot %>%
  group_by(Consequence_GENCODE, Consequence_FLAIR) %>%
  filter(n() > 1) %>%
  ungroup() %>%
  group_by(Consequence_GENCODE) %>%
  t_test(data =., PHRED ~ Consequence_FLAIR) %>%
  filter(group1 == "same" | group2 == "same") %>%
  mutate(Consequence_FLAIR = ifelse(group1 == "same", group2, group1))

mypvalues$Consequence_FLAIR <- gsub("_variant","",mypvalues$Consequence_FLAIR)
mypvalues$Consequence_GENCODE <- gsub("_variant","",mypvalues$Consequence_GENCODE)

### MAKE HEATMAP
reassign_summ <- reassign %>%
  group_by(Consequence_FLAIR, Consequence_GENCODE, Severity) %>%
  summarise(mean_CADD = mean(PHRED), Count = n()) %>%
  ungroup()
same_assign1_summ <- same_assign1 %>%
  group_by(Consequence_GENCODE, Severity) %>%
  summarise(null_CADD = mean(PHRED), Null_count = n()) %>%
  ungroup()

ready_plot_summ <- merge(reassign_summ, same_assign1_summ, by="Consequence_GENCODE", suffixes = c("_reass","_same"))
ready_plot_summ <- ready_plot_summ %>%
  mutate(Delta = mean_CADD - null_CADD)

ready_plot_summ$Consequence_FLAIR <- gsub("_variant","",ready_plot_summ$Consequence_FLAIR)
ready_plot_summ$Consequence_GENCODE <- gsub("_variant","",ready_plot_summ$Consequence_GENCODE)

ready_plot_summ <- left_join(ready_plot_summ, mypvalues)
ready_plot_summ <- ready_plot_summ %>% dplyr::mutate(p = replace_na(p, 1))
ready_plot_summ$sig <- ifelse(ready_plot_summ$p < 0.01, "yes", "no")

ggplot(ready_plot_summ, aes(x=reorder(Consequence_GENCODE,Severity_same), y=reorder(Consequence_FLAIR, Severity_reass), size=log10(Count))) +
  geom_point(aes(fill=Delta, shape=sig)) +
  theme_classic(base_size=14) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Reassigned VEP Consequence") +
  xlab("Previous VEP Consequence") +
  scale_fill_gradient2(low = "#23847c", high = "#8c5322", na.value = NA, midpoint = 0) + 
  scale_shape_manual(values = c(21, 23)) 
ggsave("heatmap_cadd_delta.pdf", width = 6.7, height = 7.5)
