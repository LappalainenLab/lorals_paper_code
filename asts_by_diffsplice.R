library(tidyverse)
library(magrittr)
library(qvalue)

### Analyse ASTS events quant
asts_quant_ont_flair_uniq_df <- read_tsv("analysis/ase/flair/all_asts_quant_ont_flair_df.txt")

### ANALYSE ASTS EVENTS by AS CATEGORY
events <- read.table("analysis/diffsplice/new_noSJ_promoters_all_events.tsv", header=TRUE)
events <- events[events$transcript_id %in% asts_quant_ont_flair_uniq_df$Transcript,]
events$diffsplice_events <- gsub("AF5","A5",events$diffsplice_events)
events$diffsplice_events <- gsub("AL3","A3",events$diffsplice_events)

tracking_char_events <- merge(events, asts_quant_ont_flair_uniq_df, by.x=c("transcript_id","gene_id"), by.y=c("Transcript","Gene"), all.y=TRUE)
tracking_char_events %<>%
  dplyr::select(transcript_id, gene_id, diffsplice_events, event_id, Type, names, variantID, sum, cohen, pvalue, FDR) %>%
  unique() %>%
  group_by(gene_id, event_id, diffsplice_events, names) %>%
  filter(length(unique(Type)) > 1) %>%
  ungroup()

tracking_char_events_sum <- tracking_char_events %>%
  filter(FDR <= 0.05) %>%
  select(gene_id, event_id, diffsplice_events, names, variantID, sum, pvalue, FDR) %>%
  unique() %>%
  group_by(gene_id, diffsplice_events, names, variantID, sum, pvalue, FDR) %>%
  summarise(No_event=n()) %>%
  ungroup() %>%
  group_by(gene_id, names) %>%
  mutate(Proportion = No_event/sum(No_event)) %>%
  ungroup() %>%
  mutate(Group = "sig")
tracking_char_events_sum <- na.omit(tracking_char_events_sum)

tracking_char_events_sum_all <- tracking_char_events %>%
  select(gene_id, event_id, diffsplice_events, names, variantID, sum, pvalue, FDR) %>%
  unique() %>%
  group_by(gene_id, diffsplice_events, names, variantID, sum, pvalue, FDR) %>%
  summarise(No_event=n()) %>%
  ungroup() %>%
  group_by(gene_id, names) %>%
  mutate(Proportion = No_event/sum(No_event)) %>%
  ungroup() %>%
  mutate(Group = "all")
tracking_char_events_sum_all <- na.omit(tracking_char_events_sum_all)

tracking_char_events_sum <- rbind(tracking_char_events_sum_all, tracking_char_events_sum)

temp_plot <- tracking_char_events_sum %>%
  filter(Proportion==1) %>%
  group_by(Group, diffsplice_events) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(Group) %>%
  mutate(Proportion = (count/sum(count))*100) %>%
  ungroup()
  
ggplot(temp_plot, aes(x=diffsplice_events, fill=Group, group=Group, y=Proportion)) +
  geom_bar(stat="identity", position='dodge') +
  theme_classic(base_size=14) +
  xlab("") +
  ylab("% of genes (single AltTS events)") +
  theme(legend.position = "bottom") +
  scale_fill_manual(values=c("#a5a3a4","#d76127"))

temp_plot_format <- temp_plot %>%
  filter(Group == "all") %>%
  select(diffsplice_events, Proportion)
  
mypvalues <- temp_plot %>%
  filter(Group == "sig") %>%
  select(-Proportion) %>%
  left_join(temp_plot_format) %>%
  mutate(sum = sum(count)) %>%
  dplyr::nest_by(diffsplice_events) %>%
  dplyr::mutate(value = list(binom.test(data$count, data$sum, data$Proportion/100)$p.value)) %>%
  ungroup() %>%
  dplyr::select(diffsplice_events, value) %>%
  unnest(value)

### MAKE HEATMAP
heatmap_prep <- tracking_char_events_sum %>%
  group_by(gene_id, names, Group) %>%
  mutate(combined_ats = paste0(sort(unique(diffsplice_events)), collapse = "_")) %>%
  ungroup() %>%
  group_by(combined_ats, Group) %>%
  summarise(Freq = n()) %>%
  mutate(names = ifelse(combined_ats == "SE", "SE_SE",
                        ifelse(combined_ats == "RI", "RI_RI",
                               ifelse(combined_ats == "MX", "MX_MX",
                                      ifelse(combined_ats == "AL", "AL_AL",
                                             ifelse(combined_ats == "AF", "AF_AF",
                                                    ifelse(combined_ats == "ALT5UTR", "ALT5UTR_ALT5UTR",
                                                           ifelse(combined_ats == "ALT3UTR", "ALT3UTR_ALT3UTR",
                                                                  ifelse(combined_ats == "A5", "A5_A5",
                                                                         ifelse(combined_ats == "A3", "A3_A3",combined_ats)))))))))) %>%
  group_by(names, Group) %>%
  summarise(Freq = sum(Freq)) %>%
  separate(names, sep="_",into =c("ats1", "ats2"), extra = "merge")

heatmap_prep <- tracking_char_events_sum %>%
  #filter(Group == "all") %>%
  group_by(gene_id, names, Group) %>%
  mutate(combined_ats = paste0(sort(unique(diffsplice_events)), collapse = "_")) %>%
  ungroup() %>%
  group_by(combined_ats, Group) %>%
  summarise(Freq = n()) %>%
  mutate(names = ifelse(combined_ats == "SE", "SE_SE",
                        ifelse(combined_ats == "RI", "RI_RI",
                               ifelse(combined_ats == "MX", "MX_MX",
                                      ifelse(combined_ats == "AL", "AL_AL",
                                             ifelse(combined_ats == "AF", "AF_AF",
                                                    ifelse(combined_ats == "ALT5UTR", "ALT5UTR_ALT5UTR",
                                                           ifelse(combined_ats == "ALT3UTR", "ALT3UTR_ALT3UTR",
                                                                  ifelse(combined_ats == "A5", "A5_A5",
                                                                         ifelse(combined_ats == "A3", "A3_A3",combined_ats)))))))))) %>%
  group_by(names, Group) %>%
  summarise(Freq = sum(Freq)) %>%
  separate(names, sep="_",into =c("ats1", "ats2"), extra = "merge")


heatmap_prep1 <- heatmap_prep[heatmap_prep$Group=="all",]
heatmap_prep2 <- heatmap_prep[heatmap_prep$Group=="sig",]
heatmap_prep <- merge(heatmap_prep1, heatmap_prep2, by=c("ats1","ats2"), suffixes = c("_all","_sig"), all=TRUE)
heatmap_prep$Freq_sig[is.na(heatmap_prep$Freq_sig)] <- 0
heatmap_prep$Proportion <- heatmap_prep$Freq_sig/heatmap_prep$Freq_all

#heatmap_prep_temp <- heatmap_prep[heatmap_prep$ats1 != heatmap_prep$ats2,]
heatmap_prep_temp <- heatmap_prep
heatmap_prep_temp %<>%
  mutate(Sum_all = sum(Freq_all)) %>%
  mutate(Sum_sig = sum(Freq_sig)) %>%
  mutate(ratio_sig = Freq_sig/Sum_sig) %>%
  mutate(ratio_all = Freq_all/Sum_all)
bt2 <- function(a, b, p) {binom.test(a, b, p, alternative=c("two.sided"), conf.level = 0.95)$p.value}
heatmap_prep_temp$pvalue <- mapply(bt2, heatmap_prep_temp$Freq_sig, heatmap_prep_temp$Sum_sig, heatmap_prep_temp$Freq_all/heatmap_prep_temp$Sum_all)
heatmap_prep_temp$qvalue <- p.adjust(heatmap_prep_temp$pvalue)

heatmap_prep <- heatmap_prep[!grepl("_",heatmap_prep$ats2),]
ggplot(heatmap_prep, aes(ats1, ats2, fill=Proportion)) +
  geom_tile(color="grey") +
  theme_classic(base_size=14) +
  scale_fill_gradient2(low = "#23847c", high = "#a8469a", na.value = NA) +
  geom_text(aes(label = round(Freq_sig, 1)))


### ANALYSE ASTS EVENTS by AS CATEGORY TO LOOK AT LENGTH
gff_class_sub <- gff_class %>%
  filter(transcript_id %in% tracking_char_events$transcript_id) %>%
  mutate(len = as.integer(len)) %>%
  group_by(gene_id) %>%
  summarise(mean_len = mean(len)) %>%
  ungroup()

tracking_char_events_sum2 <- merge(tracking_char_events_sum, gff_class_sub, by.x="gene_id", by.y="gene_id")
ggplot(tracking_char_events_sum2, aes(x=mean_len, color=Group)) +
  geom_density() +
  facet_grid(.~diffsplice_events)+
  theme_classic(base_size=14) +
  xlab("Mean gene length") +
  theme(legend.position = "bottom") +
  scale_color_manual(values=c("#a5a3a4","#d76127"))

tracking_char_events_sum2 %>% 
  select(Group, mean_len, diffsplice_events) %>% 
  group_by(diffsplice_events, Group) %>%
  summarise(mean_len = list(mean_len)) %>% 
  spread(Group, mean_len) %>%
  group_by(diffsplice_events) %>% 
  mutate(p_value = wilcox.test(unlist(all), unlist(sig))$p.value,
         t_value = wilcox.test(unlist(all), unlist(sig))$statistic)
