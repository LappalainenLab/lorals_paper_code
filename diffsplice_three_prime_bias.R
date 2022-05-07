library(tidyverse)

median_lengths_sum <- read_tsv("analysis/three_prime_bias/median_lengths_sum.txt")

top10 <- median_lengths_sum %>% top_frac(0.1, wt = median_read_length2)
bottom10 <- median_lengths_sum %>% top_frac(-0.1, wt = median_read_length2)

## Look into ALT TRANSCRIPT STRUCTURE EVENTS IN FLAIR 
event_char_flair <- read.table("analysis/diffsplice/new_noSJ_promoters_all_events.tsv", header=TRUE)
event_char_flair$diffsplice_events <- gsub("AF5","A5",event_char_flair$diffsplice_events)
event_char_flair$diffsplice_events <- gsub("AL3","A3",event_char_flair$diffsplice_events)
event_char_flair$fish <- paste0(event_char_flair$transcript_id,"_",event_char_flair$gene_id)

## GET HIGHLY EXPRESSED
ont_flair_tpm_df <- read.table("analysis/count_tables/trans_flair.tpm.txt", header = TRUE)
rownames(ont_flair_tpm_df) <- ont_flair_tpm_df$transcript_id 
ont_flair_tpm_df <- ont_flair_tpm_df[,-1]

top10_ont_flair_tpm_df <- ont_flair_tpm_df[,colnames(ont_flair_tpm_df) %in% top10$sample_id,]
top10_ont_flair_tpm_df <- top10_ont_flair_tpm_df[ rowSums(top10_ont_flair_tpm_df >= 1) >= 1, ]

bottom10_ont_flair_tpm_df <- ont_flair_tpm_df[,colnames(ont_flair_tpm_df) %in% bottom10$sample_id,]
bottom10_ont_flair_tpm_df <- bottom10_ont_flair_tpm_df[ rowSums(bottom10_ont_flair_tpm_df >= 1) >= 1, ]

top10_event_char_flair <- event_char_flair[event_char_flair$fish %in% rownames(top10_ont_flair_tpm_df),]
top10_event_char_flair %<>%
  dplyr::select(gene_id, diffsplice_events, event_id) %>%
  unique() %>%
  group_by(diffsplice_events) %>%
  summarise(Freq=n()) %>%
  ungroup() %>%
  mutate(Group = "TOP10") %>%
  mutate(Prop = Freq/sum(Freq))

bottom10_event_char_flair <- event_char_flair[event_char_flair$fish %in% rownames(bottom10_ont_flair_tpm_df),]
bottom10_event_char_flair %<>%
  dplyr::select(gene_id, diffsplice_events, event_id) %>%
  unique() %>%
  group_by(diffsplice_events) %>%
  summarise(Freq=n()) %>%
  ungroup() %>%
  mutate(Group = "BOTTOM10") %>%
  mutate(Prop = Freq/sum(Freq))

final_table <- rbind(bottom10_event_char_flair, top10_event_char_flair)
ggplot(final_table, aes(x=diffsplice_events, y=Prop, fill=Group)) +
  geom_bar(stat="identity", position="dodge") +
  theme_classic(base_size=14) +
  scale_fill_manual(values=c("#a8489b","#99a644")) +
  ylab("Proportion of AltTS events") + xlab("")
ggsave(filename = "three_prime_bias_altts.pdf", height = 4, width = 5.5)

mypvalues <- final_table %>%
  group_by(Group) %>%
  mutate(Freq_all = sum(Freq)) %>%
  ungroup() %>%
  dplyr::select(-Group, -Prop) %>%
  nest_by(diffsplice_events) %>%
  dplyr::mutate(value = list(prop.test(data$Freq, data$Freq_all)$p.value)) %>%
  ungroup() %>%
  dplyr::select(diffsplice_events, value) %>%
  unnest(value) %>%
  mutate(padj = p.adjust(value))
