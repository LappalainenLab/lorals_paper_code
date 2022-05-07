library(tidyverse)
require(RColorBrewer)
library(magrittr)

#Generates Figure 2B

#explicitly calculating ‘out of SE events seen, what percent are novel’, and not doing it in an indirect way through novel transcript identification.
event_char_flair <- read.table("analysis/diffsplice/new_noSJ_promoters_all_events.tsv", header=TRUE)
event_char_flair$diffsplice_events <- gsub("AF5","A5",event_char_flair$diffsplice_events)
event_char_flair$diffsplice_events <- gsub("AL3","A3",event_char_flair$diffsplice_events)

event_char_flair_sum <- event_char_flair %>%
  select(gene_id, diffsplice_events, event_id) %>%
  unique() %>%
  group_by(diffsplice_events) %>%
  summarise(Freq=n()) %>%
  ungroup()

event_char <- read.table("analysis/diffsplice/gencodev26_all_events.tsv", header=TRUE)

tracking <- read_tsv("analysis/flair/new_noSJ_promoters.isoforms_annotation.txt")

event_char_flair$Annotation <- tracking$Annotation[match(event_char_flair$transcript_id, tracking$transcript_id)]
event_char$Annotation <- "Annotated"

temp <- unique(rbind(event_char_flair, event_char))
event_char_all_sum <- temp %>%
  filter(Type == "inclusion") %>%
  group_by(event_id) %>%
  filter(length(unique(Annotation))==1) %>%
  ungroup() %>%
  filter(Annotation == "Novel") %>%
  select(gene_id, diffsplice_events, event_id) %>%
  unique() %>%
  group_by(diffsplice_events) %>%
  summarise(Freq=n()) %>%
  ungroup()

toplot <- merge(event_char_all_sum, event_char_flair_sum, by="diffsplice_events", suffixes = c("_all","_flair"))
toplot$ratio <- toplot$Freq_all/toplot$Freq_flair
ggplot(toplot, aes(x=diffsplice_events, y=ratio, label=Freq_flair)) +
  geom_bar(position="dodge", stat="identity", color="grey30", fill="#c2524f") +
  theme_classic(base_size=14) +
  ylab("Proportion of novel AltTS events") +
  xlab("") +
  theme(axis.title.y=element_blank(),
        panel.grid.major.x = element_line(color = "grey90")) +
  coord_flip() 

ggsave(filename = "proportion_altts_events.psf", width = 5.5, height = 4)




  geom_label(position="dodge", stat="identity")