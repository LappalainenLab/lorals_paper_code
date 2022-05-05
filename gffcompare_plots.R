library(tidyverse)
require(RColorBrewer)
library(magrittr)

#Generates Figure2A and extended Figure 4B
tracking <- read_tsv("analysis/flair/new_noSJ_promoters.isoforms_annotation.txt")

mycolors = c(brewer.pal(name="Set3", n = 5))

tracking$Class_code_sum <- ifelse(tracking$Class_code %in% c("c","k"),"Contained in reference",
                                  ifelse(tracking$Class_code == "j","Multi-exon with partial junction match",
                                         ifelse(tracking$Class_code == "=", "Exact match",
                                                ifelse(tracking$Class_code %in% c("p","u","s"), "Other","Intron retention"
                                                ))))
mean_len <- tracking %>%
  group_by(Class_code_sum) %>%
  summarise(mean = mean(len), sd = sd(len), mean_log2 = mean(log2(len)), n = n())

p1 <- ggplot(tracking, aes(x=log2(len), color=Class_code_sum)) +
  geom_density(size=2, alpha=0.5) +
  geom_vline(data = mean_len, aes(xintercept = mean_log2, color = Class_code_sum), size = 2) +
  geom_text(data = mean_len, aes(x = mean_log2, y = 0.4, label = round(mean,0)), 
            size = 4, color = "black") +
  theme_classic(base_size=14) +
  xlab("LOG2 Transcript length (bps)") +
  theme(panel.grid.major.y = element_line(colour = "grey90"),
        legend.position="none") +
  scale_color_manual(values=rev(mycolors), 
                     name="Category",
                     breaks=c("Exact match","Multi-exon with partial junction match",
                              "Intron retention","Contained in reference","Other"))

tracking$num_exons <- ifelse(tracking$num_exons >= 30, 30, tracking$num_exons)

p2 <- ggplot(tracking, aes(x=num_exons, color=Class_code_sum)) +
  geom_density(size=2, alpha=0.5) +
  theme_classic(base_size=14) +
  xlab("Number exons per transcript") +
  theme(panel.grid.major.y = element_line(colour = "grey90"),
        legend.position="none") +
  scale_color_manual(values=rev(mycolors), 
                     name="Category",
                     breaks=c("Exact match","Multi-exon with partial junction match",
                              "Intron retention","Contained in reference","Other"))

pdf(file = "general_characteristics_selected_class_codes.pdf", width = 15, height = 4)
multiplot(p1, p2, cols=2)
dev.off()

tracking_sum <- tracking %>%
  group_by(gene_id, Annotation) %>%
  dplyr::summarise(Number = n()) %>%
  ungroup() %>%
  select(gene_id, Annotation, Number) %>%
  pivot_wider(names_from = "Annotation", values_from = "Number", values_fill=0)

tracking_sum$group_annotated <- ifelse(tracking_sum$Annotated == 0, '0',
                                       ifelse(tracking_sum$Annotated ==1, '1',
                                              ifelse(tracking_sum$Annotated ==2, '2',
                                                     ifelse(tracking_sum$Annotated ==3, '3',
                                                            ifelse(tracking_sum$Annotated < 7, '4-6',"7+")))))
tracking_sum$group_novel <- ifelse(tracking_sum$Novel == 0, '0',
                                   ifelse(tracking_sum$Novel ==1, '1',
                                          ifelse(tracking_sum$Novel ==2, '2',
                                                 ifelse(tracking_sum$Novel ==3, '3',
                                                        ifelse(tracking_sum$Novel < 8, '4-7',
                                                               ifelse(tracking_sum$Novel < 12, '8-11',
                                                                      ifelse(tracking_sum$Novel < 16, '12-15',
                                                                             ifelse(tracking_sum$Novel < 26, '16-25',
                                                                                    ifelse(tracking_sum$Novel < 36, '26-35','36+')))))))))
tracking_sum %>%
  group_by(group_novel, group_annotated) %>%
  dplyr::summarise(sum = n()) %>%
  ungroup()%>%
  mutate(group_novel = factor(group_novel, levels=c("0","1","2","3","4-7","8-11","12-15","16-25",
                                                    "26-35","36+"))) %>%
  mutate(group_annotated = factor(group_annotated, levels=c("0","1","2","3","4-6","7+"))) %>%
  ggplot(aes(x=group_annotated, y=group_novel, fill=log2(sum))) +
  geom_tile() +
  theme_classic(base_size=14) +
  theme(legend.position = "bottom") +
  ylab("# of novel transcripts per gene") +
  xlab("# of annotated transcripts per gene") +
  scale_fill_gradient(low = "#ffb441", high = "#0b5971", na.value = "black")
ggsave(filename = "num_novel_vs_annotated_summary_symbol.pdf", width = 5, height = 3)