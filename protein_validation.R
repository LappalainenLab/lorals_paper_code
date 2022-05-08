library('tidyverse')
library('magrittr')
require('RColorBrewer')

#Generates Figure 2C and extended figure 5A-D

# Prepare metadata table
samps <- read.table("data/metadata.txt",
                    header=TRUE, fill = TRUE, sep = "\t")
samps <- samps %>%
  filter(name != "GTEX_T5JC_brainfrontalcortex") %>%
  filter(!grepl("CVD", name)) %>%
  filter(tissue != "K562") %>%
  filter(tissue != "Cells - Cultured fibroblasts") %>%
  filter(!grepl("rep", name))
rownames(samps) <- samps$name

## Read in tracking info to determine novel transcripts
gff_class <- read_tsv("analysis/flair/new_noSJ_promoters.isoforms_annotation.txt")
gff_class_novel <- gff_class[gff_class$Class_code != "=",]
gff_class_annot <- gff_class[gff_class$Class_code == "=",]

# Prepare counts table
tpm_flair <- read_tsv("analysis/count_tables/trans_flair_reps.tpm.txt")
tpm_flair$transcript_id <- rownames(tpm_flair)
tpm_flair <- tpm_flair[tpm_flair$transcript_id %in% gff_class$fish, ]
gff_class <- gff_class[gff_class$fish %in% tpm_flair$transcript_id,]
gff_class <- gff_class[order(gff_class$fish),]
tpm_flair <- tpm_flair[order(tpm_flair$transcript_id),]

tpm_flair <- cbind(gff_class$gene_id, gff_class$transcript_id, gff_class$fish, gff_class$Class_code, tpm_flair[,-1])

colnames(tpm_flair)[1] <- "gene_id"
colnames(tpm_flair)[2] <- "transcript_id"
colnames(tpm_flair)[3] <- "fish"
colnames(tpm_flair)[4] <- "Class_code"
tpm_flair <- na.omit(tpm_flair)

### READ IN AND PROCESS PROTEIN DATA
read_protein <- function(tissue_name) {
  tissue <- read_csv(file = paste0("analysis/protein_valid/GTEx_",tissue_name,"_Newrarepep_proteinprophet_output_infosum.csv"), col_names = TRUE)
  tissue <- tissue %>%
    mutate(tissue = tissue_name) %>%
    separate(protein_name, sep = "\\|", into = c("protein_name","transcript_id","gene_id"), extra = "drop") %>%
    mutate(protein_name = gsub(".p1|.p2|.p3|.p4|.p5|.p6|.p7|.p8","",protein_name)) %>%
    filter(probability >= 0.99) %>%
    filter(group_sibling_id == "a")
  
  tissue_filter_protein1 <- list_by_tissue[[tissue_name]][gsub("\\..*","",list_by_tissue[[tissue_name]]$fish) %in% gsub("\\..*","",tissue_filter$protein_name),]
  tissue_filter_protein2 <- list_by_tissue[[tissue_name]][gsub("\\..*","",list_by_tissue[[tissue_name]]$transcript_id) %in% gsub("\\..*","",tissue_filter$transcript_id),]
  tissue_filter_protein <- unique(rbind(tissue_filter_protein1, tissue_filter_protein2))
}

## FIND transcripts expressed at least 5TPM
mytissues <- c("lung","liver","heart","muscle","brain","breast","pancreas")
class_codes_list <- list()
for (j in c(5,10,25,50,100)) {
  list_by_sample <- list()
  for (i in 5:57){
    list_by_sample[[i]] <- tpm_flair[,c(1:4,i)]
    list_by_sample[[i]] <- list_by_sample[[i]][list_by_sample[[i]][,5] >= 5,]
    names(list_by_sample)[[i]] <- colnames(list_by_sample[[i]])[5]
    list_by_sample[[i]] <- list_by_sample[[i]][,c(1:4)]
  }
  list_by_sample <- list_by_sample[lapply(list_by_sample,length)>0]
  list_by_tissue <- list()
  for (i in mytissues) {
    list_by_tissue[[i]] <- unique(do.call(rbind, list_by_sample[grepl(i, names(list_by_sample))]))
  }
  
  tissue_filter_protein <- list()
  for (i in mytissues) {
    tissue_filter_protein[[i]] <- read_protein(i)
  }
  
  #### LOOK INTO ALL  TOGETHER
  class_codes1 <- data.frame(table(list_by_tissue[["lung"]]$Class_code),
                             table(list_by_tissue[["liver"]]$Class_code),
                             table(list_by_tissue[["heart"]]$Class_code),
                             table(list_by_tissue[["muscle"]]$Class_code),
                             table(list_by_tissue[["brain"]]$Class_code),
                             table(list_by_tissue[["breast"]]$Class_code),
                             table(list_by_tissue[["pancreas"]]$Class_code))
  
  class_codes1 <- class_codes1[,-c(3,5,7,9,11,13)]
  colnames(class_codes1) <- c("Class_code",mytissues)
  
  class_codes1 <- class_codes1 %>%
    gather(Tissue, no_transcripts, -Class_code)
  
  class_codes2 <- data.frame(table(tissue_filter_protein$lung$Class_code),
                             table(tissue_filter_protein$liver$Class_code),
                             table(tissue_filter_protein$heart$Class_code),
                             table(tissue_filter_protein$muscle$Class_code),
                             table(tissue_filter_protein$brain$Class_code),
                             table(tissue_filter_protein$breast$Class_code),
                             table(tissue_filter_protein$pancreas$Class_code))
  
  class_codes2 <- class_codes2[,-c(3,5,7,9,11,13)]
  colnames(class_codes2) <- c("Class_code",mytissues)
  
  class_codes2 <- class_codes2 %>%
    gather(Tissue, no_transcripts, -Class_code)
  
  class_codes100 <- merge(class_codes1, class_codes2, by=c("Class_code","Tissue"), suffixes = c("_all","_match"))
  class_codes100 <- as_tibble(class_codes100)
  class_codes100 %<>%
    dplyr::filter(!Class_code %in% c("u","p")) %>%
    dplyr::mutate(Group = ifelse(Class_code == "=", "Annotated", "Novel")) %>%
    dplyr::group_by(Group, Tissue) %>%
    dplyr::summarise(no_transcripts_all1 = sum(no_transcripts_all), no_transcripts_match1 = sum(no_transcripts_match)) %>%
    ungroup() %>%
    unique() %>%
    mutate(Proportion = round(no_transcripts_match1/no_transcripts_all1*100, 3))
  class_codes_list[[j]] <- class_codes100
  ggplot(class_codes100, aes(x=Tissue, y=Proportion, fill=Group)) +
    geom_bar(stat="identity", position="dodge") +
    theme_classic(base_size=14) +
    scale_fill_manual(values = c("#56ba5d","#c1514e")) +
    ylab("% Protein validation rate") +
    coord_flip()
  ggsave("percent_validation_protein.pdf", width = 6, height=5)
}

class_codes_df <- bind_rows(class_codes_list, .id = "TPM")
class_codes_df$TPM <- ifelse(class_codes_df$TPM == 1, 5,
                             ifelse(class_codes_df$TPM == 2, 10,
                                    ifelse(class_codes_df$TPM == 3, 25,
                                           ifelse(class_codes_df$TPM == 4, 50, 100))))
class_codes_all <- data_summary(class_codes_all, varname=c("Proportion"),
                                groupnames=c("TPM","Group"))

ggplot(class_codes_all, aes(x=TPM, y=Proportion, color=Group)) +
  geom_point(size=3) +
  geom_line()+
  geom_errorbar(aes(ymin=Proportion-sd, ymax=Proportion+sd), width=.2,
                position=position_dodge(0.05)) +
  theme_classic(base_size=14) +
  facet_grid(Group~., scales = "free_y") +
  ylab("% Protein validation rate") +
  scale_color_manual(values = c("#56ba5d","#c1514e"))
ggsave("tpm_threshold_validation.pdf", width = 8, height=5)

## Decide which threshold to use and run this command and re-run with only that parameter
tissue_validated_transcripts <- list()
tissue_nonvalidated_transcripts <- list()
for (i in mytissues) {
  tissue_validated_transcripts[[i]] <- gff_class[gff_class$transcript_id %in% tissue_filter_protein[[i]]$transcript_id,]
  tissue_validated_transcripts[[i]]$Tissue <- i
  
  tissue_nonvalidated_transcripts[[i]] <-  gff_class[gff_class$transcript_id %in% list_by_tissue[[i]]$transcript_id,]
  tissue_nonvalidated_transcripts[[i]] <- tissue_nonvalidated_transcripts[[i]][!tissue_nonvalidated_transcripts[[i]]$fish %in% tissue_validated_transcripts[[i]]$fish,]
  tissue_nonvalidated_transcripts[[i]]$Tissue <- i
  
}

tissue_validated_transcripts <- unique(do.call(rbind, tissue_validated_transcripts))
tissue_validated_transcripts$Category <- "Validated"
validated_transcripts <- unique(tissue_validated_transcripts[,-8])

tissue_nonvalidated_transcripts <- unique(do.call(rbind, tissue_nonvalidated_transcripts))
tissue_nonvalidated_transcripts$Category <- "Non-validated"
nonvalidated_transcripts <- unique(tissue_nonvalidated_transcripts[,-8])
nonvalidated_transcripts <- nonvalidated_transcripts[!nonvalidated_transcripts$transcript_id %in% validated_transcripts$transcript_id,]
nonvalidated_transcripts$Category <- "Non-validated"

tissue_bothvalidated_transcripts <- rbind(tissue_validated_transcripts, tissue_nonvalidated_transcripts)

mt_by_tissue <- list()
for (i in mytissues) {
  if (i %in% c("breast","pancreas")) {
    mt_by_tissue[[i]] <- cbind(tpm_flair[,c(1:4)], tpm_flair[, grepl(i, colnames(tpm_flair))])
  }
  else {
  mt_by_tissue[[i]] <- cbind(tpm_flair[,c(1:4)], rowMeans(tpm_flair[, grepl(i, colnames(tpm_flair))]))
  }
  colnames(mt_by_tissue[[i]])[5] <- "TPM"
  mt_by_tissue[[i]]$Tissue <- i
}

#mt_by_tissue <- reduce(mt_by_tissue, full_join, by = c("gene_id","transcript_id","fish","Class_code","no_transcripts"))
mt_by_tissue <- bind_rows(mt_by_tissue, .id = "column_label")

tissue_bothvalidated_transcripts_tpm <- merge(tissue_bothvalidated_transcripts, mt_by_tissue,
                                              by.x=c("Tissue","gene_id","transcript_id","fish","Class_code"),
                                              by.y=c("Tissue","gene_id","transcript_id","fish","Class_code"))
ggplot(tissue_bothvalidated_transcripts_tpm, aes(x=log2(TPM), color=Category)) +
  geom_density(size=2) +
  geom_vline(xintercept = log2(5),linetype = "dashed",color="black") +
  theme_classic(base_size=14) +
  facet_grid(Annotation~.) +
  scale_color_manual(values = c("#7f8080","#fcb343"))
ggsave("expression_protein_validation.pdf", width = 6, height=6)

unique_valid_transcripts <- tissue_bothvalidated_transcripts %>%
  dplyr::select(-Tissue) %>%
  dplyr::filter(Category == "Validated") %>%
  dplyr::distinct() %>%
  dplyr::group_by(gene_id) %>%
  dplyr::mutate(Freq =  dplyr::n()) %>%
  filter(Freq > 1) %>%
  filter(grepl("Novel",toString(unique(Annotation)))) %>%
  ungroup() %>%
  arrange(gene_id, transcript_id)

unique_bothvalid_transcripts <- tissue_bothvalidated_transcripts %>%
  filter(transcript_id %in% unique_valid_transcripts$transcript_id) %>%
  arrange(gene_id, transcript_id)

write.table(unique_bothvalid_transcripts, "analysis/protein_valid/validated_multitranscript_genes.txt", col.names = 
              TRUE, row.names = FALSE, quote = FALSE, sep="\t")  

bothvalidated_transcripts <- rbind(validated_transcripts, nonvalidated_transcripts)
bothvalidated_transcripts$num_exons <- as.numeric(bothvalidated_transcripts$num_exons)
bothvalidated_transcripts$len <- as.numeric(bothvalidated_transcripts$len)
bothvalidated_transcripts$Annotation <- ifelse(bothvalidated_transcripts$Class_code=="=", "Annotated", "Novel")

class_codes2 <- bothvalidated_transcripts %>%
  dplyr::group_by(Category, Annotation) %>%
  dplyr::summarise(Freq =  dplyr::n())
M <- as.table(rbind(c(class_codes2[class_codes2$Category == "Non-validated",]$Freq),
                    c(class_codes2[class_codes2$Category == "Validated",]$Freq)))
tt <- chisq.test(M)

## Examine effect of splicing on protein validation
## Look into ALT TRANSCRIPT STRUCTURE EVENTS IN FLAIR 
event_char_flair <- read.table("analysis/diffsplice/new_noSJ_promoters_all_events.tsv", header=TRUE)
## Look into ALT TRANSCRIPT STRUCTURE EVENTS IN GENCODE QUANT 
event_char <- read.table("analysis/diffsplice/gencodev26_all_events.tsv", header=TRUE)

event_char_flair$diffsplice_events <- gsub("AF5","A5",event_char_flair$diffsplice_events)
event_char_flair$diffsplice_events <- gsub("AL3","A3",event_char_flair$diffsplice_events)
tracking_char_events <- merge(tissue_bothvalidated_transcripts, event_char_flair, by=c("transcript_id","gene_id"))

tracking_char_events$event_annotation <- ifelse(tracking_char_events$event_id %in% event_char$event_id, "annotated", "novel")

## Look into ALT SPLICING EVENTS
event_char_sum_bg <- tracking_char_events %>%
  dplyr::filter(!Tissue %in% c("pancreas","breast")) %>%
  dplyr::filter(Type == "inclusion") %>%
  dplyr::mutate(Tissue="All") %>%
  dplyr::select(gene_id, Tissue, event_id, diffsplice_events) %>%
  unique() %>%
  dplyr::group_by(diffsplice_events, Tissue) %>%
  dplyr::summarise(Freq=n()) %>%
  dplyr::ungroup()

event_char_sum <- tracking_char_events %>%
  dplyr::filter(!Tissue %in% c("pancreas","breast")) %>%
  dplyr::filter(Type == "inclusion") %>%
  filter(Category == "Validated") %>%
  dplyr::mutate(Tissue="All") %>%
  dplyr::select(gene_id, Tissue, event_id, diffsplice_events) %>%
  unique() %>%
  dplyr::group_by(diffsplice_events, Tissue) %>%
  dplyr::summarise(Freq=n()) %>%
  dplyr::ungroup()

tracking_char_events_sum  <- merge(event_char_sum, event_char_sum_bg, suffixes = c("_test","_bg"), by=c("diffsplice_events", "Tissue"))

## Look into ALT SPLICING EVENTS across tissues
event_char_sum_bg <- tracking_char_events %>%
  dplyr::filter(!Tissue %in% c("pancreas","breast")) %>%
  dplyr::filter(Type == "inclusion") %>%
  dplyr::select(gene_id, Tissue, event_id, diffsplice_events) %>%
  unique() %>%
  dplyr::group_by(diffsplice_events, Tissue) %>%
  dplyr::summarise(Freq=n()) %>%
  dplyr::ungroup()

event_char_sum <- tracking_char_events %>%
  dplyr::filter(!Tissue %in% c("pancreas","breast")) %>%
  dplyr::filter(Type == "inclusion") %>%
  filter(Category == "Validated") %>%
  dplyr::select(gene_id, Tissue, event_id, diffsplice_events) %>%
  unique() %>%
  dplyr::group_by(diffsplice_events, Tissue) %>%
  dplyr::summarise(Freq=n()) %>%
  dplyr::ungroup()

tracking_char_events_sum_tissue <- merge(event_char_sum, event_char_sum_bg,
                                         suffixes = c("_test","_bg"), by=c("diffsplice_events","Tissue"))
tracking_char_events_sum_tissue <- tracking_char_events_sum_tissue[!tracking_char_events_sum_tissue$Tissue %in% c("Pancreas","Breast"),]

tracking_char_events_sum <- rbind(tracking_char_events_sum, tracking_char_events_sum_tissue)
ggplot(tracking_char_events_sum, aes(x=diffsplice_events, y=Freq_test/Freq_bg, fill=Tissue)) +
  geom_bar(position="dodge", stat = "identity") +
  theme_classic(base_size=14) +
  theme(panel.grid.major = element_line(color = "grey90")) +
  ylab("Proportion of AltTS events validated") +
  xlab("") +
  scale_fill_manual(values = c("grey50","#b674b1","#c4625d","#815375",
                               "#0072b5", "#1f854e"))
ggsave("altTS_protein_validation_tissue.pdf", width = 8, height=8)

### Repeat analysis by novel/anno
temp_bg <- tracking_char_events %>%
  dplyr::filter(Type == "inclusion") %>%
  dplyr::select(gene_id, event_id, event_annotation, diffsplice_events) %>%
  unique() %>%
  dplyr::group_by(diffsplice_events, event_annotation) %>%
  dplyr::summarise(Freq=n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(event_annotation) %>%
  dplyr::mutate(Prop=Freq/sum(Freq)) %>%
  dplyr::ungroup() 

temp_ts <- tracking_char_events %>%
  #dplyr::filter(!Tissue %in% c("pancreas","breast")) %>%
  dplyr::filter(Type == "inclusion") %>%
  dplyr::filter(Category=="Validated") %>%
  dplyr::select(gene_id, event_id, event_annotation, diffsplice_events) %>%
  unique() %>%
  dplyr::group_by(diffsplice_events, event_annotation) %>%
  dplyr::summarise(Freq=n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(event_annotation) %>%
  dplyr::mutate(Sum=sum(Freq)) %>%
  dplyr::mutate(Prop=Freq/sum(Freq)) %>%
  dplyr::ungroup() 

tracking_char_events_sum  <- merge(temp_ts, temp_bg, suffixes = c("_test","_bg"), by=c("diffsplice_events","event_annotation"))

ggplot(tracking_char_events_sum, aes(x=diffsplice_events, y=Freq_test/Freq_bg, fill=event_annotation)) +
  geom_bar(position="dodge", stat = "identity") +
  geom_text(aes(label=Freq_test), position=position_dodge(width=0.9), vjust=-0.25) +
  theme_classic(base_size=14) +
  theme(panel.grid.major = element_line(color = "grey90")) +
  ylab("Proportion of AltTS events validated") +
  xlab("") +
  scale_fill_manual(values = c("#56ba5d","#c1514e"))
ggsave("altTS_protein_validation.pdf", width = 8, height=8)

pvalues_prop <- tracking_char_events_sum %>%
  dplyr::select(diffsplice_events, Freq_test, Freq_bg) %>%
  dplyr::nest_by(diffsplice_events) %>%
  dplyr::mutate(value = list(prop.test(data$Freq_test, data$Freq_bg)$p.value)) %>%
  ungroup() %>%
  dplyr::select(diffsplice_events, value) %>%
  unnest(value)
