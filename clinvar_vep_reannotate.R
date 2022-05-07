library(dplyr)
library(seqminer)
library(readr)
library(tidyr)
library(purrr)
library(ggplot2)
library(forcats)

vep_severity <- read.table("data/vep_severity.txt")

gencode <- read_tsv("analysis/vep/clinvar.most_severe_gencode.txt", col_names = TRUE)
flair <- read_tsv("analysis/vep/clinvar.most_severe_flair.txt", col_names = TRUE)

vep_anno <- merge(gencode[,c(1,2,7)], flair[,c(1,2,7)], suffixes = c("_gencode","_flair"),by=c("Uploaded_variation","Location"))

#Data can be downloaded from here: https://www.ncbi.nlm.nih.gov/clinvar/docs/maintenance_use/#download (accessed on 2020/12/02)
clinvar_vcf <- read_tsv("data/clinvar.vcf.gz",
                        comment = "#", col_names = c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO"))
clinvar <- read_tsv("data/variant_summary.txt.gz", col_names = TRUE)
clinvar$CHR_POS <- ifelse(clinvar$Start == clinvar$Stop, paste0(clinvar$Chromosome,":",clinvar$Start), paste0("chr",clinvar$Chromosome,":",clinvar$Start,"-",clinvar$Stop))

clinvar_vcf_vep <- merge(clinvar_vcf, vep_anno, by.x="ID", by.y="Uploaded_variation")

tab_clinvar <- merge(clinvar_vcf_vep,clinvar, by.x=c("ID"), by.y=c("#AlleleID"))
tab_clinvar$ClinicalSignificance2 <- ifelse(grepl("Conflicting interpretations of pathogenicity",tab_clinvar$ClinicalSignificance), "Conflicting interpretations of pathogenicity",
                                            ifelse(grepl("Benign/Likely benign",tab_clinvar$ClinicalSignificance), "Benign/Likely benign",
                                                   ifelse(grepl("Benign",tab_clinvar$ClinicalSignificance), "Benign",
                                                          ifelse(grepl("Likely benign",tab_clinvar$ClinicalSignificance), "Likely benign",
                                                                 ifelse(grepl("Pathogenic/Likely pathogenic",tab_clinvar$ClinicalSignificance), "Pathogenic/Likely pathogenic",
                                                                        ifelse(grepl("Pathogenic",tab_clinvar$ClinicalSignificance), "Pathogenic",
                                                                               ifelse(grepl("Likely pathogenic",tab_clinvar$ClinicalSignificance), "Likely pathogenic",
                                                                                      ifelse(grepl("Uncertain significance",tab_clinvar$ClinicalSignificance), "Uncertain significance",
                                                                                             "Other"))))))))

tab_clinvar_filter_star <- tab_clinvar %>%
  group_by(ClinicalSignificance2, ReviewStatus) %>%
  summarise(counts_all = n()) 

tab_clinvar_reassigned_filter_star <- tab_clinvar %>%
  filter(Consequence_gencode != Consequence_flair) %>%
  group_by(ClinicalSignificance2, ReviewStatus) %>%
  summarise(counts_reassigned = n()) 

tab_clinvar_summary <- merge(tab_clinvar_filter_star, tab_clinvar_reassigned_filter_star)

toplot <- tab_clinvar_summary %>%
  mutate(Proportion = (counts_reassigned/counts_all)*100) %>%
  mutate(gold_stars = ifelse(ReviewStatus=="reviewed by expert panel",3,
                             ifelse(ReviewStatus=="criteria provided, multiple submitters, no conflicts",2,
                                    ifelse(ReviewStatus=="no assertion criteria provided",0,1)))) %>%
    mutate(name = fct_relevel(ClinicalSignificance2, 
                            "Uncertain significance", "Likely benign", "Benign/Likely benign", 
                            "Benign", "Conflicting interpretations of pathogenicity", "Likely pathogenic", 
                            "Pathogenic/Likely pathogenic", "Pathogenic","Other"))

ggplot(toplot, aes(x=name, y=Proportion, fill=as.factor(gold_stars), label=counts_reassigned)) +
  geom_bar(position = "dodge", stat="identity") +
  theme_classic(base_size=14) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Percentage of re-assigned variants") + xlab("") +
  scale_fill_brewer(palette = "Dark2") +
  geom_text(
    aes(label = counts_reassigned),
    size = 3, vjust = 1.5, position = position_dodge(.9))
ggsave(filename = "clivar_reassigned.pdf", height = 6, width = 6)
ggsave("clinvar_flair_annot.pdf", height = 7.5, width = 8)