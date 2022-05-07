library("tidyverse")
library("magrittr")
library("optparse")

#Parse command-line options
option_list <- list(
  make_option(c("--suppa_dir"), type="character", default=NULL,
              help="Directory of SUPPA output files without suffix", metavar = "type"))

message(" ## Parsing options")
opt <- optparse::parse_args(OptionParser(option_list=option_list))

suppa_dir = opt$suppa_dir

process_suppa_output <- function(dir, as_event) {
  myas <- read.table(paste0(dir,as_event,".ioe"), header=TRUE)
  myas$alternative_transcripts <- as.character(myas$alternative_transcripts)
  s <- strsplit(myas$alternative_transcripts, split = ",")
  myas <- data.frame(diffsplice_events = as_event,
                   seqname = rep(myas$seqname, sapply(s, length)),
                   gene_id = rep(myas$gene_id, sapply(s, length)),
                   event_id = rep(myas$event_id, sapply(s, length)),
                   total_transcripts = rep(myas$total_transcripts, sapply(s, length)),
                   inclusion = unlist(s))
  myas$total_transcripts <- as.character(myas$total_transcripts)
  s <- strsplit(myas$total_transcripts, split = ",")
  myas <- data.frame(diffsplice_events = as_event,
                   seqname = rep(myas$seqname, sapply(s, length)),
                   gene_id = rep(myas$gene_id, sapply(s, length)),
                   event_id = rep(myas$event_id, sapply(s, length)),
                   inclusion = rep(myas$inclusion, sapply(s, length)),
                   exlusion = unlist(s))
  myas$inclusion <- as.character(myas$inclusion)
  myas$exlusion <- as.character(myas$exlusion)
  myas <- subset(myas, inclusion != exlusion)
  myas <- reshape2::melt(myas, measure.vars = c("inclusion","exlusion"), variable.name = "Type", value.name = "transcript_id")
  myas %<>%
    unique() %>%
    arrange(seqname, gene_id, event_id, transcript_id, Type) %>%
    distinct(event_id, transcript_id, .keep_all= TRUE, fromLast=T) %>%
    dplyr::select(-fromLast)
  }

A3 <- process_suppa_output(suppa_dir, "A3")
A5 <- process_suppa_output(suppa_dir, "A5")
AL3 <- process_suppa_output(suppa_dir, "AL3")
AF5 <- process_suppa_output(suppa_dir, "AF5")
AF <- process_suppa_output(suppa_dir, "AF")
AL <- process_suppa_output(suppa_dir, "AL")
MX <- process_suppa_output(suppa_dir, "MX")
RI <- process_suppa_output(suppa_dir, "RI")
SE <- process_suppa_output(suppa_dir, "SE")

alt3utr <- read.table(paste0(suppa_dir, "alt3utr.events.tsv"), header=TRUE)
alt3utr$excluded_isoform <- as.character(alt3utr$excluded_isoform)
s <- strsplit(alt3utr$excluded_isoform, split = ",")
alt3utr <- data.frame(seqname = rep(alt3utr$seqnames, sapply(s, length)),
                      gene_id = rep(alt3utr$feature_id, sapply(s, length)),
                   event_id = rep(alt3utr$coordinate, sapply(s, length)),
                   included_isoform = rep(alt3utr$included_isoform, sapply(s, length)),
                   excluded_isoform = unlist(s))

alt3utr$included_isoform <- as.character(alt3utr$included_isoform)
s <- strsplit(alt3utr$included_isoform, split = ",")
alt3utr <- data.frame(diffsplice_events = "A3UTR",
                      seqname = rep(alt3utr$seqname, sapply(s, length)),
                   gene_id = rep(alt3utr$gene_id, sapply(s, length)),
                   event_id = rep(alt3utr$event_id, sapply(s, length)),
                   inclusion = unlist(s),
                   exlusion = rep(alt3utr$excluded_isoform, sapply(s, length)))
alt3utr <- unique(reshape2::melt(alt3utr, measure.vars = c("inclusion","exlusion"), variable.name = "Type", value.name = "transcript_id"))

alt5utr <- read.table(paste0(suppa_dir, "alt5utr.events.tsv"), header=TRUE)
alt5utr$excluded_isoform <- as.character(alt5utr$excluded_isoform)
s <- strsplit(alt5utr$excluded_isoform, split = ",")
alt5utr <- data.frame(seqname = rep(alt5utr$seqnames, sapply(s, length)),
                      gene_id = rep(alt5utr$feature_id, sapply(s, length)),
                      event_id = rep(alt5utr$coordinate, sapply(s, length)),
                      included_isoform = rep(alt5utr$included_isoform, sapply(s, length)),
                      excluded_isoform = unlist(s))

alt5utr$included_isoform <- as.character(alt5utr$included_isoform)
s <- strsplit(alt5utr$included_isoform, split = ",")
alt5utr <- data.frame(diffsplice_events = "A5UTR",
                      seqname = rep(alt5utr$seqname, sapply(s, length)),
                      gene_id = rep(alt5utr$gene_id, sapply(s, length)),
                      event_id = rep(alt5utr$event_id, sapply(s, length)),
                      inclusion = unlist(s),
                      exlusion = rep(alt5utr$excluded_isoform, sapply(s, length)))
alt5utr <- unique(reshape2::melt(alt5utr, measure.vars = c("inclusion","exlusion"), variable.name = "Type", value.name = "transcript_id"))

all_events <- rbind(A3, A5, AF, AL, MX, RI, SE, alt3utr, alt5utr)
write.table(all_events, paste0(suppa_dir,"_all_events.tsv"), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
