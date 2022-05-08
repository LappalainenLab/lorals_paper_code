#' Binomial test function
bt <- function(a, b, p = 0.5) {binom.test(a, b, 0.5, alternative=
                                            c("two.sided"), conf.level = 0.95)$p.value}

#' Process the output of calculate_asts --quant to get ASE
#' 
#' @param path The path to the table output of calculate_asts --quant.
#' @param min_reads_per_trans Minimum reads originating from either allele for a gene to be kept.
#' @param mytranscripts An SE object that matches the transcript name to a gene name. 
#' @return A processed table with p-values and FDR corrected statistics per gene.
#' 

process_ase_quant <- function(path, min_reads, mytranscripts){
  tab <- read_tsv(path, col_names = TRUE)
  tab$Gene <- mytranscripts@elementMetadata$GENCODE.ID[match(tab$transcript, mytranscripts@elementMetadata$TXNAME)]
  tab %<>%
    mutate(variantID = paste0(contig,"_",position,"_",refAllele,"_",altAllele)) %>%
    group_by(contig,position,refAllele,altAllele,Gene,variantID) %>%
    summarise(refCount = sum(refCount), altCount = sum(altCount)) %>%
    ungroup() %>%
    filter(!contig %in% c("chrM", "chrX", "chrY")) %>%
    mutate(Total_counts = refCount+altCount) %>%
    filter(Total_counts >= min_reads) %>%
    group_by(Gene) %>%
    arrange(Total_counts, .by_group = TRUE) %>%
    filter(row_number(Gene) == 1)  %>% # take the first row within each id
    ungroup() %>%
    mutate(refRatio = round(refCount/Total_counts, digits = 8)) %>%
    mutate(afc = round(log2(altCount/refCount), digits = 8)) %>%mutate(pvalue = mapply(bt, refCount, Total_counts)) %>%
    mutate(qvalue = qvalue(pvalue)$qvalue)
  return(tab)
}


#' Process the output of calculate_asts --length
#' 
#' @param path The path to the table output of calculate_asts --length
#' @param min_reads_ref Minimum reads for the REF allele
#' @param min_reads_alt Minimum reads for the ALT allele
#' @return A processed table with p-values and FDR corrected statistics.
#' 

process_asts_lengths <- function(path, min_reads_ref,min_reads_alt){
  tab <- read_tsv(path, col_names = c("Chrom","Position","Ref","Alt","Ref_count","Alt_count","D","pvalue"))
  tab <- tab %>%
    filter(!contig %in% c("chrM", "chrX", "chrY")) %>%
    mutate(variantID = paste0(Chrom,"_",Position,"_",Ref,"_",Alt)) %>%
    filter(Ref_count >= min_reads_ref | Alt_count >= min_reads_alt)
  tab <- na.omit(tab)
  if (length(tab$Ref_count) > 0)
    tab$FDR <- qvalue(tab$pvalue)$qvalue
  else
    tab$FDR <- NA
  return(tab)
}



#' Process the output of calculate_asts --quant
#' 
#' @param path The path to the table output of calculate_asts --quant.
#' @param min_reads_per_trans Minimum reads originating from either allele for a transcript to be kept.
#' @param min_total_reads Minimum reads in total for a transcript to be kept.
#' @param mytranscripts An SE object that matches the transcript name to a gene name. 
#' @return A processed table with p-values and FDR corrected statistics.
#' 

process_asts_quant <- function(path, min_reads_per_trans, min_total_reads, mytranscripts){
  tab <- read_tsv(path, col_names = TRUE)
  tab$Gene <- mytranscripts@elementMetadata$GENCODE.ID[match(tab$transcript, mytranscripts@elementMetadata$TXNAME)]
  tab %<>%
    filter(!contig %in% c("chrM", "chrX", "chrY")) %>%
    mutate(variantID = paste0(contig,"_",position,"_",refAllele,"_",altAllele)) %>%
    filter(refCount >= min_reads_per_trans | altCount >= min_reads_per_trans) %>%
    group_by(Gene, variantID) %>%
    filter(n() > 1) %>%
    mutate(by_reference = sum(refCount), by_alternative = sum(altCount)) %>%
    ungroup() %>%
    mutate(sum = by_reference + by_alternative) %>%
    filter(sum >= min_total_reads) %>%
    group_by(Gene) %>%
    filter(sum == max(sum)) %>%
    group_by(Gene, variantID) %>%
    filter(n() > 1) %>%
    ungroup() %>%
    dplyr::select(-by_reference, -by_alternative)  
  Mtab <- tab %>%
    gather(Allele, Count, -Gene, -variantID, -transcript, -contig, -position, -refAllele, -altAllele , -sum) %>%
    distinct()
  frequencies <- Mtab %>%
    group_by(Gene) %>%
    nest() %>%
    mutate(M = map(data, function(dat){
      dat2 <- dat %>% spread(transcript, Count)
      M <- as.matrix(dat2[, -c(1:7)])
      M[is.na(M)] <- 0
      row.names(M) <- dat2$Allele
      return(M)
    }))
  frequencies2 <- frequencies %>%
    mutate(pvalue = map_dbl(M, ~chisq.test(.x)$p.value)) %>%
    mutate(cohen = map_dbl(M, ~sqrt(chisq.test(.x)$statistic/sum(rowSums(.x))))) %>%
    dplyr::select(-data, -M) %>%
    ungroup() %>%
    mutate(FDR = qvalue(pvalue, lambda=0)$qvalue)
  tab$pvalue <- frequencies2$pvalue[match(tab$Gene,frequencies2$Gene)]
  tab$cohen <- frequencies2$cohen[match(tab$Gene,frequencies2$Gene)]
  tab$FDR <- frequencies2$FDR[match(tab$Gene,frequencies2$Gene)]
  retab2 <- tab %>%
    group_by(Gene, variantID) %>%
    mutate(transcript_number = length(unique(transcript))) %>%
    ungroup() %>%
    dplyr::select(Gene, variantID, sum, transcript_number, cohen, pvalue, FDR) %>%
    unique()
  return(tab)
}


#' Plot the transcript count and transcript ratios for a gene of interest
#' 
#' @param gene_name The HGNC name of the gene to plot
#' @return A ggplot 
#' 

plot_asts <- function(table, gene_name){
  gene_table <- table[table$Gene==gene_name,]
  gene_table$names <- gsub("GTEX_","",gene_table$names)
  gene_table_m <- gather(gene_table, Allele, Counts, refCount, altCount)
  gene_table_m$Significant <- ifelse(gene_table_m$FDR <= 0.05, "Yes", "No")
  p1 <- ggplot(gene_table_m, aes(x=Transcript, y=log2(Counts), fill=Allele, alpha=Significant)) +
    geom_bar(stat="identity", position = "dodge") +
    theme_classic(base_size=14) +
    facet_grid(variantID~names) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title.x = element_blank(),
          legend.position="bottom") +
    ylab("LOG2 Transcript Counts") +
    scale_fill_manual(name="Allele",
                      breaks=c("refCount", "altCount"),
                      labels=c("REF", "ALT"),
                      values=c("#155FA6","#F3723A")) +
    scale_alpha_manual(name="Significant",
                       values=c("No" = 0.5, "Yes" = 1))
  gene_table_m$Ratio <- gene_table_m$Counts/gene_table_m$sum
  p2 <- ggplot(gene_table_m, aes(x=Transcript, y=Ratio, fill=Allele, alpha=Significant)) +
    geom_bar(stat="identity", position = "dodge") +
    theme_classic(base_size=14) +
    facet_grid(~names) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title.x = element_blank()) +
    ylab("Transcript ratio") +
    scale_fill_manual(name="Allele",
                      breaks=c("refCount", "altCount"),
                      labels=c("REF", "ALT"),
                      values=c("#155FA6","#F3723A")) +
    scale_alpha_manual(name="Significant",
                       values=c("No" = 0.5, "Yes" = 1))
  
  g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)}
  
  mylegend<-g_legend(p1)
  
  p3 <- grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                                 p2 + theme(legend.position="none"),
                                 nrow=2),
                     mylegend, nrow=2,heights=c(10, 1),
                     top = gene_name)
  return(p1)
  #top = grid::textGrob(gene_name,gp=gpar(fontsize=20,font=3)))
}


#' Plot the ase count and transcript ratios for a gene of interest
#' 
#' @param table ase table output of function XX
#' @param gene_name The HGNC name of the gene to plot
#' @return A ggplot 
#' 
plot_ase <- function(table, gene_name) {
  gene_table <- table[,c("names", "variantID", "Gene", "pvalue", "qvalue", "Ref_counts", "Alt_counts")]
  temp <- gene_table[gene_table$Gene == gene_name,]
  temp <- gather(temp, Allele, Counts, Ref_counts, Alt_counts)
  temp$Significant <- ifelse(temp$qvalue <= 0.05, "Yes", "No")
  ggplot(temp, aes(x=variantID, y=log2(Counts), fill=Allele, alpha=Significant)) +
    geom_bar(stat="identity", position = "dodge") +
    theme_classic(base_size=14) +
    facet_grid(~names, scales = "free_x") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title.x = element_blank()) +
    ylab("LOG2 Gene counts") +
    ggtitle(gene_name) +
    scale_fill_manual(name="Allele",
                      breaks=c("Ref_counts", "Alt_counts"),
                      labels=c("REF", "ALT"),
                      values=c("#155FA6","#F3723A")) +
    scale_alpha_manual(name="Significant",
                       values=c("No" = 0.5, "Yes" = 1))
}



# Normalize ONT RNA data
#
# Args:
#   counts: gene-by-sample matrix
#   log
# Returns a matrix normalised by library size

tpm_ont <- function(counts, log=TRUE) {
  rate <- counts / sum(counts) * 1e6
  if (log == TRUE) {
    return(log2(rate+1))
  } else {
    return(rate)
  } 
}



#' Process the output of calculate_asts --quant to get ASE from FLAIR
#' 
#' @param path The path to the table output of calculate_asts --quant.
#' @param min_reads_per_trans Minimum reads originating from either allele for a gene to be kept.
#' @param mytranscripts An SE object that matches the transcript name to a gene name. 
#' @return A processed table with p-values and FDR corrected statistics per gene.
#' 
process_ase_quant_flair <- function(path, min_reads, gff_table){
  tab <- read_tsv(path, col_names = TRUE)
  tab$Transcript <- gff_table$transcript_id[match(tab$transcript, gff_table$fish)]
  tab$Gene <- gff_table$gene_id[match(tab$transcript, gff_table$fish)]
  tab %<>%
    mutate(variantID = paste0(contig,"_",position,"_",refAllele,"_",altAllele)) %>%
    group_by(contig,position,refAllele,altAllele,Gene,variantID) %>%
    summarise(refCount = sum(refCount), altCount = sum(altCount)) %>%
    ungroup() %>%
    filter(!contig %in% c("chrM","chrX","chrY")) %>%
    mutate(Total_counts = refCount+altCount) %>%
    filter(Total_counts >= min_reads) %>%
    group_by(Gene) %>%
    arrange(Total_counts, .by_group = TRUE) %>%
    filter(row_number(Gene) == 1)  %>% # take the first row within each id
    ungroup() %>%
    mutate(refRatio = round(refCount/Total_counts, digits = 8)) %>%
    mutate(afc = round(log2(altCount/refCount), digits = 8)) %>%
    mutate(pvalue = mapply(bt, refCount, Total_counts)) %>%
    mutate(qvalue = qvalue(pvalue)$qvalue)
  #write.table(tab[,-c(1:4)], paste0(path,"_results.txt"), quote = FALSE, col.names = TRUE, row.names = FALSE, sep="\t")
  return(tab)
}


#' Process the output of calculate_asts --quant from FLAIR
#' 
#' @param path The path to the table output of calculate_asts --quant.
#' @param min_reads_per_trans Minimum reads originating from either allele for a transcript to be kept.
#' @param min_total_reads Minimum reads in total for a transcript to be kept.
#' @param mytranscripts An SE object that matches the transcript name to a gene name. 
#' @return A processed table with p-values and FDR corrected statistics.
#' 

process_asts_quant_flair <- function(path, min_reads_per_trans, min_total_reads, gff_table){
  tab <- read_tsv(path, col_names = TRUE)
  tab$Transcript <- gff_table$transcript_id[match(tab$transcript, gff_table$fish)]
  tab$Gene <- gff_table$gene_id[match(tab$transcript, gff_table$fish)]
  tab %<>%
    filter(!contig %in% c("chrM", "chrX", "chrY")) %>%
    dplyr::select(-transcript) %>%
    mutate(variantID = paste0(contig,"_",position,"_",refAllele,"_",altAllele)) %>%
    filter(refCount >= min_reads_per_trans | altCount >= min_reads_per_trans) %>%
    group_by(Gene, variantID) %>%
    filter(n() > 1) %>%
    mutate(by_reference = sum(refCount), by_alternative = sum(altCount)) %>%
    ungroup() %>%
    mutate(sum = by_reference + by_alternative) %>%
    filter(sum >= min_total_reads) %>%
    group_by(Gene) %>%
    filter(sum == max(sum)) %>%
    group_by(Gene, variantID) %>%
    filter(n() > 1) %>%
    ungroup() %>%
    dplyr::select(-by_reference, -by_alternative)  
  Mtab <- tab %>%
    gather(Allele, Count, -Gene, -variantID, -Transcript, -contig, -position, -refAllele, -altAllele , -sum) %>%
    distinct()
  frequencies <- Mtab %>%
    group_by(Gene) %>%
    nest() %>%
    mutate(M = map(data, function(dat){
      dat2 <- dat %>% spread(Transcript, Count)
      M <- as.matrix(dat2[, -c(1:7)])
      M[is.na(M)] <- 0
      row.names(M) <- dat2$Allele
      return(M)
    }))
  frequencies2 <- frequencies %>%
    mutate(pvalue = map_dbl(M, ~chisq.test(.x)$p.value)) %>%
    mutate(cohen = map_dbl(M, ~sqrt(chisq.test(.x)$statistic/sum(rowSums(.x))))) %>%
    dplyr::select(-data, -M) %>%
    ungroup() %>%
    mutate(FDR = qvalue(pvalue, lambda=0)$qvalue)
  tab$cohen <- frequencies2$cohen[match(tab$Gene,frequencies2$Gene)]
  tab$pvalue <- frequencies2$pvalue[match(tab$Gene,frequencies2$Gene)]
  tab$FDR <- frequencies2$FDR[match(tab$Gene,frequencies2$Gene)]
  retab2 <- tab %>%
    group_by(Gene, variantID) %>%
    mutate(transcript_number = length(unique(Transcript))) %>%
    ungroup() %>%
    dplyr::select(Gene, variantID, sum, transcript_number, cohen, pvalue, FDR) %>%
    unique()
  return(tab)
}


# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


#' Remove rows from dataframe if they contain NA in specific column
#' 
#' @param data A dataframe
#' @param desiredCols The Hname of the columns you want to filter on
#' @return A dataframe 
#' 

completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}


# Save pheatmap function
#
save_pheatmap_pdf <- function(x, filename, width=12.9, height=8.4) {
  pdf(filename, width = width, height = height, useDingbats = FALSE)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}


#' Perform fisher's test using four columns of a dataframe and output the results for each row
#' 
#' @param df The dataframe
#' @param col1 Column 1
#' @param col2 Column 2
#' @param col3 Column 3
#' @param col4 Column 4
#' @return A vector of fisher's t-test results
#' 

get_fisher <- function(df,col1,col2,col3,col4){
  mat <- matrix(as.numeric(df[c(col1,col2,col3,col4)]), ncol=2)
  f <- fisher.test(as.table(mat), alt="two.sided")
  return(f$p.value)}


# read files that follow a similar pattern and perform some predefined calculation
readIn <- function(table){
  tab <- read.delim(table, header=TRUE)
  tab[,table] <- log2(tab$ALT_COUNT/tab$REF_COUNT)
  return(tab)}
