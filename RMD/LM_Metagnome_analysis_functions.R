library(ggpubr)
library(ggrepel)
library(patchwork)
library(tidyverse)
library(phyloseq)
library(ggplot2)
library(reshape2)
library(nlme)
library(compositions)
library(vegan)
library(ape)
library(dplyr)
library(psych)
library (gridExtra)
library(ecodist)
library(rstatix)
library(tidyr)
library(ggprism)
library(RColorBrewer)

import_buglist_phyloseq <- function(path_buglist, metadata, delimiter,tmm) {
  library(phyloseq)
  buglist <- read.csv(path_buglist,header = T, sep = delimiter, row.names = 1, check.names = F)
  buglist[is.na(buglist)] <- 0
  metadata[is.na(metadata)] <- "NA"
  rm(buglist_species)
  buglist_species <- extract_species(buglist)  
  if (tmm==F) {
    bug_norm <- buglist_species %>%  sweep(2,colSums(.),"/")*100
  } else if (tmm==T) {
    library(edgeR)
    bug_dge <- DGEList(buglist_species,lib.size = colSums(buglist_species))
    bug_dge <- calcNormFactors(bug_dge, method = "TMM") #between sample norm
    bug_norm <- as.data.frame(edgeR::cpm(bug_dge, log = F))
  }
  rm(taxa_table)
  taxa_Table <- extract_taxaTable(bug_norm)
  phy_species <- import_phyloseq(bug_norm,taxa_Table,metadata)
  return(phy_species)
}
# import_buglist_phyloseq <- function(path_buglist, metadata, delimiter) {
#   library(phyloseq)
#   buglist <- read.csv(path_buglist,header = T, sep = delimiter, row.names = 1, check.names = F)
#   buglist[is.na(buglist)] <- 0
#   metadata[is.na(metadata)] <- "NA"
#   rm(buglist_species)
#   buglist_species <- extract_species(buglist)
#   buglist_species <- buglist_species %>%  sweep(2,colSums(.),"/")
#   buglist_species <- buglist_species*100
#   rm(taxa_table)
#   taxa_Table <- extract_taxaTable(buglist_species)
#   phy_species <- import_phyloseq(buglist_species,taxa_Table,metadata)
#   return(phy_species)
# }
import_pwylist_phyloseq <- function(path_pwylist, metadata, delimiter){
  library('phyloseq')
  library(edgeR)
  pwylist <- read.csv(path_pwylist,header = T, sep = delimiter, row.names = 1, check.names = F)
  pwylist[is.na(pwylist)] <- 0
  flt_pwy <- filter_prevalence_abundance_dataframe(pwylist)
  pwy_dge <- DGEList(flt_pwy,lib.size = colSums(flt_pwy))
  pwy_dge <- calcNormFactors(pwy_dge, method = "TMM") #between sample norm
  pwy_tmm <- as.data.frame(edgeR::cpm(pwy_dge, log = F))
  metadata[is.na(metadata)] <- "NA"
  sampdata <- sample_data(metadata)
  ptable <- otu_table(pwy_tmm, taxa_are_rows = T)
  phy_protein <- phyloseq(ptable, sampdata)
  random_tree = rtree(ntaxa(phy_protein), rooted=TRUE, tip.label=taxa_names(phy_protein))
  phy_protein = merge_phyloseq(phy_protein, sampdata, random_tree)
  return(phy_protein)
}
extract_species <- function(buglist){
  to_look = "s__"
  #  for (i in row.names(buglist)) {
  #   if (str_detect(i, to_look) && (! str_detect(i, "c__"))) {
  #       if (!exists("extracted_taxa")){
  #       extracted_taxa <- data.frame(buglist[i,])
  #     }
  #     else if (exists("extracted_taxa")){
  #       extracted_taxa <- rbind(extracted_taxa, buglist[i,])
  #     }
  #   }
  # }
  for (i in row.names(buglist)) {
    if (grepl(to_look,i)) {
      if (!exists("extracted_taxa")){
        extracted_taxa <- data.frame(buglist[i,],check.names = F)
      }
      else if (exists("extracted_taxa")){
        extracted_taxa <- rbind(extracted_taxa, buglist[i,])
      }
    }
  }
  extracted_taxa[is.na(extracted_taxa)] <- 0
  return(extracted_taxa)
}
extract_taxaTable <- function(only_species){
  library(stringr)
  library(phyloseq)
  for (i in row.names(only_species)) {
    tax_list <- str_split(i,"\\|", simplify = T)
    kingdom <- str_split(tax_list[1], pattern = "k__", simplify = T)[2]
    phylum<- str_split(tax_list[2], pattern = "p__", simplify = T)[2]
    class<- str_split(tax_list[3], pattern = "c__", simplify = T)[2]
    order<- str_split(tax_list[4], pattern = "o__", simplify = T)[2]
    family<- str_split(tax_list[5], pattern = "f__", simplify = T)[2]
    genus <- str_split(tax_list[6], pattern = "g__", simplify = T)[2]
    specie <- str_split(tax_list[7], pattern = "s__", simplify = T)[2]
    
    if (!exists("taxa_table")){
      taxa_table <- data.frame(kingdom, phylum, class, order, family, genus, specie, stringsAsFactors = F)
    } 
    else if (exists("taxa_table")){
      taxa_table <- rbind(taxa_table, c(kingdom, phylum, class, order, family, genus, specie))
    }
  }
  row.names(taxa_table) <- row.names(only_species)
  return(taxa_table)
}
import_phyloseq <- function(only_species,tax,meta){
  library(phyloseq)
  tax = tax_table(tax)
  colnames(tax) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  taxa_names(tax) <- row.names(only_species)
  otable <- otu_table(only_species,taxa_are_rows = T)
  mdata <- sample_data(meta)
  phyObject <- phyloseq(otable, mdata, tax)
  random_tree = rtree(ntaxa(phyObject), rooted=TRUE, tip.label=taxa_names(phyObject))
  phyObject = merge_phyloseq(phyObject,random_tree)
  return(phyObject)
}
filter_prevalence_abundance_dataframe <- function(unfilter_dataframe, abundance = 0.00, prevalence = 0.1){
  rows_to_keep=c()
  number_of_samples=dim(unfilter_dataframe)[2]
  colsums_vector = colSums(unfilter_dataframe)
  for (i in 1:dim(unfilter_dataframe)[1]) {
    row_vector = unfilter_dataframe[i,]
    relabun_row_vector = row_vector/colsums_vector
    num_over_abundance = sum((relabun_row_vector > abundance) == TRUE, na.rm = TRUE)
    if (num_over_abundance/number_of_samples > prevalence) {
      rows_to_keep = c(rows_to_keep, i)
    }
  }
  filtered_dataframe <- unfilter_dataframe[rows_to_keep,]
  filtered_dataframe <- filtered_dataframe[colSums(abs(filtered_dataframe), na.rm = TRUE) > 0]
  return(filtered_dataframe)
}
residual_Eadjusted <- function(table){
  Residual_adjusted_table <- matrix(nrow=nrow(table), ncol=ncol(table))
  for (i in 1:nrow(table)){
    linear<-lm(as.numeric(table[i,])~as.numeric(table[1,]))
    Residual_adjusted_table[i,]<-linear$residuals+ mean(as.numeric(table[i,]))
  }
  colnames(Residual_adjusted_table)<- colnames(table)
  rownames(Residual_adjusted_table)<-rownames(table)
  Residual_adjusted_table <- Residual_adjusted_table[-1,]
  return(Residual_adjusted_table)
}
density_Eadjusted <- function(table){
  Density_adjusted_table <- matrix(nrow=nrow(table), ncol=ncol(table))
  Density_adjusted_table <- table/table[1,]*100
  Density_adjusted_table <- Density_adjusted_table[-1,]
  return(Density_adjusted_table)
}
wilcoxon_calculate <- function(table1, table2){
  pval_wilcoxon <- matrix(ncol=1, nrow=nrow(table1))
  for (i in 1:nrow(table1)) {
    t <- wilcox.test(table1[i,],table2[i,])
    pval_wilcoxon[i,1] <- t$p.value
  }
  rownames(pval_wilcoxon) <- rownames(table1)
  colnames(pval_wilcoxon) <- "Wilcoxon pval"
  return(pval_wilcoxon)
}
icc_calculate <- function(table1, table2){
  icc.table <- matrix(ncol=1, nrow=nrow(table1))
  for (i in 1:nrow(table1)){
    corr_table1 <- data.frame(table1[i,], table2[i,])
    icc1<-ICC(corr_table1)
    icc.table[i,1]<-icc1$results$ICC[3]
  }  
  rownames(icc.table) <- rownames(table1)
  colnames(icc.table) <- "ICC number"
  return(icc.table)
}
## Format results
formatted_cor <- function(cor,pval,sig_level){
  res <- 
    cbind(
      # correlation coefficients
      cor %>% melt() %>% dplyr::rename("r" = value),
      # p-values
      pval %>% melt() %>% dplyr::rename("pval" = value) %>% dplyr::select(pval)
      # # number of observations -> Observed missing values
      # cor_res$n %>% melt() %>% rename("n" = value) %>% select(n)
    ) %>%
    mutate( sig_pval = ifelse(pval < sig_level, T, F),
            # pval_if_sig = ifelse(sig_pval, r, NA),
            r_if_sig = ifelse(sig_pval, r, NA),
            r=ifelse(sig_pval,r,0))
  
  return(res)
}
sp_cor <- function(otu,mdata,cor_method){
  # COMBINATION PAIRS
  tblcols <- expand.grid(row.names(otu), row.names(mdata))
  
  # COR TEST FUNCTION
  cfunc <- function(var1, var2) {
    cor.test(as.numeric(otu[var1,]), as.numeric(mdata[var2,]), method=cor_method,exact=FALSE)$p.value
  }
  
  # PEARSON MATRIX BUILD
  pval <- matrix(mapply(cfunc, tblcols$Var1, tblcols$Var2),
                 ncol = nrow(otu), nrow = nrow(mdata),
                 dimnames = list(rownames(mdata), rownames(otu)))
  
  cor <- cor(t(mdata),t(otu),method=cor_method)
  pval_adj <- matrix(NA,nrow = nrow(pval),ncol = ncol(pval),dimnames = dimnames(pval))
  pval_adj[] <- p.adjust(pval, method = "BH")
  
  res <- list(cor=cor,pval=pval,pval_adj=pval_adj)
  return(res)
}
filter_cor <- function(cor_table, pval_table, sig_level){
  cor_flt <- data.frame(matrix(NA, nrow = nrow(cor_table), ncol = ncol(cor_table),dimnames = dimnames(cor_table)),check.names = F)
  
  for (i in 1:nrow(pval_table)) {
    for (j in 1:ncol(pval_table)) {
      if (pval_table[i,j]>sig_level) {
        cor_flt[i,j]=NA
      } else {
        cor_flt[i,j]=cor_table[i,j]
      }
    }
  }
  
  ind <- apply(cor_flt, 1, function(x) all(is.na(x)))
  cor_flt <- cor_flt[ !ind, ]
  ind2 <- apply(cor_flt, 2, function(x) all(is.na(x)))
  cor_flt <- cor_flt[, !ind2 ]
  
  pval_flt <- data.frame(pval_table[rownames(pval_table)%in%rownames(cor_flt),colnames(pval_table)%in%colnames(cor_flt)],check.names = F)
  
  res <- list(cor_flt=cor_flt,pval_flt=pval_flt)
}
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
makeStars <- function(x){
  stars <- c("****", "***", "**", "*", "")
  vec <- c(0, 0.0001, 0.001, 0.01, 0.05, 1)
  i <- findInterval(x, vec)
  stars[i]
}
get_core <- function(table, abundance = 0.001, prevalence = 0.9){
 
 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 # FUNCTION TO get core taxa/function
 #
 # parameters:
 #   - table: abundance table
 #   - abundance: minimum relative abundance, default is 0.001
 #   - prevalence: minimum prevalence, default is 0.9
 #
 # returns: three tables: shannon index, chao1 index, observed species numbers
 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
 rows_to_keep=c()
 number_of_samples=dim(table)[2]
 colsums_vector = colSums(table)
 for (i in 1:dim(table)[1]) {
  row_vector = table[i,]
  relabun_row_vector = row_vector/colsums_vector
  num_over_abundance = sum((relabun_row_vector > abundance) == TRUE, na.rm = TRUE)
  if (num_over_abundance/number_of_samples > prevalence) {
   rows_to_keep = c(rows_to_keep, i)
  }
  table[i,"prevalence"] <- num_over_abundance/number_of_samples*100
 }
 filtered_dataframe <- table[rows_to_keep,]
 return(list(filtered_dataframe, rows_to_keep,table))
}

plot_high_abundance <- function (physeq, x = "Sample", y = "Abundance", fill = NULL, 
                                 title = NULL, lab = FALSE, ab=1, legend = c("right", "bottom")) {
 
 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 # FUNCTION TO CREATE A TAXONOMY BARPLOT FROM A PHYLOSEQ OBJECT EXCLUDING
 # LOW ABUNDANT GROUPS
 #
 # parameters:
 #   - physeq: phyloseq object
 #   - x: x axis variable. Default: "Sample"
 #   - y: y axis variable. Default: "Abundance"
 #   - fill: taxa level to display. For example: "Species"
 #   - title: plot title if needed
 #   - lab: display x axis labels. Default: FALSE
 #   - ab: relative abundance threshold. Default: 1.0
 #   - legend: legend position
 #
 # returns: ggplot barplot 
 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 palette_c25 <- c(
  "#bfd1d0","#f4d7de","#f4ceb8","#ebe1b9","#e2cbd9","#a4d4dc","#c2a2c2","#f3f3ab","#d7e5ec","#e9a7b8",
  "#fccfb8","#cabed8","#738b8a","#e0f0e3","#fcecc0","#da9076","#ffdada","#647c8b","#caa26b","#f4bbad",
  "#E3971F","#5A7C4D","#8B5357","#123B57","#3B6F7C", "#A27554",
  "#CAC3BD","#646F83","#A08D75","#CDA632","#4F7792","#3A530D",
  "#E2C48D","#88A0B5","#444349","#C4BBBE","#72939E","#C25E7B",
  "#9C918E","#887434","#D48B28","#636B83","#5B8BAF","#E5A93C",
  "#D7D6D7","#BB2649",
  "#CE9C9D", "#F6F3EE", "#DFD8AB","#EDECEB",
  "#ADB5BE","#EAE8EB",  "#FBEFE3", "#FCE4E2",
  "#FAE8E8", "#F8E6E6", "#FAE9E2","#D5A1A3",
  "#E9BDBE", "#FAF3EB", "#EAEEE0","#B5A28A", 
  "#D9E6EC", "#EBF6FA", "#A4C8D5", "#80ADBC",
  "#AED9EA", "#C0D3D8", "#E7E7DB",
  "#ECDCDC", "#EEE8E8", "#E9CCC4",
  "#EFDFD5", "#DEE9EB", "#F0EDE8", "#C2B7B1", "#E1D7CD",
  "#D7E0E5", "#D4C8B6", "#D9BD9C",
  "#CAB08B", "#DEC8BD", "#B4C6DC"
 ) #https://www.tadmint.com/blog/8-pastel-color-palettes-inspired-by-nature
 
 
 # palette_c25 <- c( "#E3971F","#5A7C4D","#8B5357","#123B57","#3B6F7C", "#A27554",
 # "#CAC3BD","#646F83","#A08D75","#CDA632","#4F7792","#3A530D",
 # "#E2C48D","#88A0B5","#444349","#C4BBBE","#72939E","#C25E7B",
 # "#9C918E","#887434","#D48B28","#636B83","#5B8BAF","#E5A93C",
 # "#D7D6D7","#BB2649")
 
 rel_ab <- function(x){
  # compute relative abundance
  if (sum(x) == 0){
   return(x)
  } else {
   return(100 * x/sum(x))
  }
 }
 
 # Transform to relative abundance
 AID_norm <- transform_sample_counts(physeq, rel_ab)
 
 #Compile taxa by rank (filtering out low abundance taxa)
 AID_Rank <- AID_norm  %>%
  tax_glom(taxrank = fill) %>%                     # agglomerate taxa at order level
  psmelt() %>%                                        # Melt phyloseq object to long format for producing graphics with ggplot2
  filter(Abundance > ab)                      # Filter out orders below 1% in each sample
 
 ## order taxa levels by abundance 
 
 taxa_ordered <- AID_Rank %>%  
  dplyr::select(Abundance, fill) %>% 
  group_by(get(fill)) %>% 
  summarise(Total=sum(Abundance)) %>% 
  ungroup() %>% 
  arrange(-Total) %>% 
  pull("get(fill)")
 
 AID_Rank[, fill] <- factor(AID_Rank[, fill], levels=taxa_ordered)
 
 
 p = ggplot(AID_Rank, aes_string(x = x, y = y, fill = fill))
 p = p + geom_bar(stat = "identity", position = position_stack(reverse = TRUE))
 p = p + theme_prism(base_fontface = "bold", base_line_size = 1, base_family = "Arial",base_size = 14)+ theme(
  legend.position = "right",
  axis.title.x = element_blank(),
  strip.text = element_text(size = 14),
  legend.spacing.x = unit(0, "pt"),
  legend.text = element_text(margin = margin(r = 20))
 ) + guides(fill = guide_legend(override.aes = list(size = 3)))+ scale_y_continuous(
  limits = c(0, 100),
  expand = c(0, 0),
  breaks = seq(0, 100, 25),
  guide = "prism_offset"
 ) + scale_fill_manual(values=rep(palette_c25, 10)) 


 if (lab == TRUE){
  p = p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
 } else {
  p = p + scale_x_discrete(labels = NULL)
 }
 if (!is.null(title)) {
  p <- p + ggtitle(title)
 }
 if (legend == "bottom"){
  p = p + theme(legend.position="bottom") + guides(fill = guide_legend(ncol=4))
 }
 
 return(p) 
}
plot_beta_diversity <- function(physeq, method, weighted=T, variable){
 
 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 # FUNCTION TO FILTER AN ABUNDANCE TABLE BASED ON RELATIVE ABUNDANCE AND PREVALANCE
 #
 # parameters:
 #   - physeq: phyloseq object
 #   - method: distance metric, "unifrac" or "bray"
 #   - weighted: when "unifrac" is chosen, you have to select weighted or no, default is TRUE
 #   - variable: the metadata variable that you want to use to color the data points
 #
 # returns: a filtered abundance table
 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
 if (method=="unifrac") {
  if (weighted == T) {
   wunifrac_dist <- phyloseq::distance(physeq, method="unifrac", weighted=T)
  } else if (weighted == F) {
   wunifrac_dist <- phyloseq::distance(physeq, method="unifrac", weighted=F)
  }
  ordination <- ordinate(physeq, method="PCoA", distance=wunifrac_dist)
 } else {
  ordination <- ordinate(physeq, method="PCoA", distance=method)
 }
 p <- plot_ordination(physeq, ordination, color=variable)
 
 return(p)
}
plot_alpha_diversity <- function(physeq, method, attribute){
 
 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 # FUNCTION TO CALCULATE ALPHA DIVERSITIES OF A PHYLOSEQ OBJECT
 #
 # parameters:
 #   - physeq: phyloseq object
 #   - method: alpha diversity index, "Shannon", "Chao1", "Obs"
 #   - attribute: the metadata variable that you want to use to color the data points
 #
 # returns: three tables: shannon index, chao1 index, observed species numbers
 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
 round_meta <- data.frame(sample_data(physeq))
 round_otu <- data.frame(otu_table(physeq))
 round_otu <- round_otu %>% mutate_if(is.numeric,round)
 round_otu <- round_otu[colSums(round_otu) > 0]
 round_otu <- round_otu[rowSums(round_otu) > 0,]
 otable <- otu_table(round_otu,taxa_are_rows = T)
 mdata <- sample_data(round_meta)
 round_taxa <- extract_taxaTable(round_otu) %>% as.matrix() %>% tax_table()
 phy_round <- phyloseq(otable,mdata,round_taxa)
 random_tree = rtree(ntaxa(phy_round), rooted=TRUE, tip.label=taxa_names(phy_round))
 phy_round = merge_phyloseq(phy_round, mdata, random_tree)
 
 diversity <- plot_richness(phy_round, measures = method)$data
 
 parse_var <- parse(text = attribute)[[1]]
 
 p <- ggplot(diversity,aes(x=reorder(eval(parse_var),value,na.rm=T),y=value)) + 
  geom_boxplot(width=0.5,lwd=1.5) +
  geom_jitter(width=0.3,aes(color=eval(parse_var))) +
  scale_color_brewer(palette="Dark2")+theme_classic()+
  theme(text = element_text(size=25),
        axis.text.x = element_text(angle=0, hjust=1)) +
  xlab(attribute) +
  ylab(method)

 return(p)
}
ND <- function(mat, beta = 0.99, alpha = 1, control = TRUE, linear_mapping_before = FALSE, linear_mapping_after = FALSE) {
  stopifnot(is.matrix(mat))
  stopifnot(nrow(mat) == ncol(mat))
  stopifnot(all(mat >= 0))
  stopifnot(beta > 0 & beta < 1)
  stopifnot(alpha > 0 & alpha <= 1)
  mat_max = max(mat)
  mat_min = min(mat)
  if (mat_max == mat_min) {
    stop("the input matrix is a constant matrix")
  }
  
  ################ preprocessing the input matrix ################
  
  # linearly mapping the input matrix to be between 0 and 1
  if (linear_mapping_before) {
    mat = (mat - mat_min) / (mat_max - mat_min)
  }
  # diagonal values are filtered, as 0
  diag(mat) <- 0.
  # filtered the edges
  y <- quantile(mat, 1 - alpha)
  mat_th <- mat
  mat_th[mat_th < y] <- 0.
  # making the matrix symetric if already not
  mat_th <- (mat_th + t(mat_th)) / 2
  
  ################ eigen decomposition ################
  eigen_fit <- eigen(mat_th)
  U <- eigen_fit$vectors
  D <- eigen_fit$values
  
  lam_n <- min(c(min(D), 0))
  lam_p <- max(c(max(D), 0))
  m1 <- lam_p * (1 - beta) / beta
  m2 <- lam_n * (1 + beta) / beta
  m <- max(m1, m2)
  
  ################ network deconvolution ################
  D <- D / (D + m)
  mat_new1 <- U %*% diag(D) %*% solve(U)
  
  ################ displying direct weights ################
  if (control) {
    m2 <- min(mat_new1)
    mat_new2 <- mat_new1 + max(-m2, 0)
  } else {
    ind_nonedges <- mat_th == 0.
    m1 <- max(mat * ind_nonedges)
    m2 <- min(mat_new1)
    mat_new2 <- mat_new1 + max(m1 - m2, 0)
    mat_new2[ind_nonedges] <- mat[ind_nonedges]
  }
  
  if (linear_mapping_after) {
    m1 <- min(mat_new2)
    m2 <- max(mat_new2)
    mat_nd <- (mat_new2 - m1) / (m2 - m1)
  } else {
    mat_nd <- mat_new2
  }
  
  return(mat_nd)
}
#https://github.com/luyiyun/Network-Deconvolution