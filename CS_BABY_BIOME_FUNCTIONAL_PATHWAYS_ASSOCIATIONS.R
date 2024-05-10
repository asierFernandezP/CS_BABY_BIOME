# Title: "Humann_functions+associations_phenotypes"
# Author: "Trishla Sinha"
# Date: "14/03/2024"
# Last update: "14/03/2024"


#Loading packages
rm (list = ls())

#load packages 
library(tidyverse)
library(stringr)
library(vegan)
library(RColorBrewer)
library(wesanderson)
library(reshape2)
library(dplyr)
library(ggpubr)
library(wesanderson)
library(ggpubr)
library(lmerTest)
library(plyr)
library(pheatmap)

# Load functions

prepCleanHumann <- function(inPath, dropUnintegrated = T, dropUnmapped = T, dropTaxonSpecific = T,
                            presenceFilter = -1,minRelativeAbundance = -1, rescaleTaxa = T, novogeneIdsClean = T) {
  inDF <- read.table(inPath,sep='\t',header=T,quote ='',comment.char = '')
  #fix pathway ID (these tend to be weird coming out of humann)
  rownames(inDF) <- inDF$X..Pathway
  inDF$X..Pathway <- NULL
  #drop "junk" from humann (unintegrated/unmapped data & taxon-specific pathways)
  if (dropTaxonSpecific) {
    inDF <- inDF[grep('\\|',rownames(inDF),invert = T),]
  }
  if (dropUnintegrated) {
    inDF <- inDF[grep('UNINTEGRATED',rownames(inDF),invert = T),]
  }
  if (dropUnmapped) {
    inDF <- inDF[grep('UNMAPPED',rownames(inDF),invert = T),]
  }
  rownames(inDF)[grep('PWY',rownames(inDF),invert = T)] <- paste0('PWY_',rownames(inDF)[grep('PWY',rownames(inDF),invert = T)])
  inDF <- as.data.frame(t.data.frame(inDF))
  # fix sample IDs, remove duplicates
  inDF$ID <- rownames(inDF)
  if (novogeneIdsClean) {
    print ('NOTE: doing cleaning of Novogene IDs (keeping format <AB>_<CD>_<EFG...>)')
    for (i in c(1:nrow(inDF))) {
      ss <- strsplit(inDF$ID[i],'_')[[1]]
      sss <- paste0(ss[1],'_',ss[2],'_',ss[3])
      inDF$ID[i] <- sss
    }
  }
  if (sum(duplicated(inDF$ID) > 0)) {
    print(paste('WARNING: found ',sum(duplicated(inDF$ID) > 0),'duplicates, dropping them!'))
  }
  inDF <- inDF[!duplicated(inDF$ID),]
  rownames(inDF) <- inDF$ID
  inDF$ID <- NULL
  # make sure columns are actually numbers (otherwise filter dies; NOTE: this should not be necessary, but... )
  for (c in colnames(inDF)) {inDF[[c]] <- as.numeric(inDF[[c]])}
  # clean, rescale and save
  inDFt2 <- filterHumannDF(inDF,presPerc = presenceFilter,minMRelAb = minRelativeAbundance,minMedRelAb = -1,rescale = T,minSum = 1,verbose = 
                             T)
  inDFt2 <- inDFt2[,colSums(inDFt2)!=0]
  inDFt2$ID <- row.names(inDFt2)
  print(paste('Done, returning ',nrow(inDFt2),'samples'))
  inDFt2
}


filterHumannDF <- function(inDF,presPerc = -1 ,minMRelAb = -1,minMedRelAb= -1, minSum = -1, rescale=T,verbose=T,type='MetaCyc') {
  
  nonPWYpwys <- c("ARG+POLYAMINE-SYN: superpathway of arginine and polyamine biosynthesis",
                  "CHLOROPHYLL-SYN: chlorophyllide a biosynthesis I (aerobic, light-dependent)",
                  "GLYCOLYSIS-E-D: superpathway of glycolysis and Entner-Doudoroff",
                  "GLYCOLYSIS-TCA-GLYOX-BYPASS: superpathway of glycolysis, pyruvate dehydrogenase, TCA, and glyoxylate bypass",
                  "GLYCOLYSIS: glycolysis I (from glucose 6-phosphate)",
                  "GLYOXYLATE-BYPASS: glyoxylate cycle",
                  "HEME-BIOSYNTHESIS-II: heme biosynthesis I (aerobic)",
                  "MANNOSYL-CHITO-DOLICHOL-BIOSYNTHESIS: protein N-glycosylation (eukaryotic, high mannose)",
                  "NAD-BIOSYNTHESIS-II: NAD salvage pathway II",                  
                  "REDCITCYC: TCA cycle VIII (helicobacter)",
                  "TCA-GLYOX-BYPASS: superpathway of glyoxylate bypass and TCA",
                  "TCA: TCA cycle I (prokaryotic)")
  
  colnames(inDF)[colnames(inDF) %in% nonPWYpwys] <- paste0('PWY_',colnames(inDF)[colnames(inDF) %in% nonPWYpwys])
  
  if (type=='MetaCyc') {
    nonPWYdf <- as.data.frame(inDF[,-grep('PWY',colnames(inDF))])
    cnsNonPWYdf <- colnames(inDF[colnames(inDF)[-grep('PWY',colnames(inDF))] ])
  } else if (type=='EC') {
    nonPWYdf <- as.data.frame(inDF[,-grep('^EC_',colnames(inDF))])
    cnsNonPWYdf <- colnames(inDF[colnames(inDF)[-grep('^EC_',colnames(inDF))] ])
  } else if (type=='RXN') {
    nonPWYdf <- as.data.frame(inDF[,-grep('RXN',colnames(inDF))])
    cnsNonPWYdf <- colnames(inDF[colnames(inDF)[-grep('RXN',colnames(inDF))] ])
  } else if (type=='PFAM') {
    nonPWYdf <- as.data.frame(inDF[,-grep('^PF[01]',colnames(inDF))])
    cnsNonPWYdf <- colnames(inDF[colnames(inDF)[-grep('^PF[01]',colnames(inDF))] ])
  } else if (type=='GO') {
    nonPWYdf <- as.data.frame(inDF[,-grep('^GO',colnames(inDF))])
    cnsNonPWYdf <- colnames(inDF[colnames(inDF)[-grep('^GO',colnames(inDF))] ])
  } else if (type=='KEGG') {
    nonPWYdf <- as.data.frame(inDF[,-grep('^K[012]',colnames(inDF))])
    cnsNonPWYdf <- colnames(inDF[colnames(inDF)[-grep('^K[012]',colnames(inDF))] ])
  }
  colnames(nonPWYdf) <- cnsNonPWYdf
  if (type=='MetaCyc') {
    yesPWYdf <- as.data.frame(inDF[,grep('PWY',colnames(inDF))])
    cnsYesPWYdf <- colnames(inDF[colnames(inDF)[grep('PWY',colnames(inDF))] ])
  } else if (type=='EC') {
    yesPWYdf <- as.data.frame(inDF[,grep('^EC_',colnames(inDF))])
    cnsYesPWYdf <- colnames(inDF[colnames(inDF)[grep('^EC_',colnames(inDF))] ])
  } else if (type=='RXN') {
    yesPWYdf <- as.data.frame(inDF[,grep('RXN',colnames(inDF))])
    cnsYesPWYdf <- colnames(inDF[colnames(inDF)[grep('RXN',colnames(inDF))] ])
  } else if (type=='PFAM') {
    yesPWYdf <- as.data.frame(inDF[,grep('^PF[01]',colnames(inDF))])
    cnsYesPWYdf <- colnames(inDF[colnames(inDF)[grep('^PF[01]',colnames(inDF))] ])
  } else if (type=='GO') {
    yesPWYdf <- as.data.frame(inDF[,grep('^GO',colnames(inDF))])
    cnsYesPWYdf <- colnames(inDF[colnames(inDF)[grep('^GO',colnames(inDF))] ])
  } else if (type=='KEGG') {
    yesPWYdf <- as.data.frame(inDF[,grep('^K[012]',colnames(inDF))])
    cnsYesPWYdf <- colnames(inDF[colnames(inDF)[grep('^K[012]',colnames(inDF))] ])
  }
  
  # replaces NAs with 0s
  for (c in colnames(yesPWYdf)) {
    yesPWYdf[,c][is.na(yesPWYdf[,c])] <- 0.0
  }
  # rescale to rel ab (if rescale = T)
  if (rescale==T) {
    if (verbose) {print ('  >> rescaling')}
    rsums <- rowSums(yesPWYdf)
    rsums[rsums==0] <- 1.0
    yesPWYdf <- yesPWYdf/rsums
  }
  
  # filter for presence
  # -----------------
  nrRemoved = 0
  toRemove = c()
  for (c in colnames(yesPWYdf)) {
    nrnZ = as.numeric(sum(yesPWYdf[,c]!=0.0))
    if (nrnZ/as.numeric(nrow(yesPWYdf)) < presPerc) {
      nrRemoved = nrRemoved + 1
      toRemove <- c(toRemove,c)
    }
  }
  if (length(toRemove) > 0) {
    yesPWYdf <- yesPWYdf[,!(colnames(yesPWYdf) %in% toRemove)]
  }
  if (verbose) {print (paste(' > presence filter: Removed',nrRemoved,'pathways!, ',length(colnames(yesPWYdf)),'pathways left!')); }
  
  # filter for abundance (mean)
  # ---------------------------
  nrRemoved = 0
  toRemove = c()
  for (c in colnames(yesPWYdf)) {
    mn = mean(yesPWYdf[,c])
    if ( mn < minMRelAb) {
      nrRemoved = nrRemoved + 1
      toRemove <- c(toRemove,c)
    }
  }
  if (length(toRemove) > 0) {
    yesPWYdf <- yesPWYdf[,!(colnames(yesPWYdf) %in% toRemove)]
  }
  if (verbose) {print (paste(' > mean abundance filter: Removed',nrRemoved,'pathways!, ',length(colnames(yesPWYdf)),'pathways left!')); }
  
  # filter for abundance (median)
  # -----------------------------
  nrRemoved = 0
  toRemove = c()
  for (c in colnames(yesPWYdf)) {
    mn = median(yesPWYdf[,c])
    if ( mn < minMedRelAb) {
      nrRemoved = nrRemoved + 1
      toRemove <- c(toRemove,c)
    }
  }
  if (length(toRemove) > 0) {
    yesPWYdf <- yesPWYdf[,!(colnames(yesPWYdf) %in% toRemove)]
  }
  if (verbose) {print (paste(' > median abundance filter: Removed',nrRemoved,'pathways!, ',length(colnames(yesPWYdf)),'pathways left!')); }
  
  # do final rescale
  if (rescale==T) {
    if (verbose) {print ('  >> rescaling')}
    rsums <- rowSums(yesPWYdf)
    rsums[rsums==0] <- 1.0
    yesPWYdf <- yesPWYdf/rsums
  }
  inDF <- cbind.data.frame(nonPWYdf,yesPWYdf)
  if (verbose) {print ('> DONE')}
  inDF
}

mixed_models <- function(metadata, ID, CLR_transformed_data, pheno_list) {
  df <- metadata
  row.names(df) <- df[,ID]
  df<-merge(df, CLR_transformed_data, by='row.names')
  row.names(df) <- df$Row.names
  df$Row.names <- NULL
  
  Prevalent= c(colnames(CLR_transformed_data))
  #pheno_list= phenotypes
  
  Overall_result_phenos =tibble() 
  
  for (Bug in Prevalent){
    if (! Bug %in% colnames(df)){ next }
    #Prevalence = sum(as.numeric(as_vector(select(df, Bug)) > 0)) / dim(df)[1]
    # print (c(Bug, Prevalence))
    Bug2 = paste(c("`",Bug, "`"), collapse="")
    for ( pheno in pheno_list){
      pheno2 = paste(c("`",pheno, "`"), collapse="")
      df[is.na(df[colnames(df) == pheno]) == F, ID] -> To_keep
      df_pheno = filter(df, !!sym(ID) %in% To_keep )
      Model0 = as.formula(paste( c(Bug2,  " ~ read_depth + DNA_concentration_ng_ul + Timepoint_numeric + (1|CS_BABY_BIOME_ID)"), collapse="" )) 
      lmer(Model0, df_pheno) -> resultmodel0
      base_model=resultmodel0
      Model2 = as.formula(paste( c(Bug2,  " ~ read_depth + DNA_concentration_ng_ul + Timepoint_numeric + ",pheno2, "+ (1|CS_BABY_BIOME_ID)"), collapse="" ))
      lmer(Model2, df_pheno, REML = F) -> resultmodel2
      M = "Mixed"
      as.data.frame(anova(resultmodel2, base_model))['resultmodel2','Pr(>Chisq)']->p_simp
      as.data.frame(summary(resultmodel2)$coefficients)[grep(pheno, row.names(as.data.frame(summary(resultmodel2)$coefficients))),] -> Summ_simple
      Summ_simple %>% rownames_to_column("Feature") %>% as_tibble() %>% mutate(P = p_simp, Model_choice = M, Bug =Bug, Pheno=pheno, Model="simple") -> temp_output
      rbind(Overall_result_phenos, temp_output) -> Overall_result_phenos
    }
  }
  
  p=as.data.frame(Overall_result_phenos)
  p$FDR<-p.adjust(p$P, method = "BH")
  
  return(p)
  
}

mixed_models_cor_feeding <- function(metadata, ID, CLR_transformed_data, pheno_list) {
  df <- metadata
  row.names(df) <- df[,ID]
  df<-merge(df, CLR_transformed_data, by='row.names')
  row.names(df) <- df$Row.names
  df$Row.names <- NULL
  
  Prevalent= c(colnames(CLR_transformed_data))
  #pheno_list= phenotypes
  
  Overall_result_phenos =tibble() 
  
  for (Bug in Prevalent){
    if (! Bug %in% colnames(df)){ next }
    #Prevalence = sum(as.numeric(as_vector(select(df, Bug)) > 0)) / dim(df)[1]
    # print (c(Bug, Prevalence))
    Bug2 = paste(c("`",Bug, "`"), collapse="")
    for ( pheno in pheno_list){
      pheno2 = paste(c("`",pheno, "`"), collapse="")
      df[is.na(df[colnames(df) == pheno]) == F, ID] -> To_keep
      df_pheno = filter(df, !!sym(ID) %in% To_keep )
      Model0 = as.formula(paste( c(Bug2,  " ~ read_depth + DNA_concentration_ng_ul + Timepoint_numeric + feeding_mode + (1|CS_BABY_BIOME_ID)"), collapse="" )) 
      lmer(Model0, df_pheno) -> resultmodel0
      base_model=resultmodel0
      Model2 = as.formula(paste( c(Bug2,  " ~ read_depth + DNA_concentration_ng_ul + Timepoint_numeric + feeding_mode+",pheno2, "+ (1|CS_BABY_BIOME_ID)"), collapse="" ))
      lmer(Model2, df_pheno, REML = F) -> resultmodel2
      M = "Mixed"
      as.data.frame(anova(resultmodel2, base_model))['resultmodel2','Pr(>Chisq)']->p_simp
      as.data.frame(summary(resultmodel2)$coefficients)[grep(pheno, row.names(as.data.frame(summary(resultmodel2)$coefficients))),] -> Summ_simple
      Summ_simple %>% rownames_to_column("Feature") %>% as_tibble() %>% mutate(P = p_simp, Model_choice = M, Bug =Bug, Pheno=pheno, Model="simple") -> temp_output
      rbind(Overall_result_phenos, temp_output) -> Overall_result_phenos
    }
  }
  
  p=as.data.frame(Overall_result_phenos)
  p$FDR<-p.adjust(p$P, method = "BH")
  
  return(p)
  
}

# Loading phenotypic data 

metadata<-read.delim("~/Desktop/CS_Baby_Biome/submission/2024_submission/analysis/Metadata_EGA_CS_BABY_BIOME.txt")
metadata[sapply(metadata, is.character)] <- lapply(metadata[sapply(metadata, is.character)],  #convert character columns to factors
                                                   as.factor)
metadata$Timepoint_categorical<-factor(metadata$Timepoint_categorical, levels = c("W01", "W02", "W03", "W04", "W05", "W06", "M06", "M12"))
metadata$growth_p_limited<-factor(metadata$growth_p_limited, levels = c("P<10", "P10-P50", "P51-90", ">P90"))
metadata$feeding_mode<-factor(metadata$feeding_mode, levels = c("breast_feeding", "mixed_feeding", "formula_feeding"))

row.names(metadata)<-metadata$bioSampleId
# Loading and cleaning humaann pathways 

pathways <- prepCleanHumann("~/Desktop/CS_Baby_Biome/submission/2024_submission/analysis/PATHWAYS_HUMANN_2024/humann_merged_pathabundance.txt",
                        dropUnintegrated = T,dropUnmapped = T, dropTaxonSpecific = T,
                        presenceFilter = -1, minRelativeAbundance = -1, rescaleTaxa = T, novogeneIdsClean = F) 

row.names(pathways)<- substr(row.names(pathways), 0, 13) 
row.names (pathways)  <- str_replace(row.names (pathways)  , "_", "")
row.names (pathways)  <- str_replace(row.names (pathways)  , "kneadd", "")

names (pathways)
pathways$ID=NULL

pathways_filt <- filterHumannDF(pathways,presPerc = 0.20,minMRelAb = 0.001,minMedRelAb = -1,rescale = T,minSum = 1,verbose 
                                          = T)

# Selecting only CS Baby Biome samples 
common_row_names <- intersect(row.names(metadata), row.names(pathways_filt))
metadata_filtered <- metadata[common_row_names, ]
pathways_filt_CS <- pathways_filt[common_row_names, ]

# Transform data 
my_pseudocount_normal=min(pathways_filt_CS[pathways_filt_CS!=0])/2
pathways_filt_CLR<-decostand(pathways_filt_CS, "clr", pseudocount=my_pseudocount_normal)

metadata_filtered$ID<-row.names(metadata_filtered)

#Associations of phenotype with metacyc pathways
pathways_mixed_all <- mixed_models(metadata_filtered, "ID", pathways_filt_CLR, c("rand_AB",  "feeding_mode"))
setwd("/Users/trishlasinha/Desktop/CS_Baby_Biome/supplementary/tables/individual")
write.table(pathways_mixed_all, "pathway_associations_mixed_all_raw_15_03_2024.txt", sep="\t", row.names=F, quote = F)

# Association of AB with metacyc pathways after correcting for feeding 
pathways_mixed_cor_feeding <- mixed_models_cor_feeding(metadata_filtered, "ID", pathways_filt_CLR, c( "rand_AB"))

setwd("/Users/trishlasinha/Desktop/CS_Baby_Biome/supplementary/tables/individual")
write.table(pathways_mixed_cor_feeding, "pathway_associations_mixed_cor_feeding_15_03_2024.txt", sep="\t", row.names=F, quote = F)


