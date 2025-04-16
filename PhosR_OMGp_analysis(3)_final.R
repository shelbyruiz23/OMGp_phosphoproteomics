#Grubisha omgp data analysis
#PhosR data analysis
#################
#08/15/2023
#################

rm(list=ls()); gc()
options(stringsAsFactors = F)

#setwd
setwd("C:/Users/Shelby_macdonaldlab/OneDrive - University of Pittsburgh/MacDonaldLab_Data/Collaborations/Grubisha_lab/Phospho_OMGp")

library(tidyverse)
library(viridis)
library(readxl)
library(ggrepel)
library(ggplot2)
require(limma)
require(sva)
require(qvalue)

library(rstatix)
library(ggpubr)
library(dendextend)
library(dplyr)
library(data.table)
library(stringr)

###################
#install PhosR
#if(!require(devtools)){
#  install.packages("devtools") # If not already installed
# }

#devtools::install_github("PYangLab/PhosR",
#                          build_opts = c("--no-resave-data", "--no-manual"),
#                         
#                        build_vignettes = TRUE,
#                       
#                      dependencies = TRUE)

library(PhosR)
suppressPackageStartupMessages(library(SummarizedExperiment))
###################

load(paste0("PhosR/WorkData/PhosR_begin_final.RData"))

#factorize covariates
Plate=factor(colData(ppe)$Plate, levels=unique(colData(ppe)$Plate))
Time=factor(colData(ppe)$Time,levels=unique(colData(ppe)$Time))
Replicate=factor(colData(ppe)$Replicate,levels=unique(colData(ppe)$Replicate))
Replicate2=factor(colData(ppe)$Replicate2,levels=unique(colData(ppe)$Replicate2))

########################################################
########################################################
########################################################
# Gene centric analysis

suppressPackageStartupMessages({
  library(calibrate)
  library(limma)
  library(directPA)
  library(org.Rn.eg.db)
  library(reactome.db)
  library(annotate)
  library(PhosR)
})

data("PhosphoSitePlus")

#extract phosphosite information from the ppe object
sites = paste(sapply(GeneSymbol(ppe), function(x)x),";",
              sapply(Residue(ppe), function(x)x),
              sapply(Site(ppe), function(x)x),
              ";", sep = "")



# construct design matrix by group
design <- model.matrix(~ Time + 0)


# fit linear model for each phosphosite
fit <- lmFit(ppe@assays@data$normalised, design)


cont.matrix <- makeContrasts(
  TimeT30 - TimeT1,
  TimeT300 - TimeT1,
  TimeT900 - TimeT1,
  levels=design)

fit2 <- contrasts.fit(fit, cont.matrix)
fit3 <- eBayes(fit2)

#extract top-differentially regulated phosphosites from each condition
table.T30=topTable(eBayes(fit3),number = Inf,coef=1)
table.T300=topTable(eBayes(fit3),number = Inf,coef=2)
table.T900=topTable(eBayes(fit3),number = Inf,coef=3)

DE2.RUV <- c(sum(table.T30[,"adj.P.Val"] < 0.05),
             sum(table.T300[,"adj.P.Val"] < 0.05),
             sum(table.T900[,"adj.P.Val"] < 0.05))


############
###
#prepare the reactome annotation
pathways = as.list(reactomePATHID2EXTID)

path_names = as.list(reactomePATHID2NAME)
name_id = match(names(pathways), names(path_names))
names(pathways) = unlist(path_names)[name_id]

pathways = pathways[which(grepl("Rattus norvegicus", names(pathways), ignore.case = TRUE))]

pathways = lapply(pathways, function(path) {
  gene_name = unname(getSYMBOL(path, data = "org.Rn.eg"))
  toupper(unique(gene_name))
})
###

############
###prepare gene set that is to be tested for enrichment
#loop through each condition separately (T30,T300,T900)
#o <- rownames(table.T30)

#or do each condition separately
#using Foldchange as a proxy but you could also use p value


####
#Tc <- cbind(table.T30[o,"logFC"])
#Tc <- cbind(table.T300[o,"logFC"])
#Tc <- cbind(table.T900[o,"logFC"])

#OR
tmp=table.T30
tmp$sign=sign(tmp$logFC)
tmp$logP=-log10(tmp$P.Value)#adjusted or unadjusted p vlaue doesnt change anything
tmp$metric=tmp$logP/tmp$sign

o=rownames(tmp)

Tc <- cbind(tmp[o,"metric"])
####

#OR
#can also use the p value and weight of the log2FC but did not perform this for everyone
#tmp=table.T30
#tmp$logP=-log10(tmp$P.Value)
#tmp$metric=tmp$logP*tmp$logFC

#o=rownames(tmp)

#Tc <- cbind(tmp[o,"metric"])
####


rownames(Tc)=sites[match(o,rownames(ppe))]
rownames(Tc) <- gsub("(.*)(;[A-Z])([0-9]+)(;)", "\\1;\\3;", rownames(Tc))


colnames(Tc) <- c("T30")


#Summarize phosphosite-level information to proteins for the downstream gene-centric analysis.
Tc.gene <- phosCollapse(Tc, id=gsub(";.+", "", rownames(Tc)), 
                        stat=apply(abs(Tc), 1, max), by = "max")


#Overrepresentation analysis (fisher's exact test) & rank-based (gsea, wilcoxon rank-sum test)
#can do greater, less, or two.sided

path2 <- pathwayRankBasedEnrichment(Tc.gene[,1], 
                                    annotation=pathways, 
                                    alter = "two.sided")

path2 <- pathwayRankBasedEnrichment(Tc.gene[,1], 
                                    annotation=pathways, 
                                    alter = "greater")

path2 <- pathwayRankBasedEnrichment(Tc.gene[,1], 
                                    annotation=pathways, 
                                    alter = "less")

#T30.LogFC.two.sided=path2
#T30.LogFC.UP=path2
#T30.LogFC.DOWN=path2

#T300.LogFC.two.sided=path2
#T300.LogFC.UP=path2
#T300.LogFC.DOWN=path2

#T900.LogFC.two.sided=path2
#T900.LogFC.UP=path2
#T900.LogFC.DOWN=path2

###
#T30.p.two.sided=path2
#T30.p.UP=path2
#T30.p.DOWN=path2

#T300.p.two.sided=path2
#T300.p.UP=path2
#T300.p.DOWN=path2

#T900.p.two.sided=path2
#T900.p.UP=path2
#T900.p.DOWN=path2


#calculate GeneRatio = #DE genes (adj.Pval < 0.05)/total # genes in pathway

#for every time point make a list of DE genes
tmp.list=rownames_to_column(table.T30)#table.T900
tmp.list=tmp.list %>% separate(rowname, c("Protein","Gene","Site","AA.centered"),sep="~")
tmp.list=filter(tmp.list, adj.P.Val < 0.05)
tmp.list=tmp.list[,"Gene"]


#turn matrix of interest into df
tmp=as.data.frame(T30.p.two.sided)#T900.p.two.sided

tmp=rownames_to_column(tmp)
names(tmp)[1]="Pathway"
tmp$`# of substrates`=as.numeric(tmp$'# of substrates')


#identify gene ratios
tmp <- tmp %>%
  rowwise() %>%
  mutate(
    matching_genes = strsplit(substrates, ";") %>%
      unlist() %>%
      intersect(tmp.list) %>%
      length()
  ) %>%
  ungroup()

tmp$Gene.ratio=tmp$matching_genes/tmp$`# of substrates`



#T30.LogFC.two.sided=tmp
#T30.LogFC.UP=tmp
#T30.LogFC.DOWN=tmp

#T300.LogFC.two.sided=tmp
#T300.LogFC.UP=tmp
#T300.LogFC.DOWN=tmp

#T900.LogFC.two.sided=tmp
#T900.LogFC.UP=tmp
#T900.LogFC.DOWN=tmp

###
#T30.p.two.sided=tmp
#T30.p.UP=tmp
#T30.p.DOWN=tmp

#T300.p.two.sided=tmp
#T300.p.UP=tmp
#T300.p.DOWN=tmp

#T900.p.two.sided=tmp
#T900.p.UP=tmp
#T900.p.DOWN=tmp

#pull all dataframes into one excel file
library(openxlsx)

# Create a new Excel workbook
wb <- createWorkbook()

# List of dataframes
rankenrichment <- list(T30.LogFC.two.sided = T30.LogFC.two.sided,
                   T30.LogFC.UP = T30.LogFC.UP,
                   T30.LogFC.DOWN = T30.LogFC.DOWN,
                   T300.LogFC.two.sided = T300.LogFC.two.sided,
                   T300.LogFC.UP = T300.LogFC.UP,
                   T300.LogFC.DOWN = T300.LogFC.DOWN,
                   T900.LogFC.two.sided = T900.LogFC.two.sided,
                   T900.LogFC.UP = T900.LogFC.UP,
                   T900.LogFC.DOWN = T900.LogFC.DOWN,
                   T30.p.two.sided = T30.p.two.sided,
                   T30.p.UP = T30.p.UP,
                   T30.p.DOWN = T30.p.DOWN,
                   T300.p.two.sided = T300.p.two.sided,
                   T300.p.UP = T300.p.UP,
                   T300.p.DOWN = T300.p.DOWN,
                   T900.p.two.sided = T900.p.two.sided,
                   T900.p.UP = T900.p.UP,
                   T900.p.DOWN = T900.p.DOWN)

# Add each dataframe to the workbook
for (name in names(rankenrichment)) {
  addWorksheet(wb, sheetName = name)
  writeData(wb, sheet = name, x = rankenrichment[[name]], startCol = 1, startRow = 1)
}

# Save the workbook
saveWorkbook(wb, "PhosR/Results/files/rankenrichment.xlsx")


##############################################################
##############################################################
#make rank enrichment figure

#use T30.p.two.sided enrichment for figure

#select enriched pathway rows (p value < 0.05)
tmp=filter(T30.p.two.sided, pvalue < 0.05)

#remove rat from beginning of pathway name
tmp=tmp %>% 
  mutate(Pathway = sub("^Rattus norvegicus: ","",Pathway))

#order rows by Gene ratio
tmp=tmp[order(-tmp$Gene.ratio), ]

#keep the top 15 pathways
tmp=tmp[c(-16,-17), ]

names(tmp)[3]="number"
tmp$pvalue=as.numeric(tmp$pvalue)

#make figure



library(ggplot2)
library(forcats)
##theme function
theme_dose <- function(font.size=14) {
  theme_bw() +
    theme(axis.text.x = element_text(colour = "black",
                                     size = font.size, vjust =1 ),
          axis.text.y = element_text(colour = "black",
                                     size = font.size, hjust =1 ),
          axis.title = element_text(margin=margin(10, 5, 0, 0),
                                    color = "black",
                                    size = font.size),
          axis.title.y = element_text(angle=90)
    )
}

ggplot(tmp, showCategory = 10,
       
       aes(Gene.ratio, fct_reorder(Pathway, Gene.ratio))) +
  
  geom_segment(aes(xend=0, yend = Pathway)) +
  
  geom_point(aes(color=pvalue, size = number)) +
  
  #scale_color_gradientn(colours=c("#9D4EDD","#5A189A","#240046"),
                        scale_color_gradientn(colours=c("#9D4EDD","grey"),
                        
                        #trans = "log10",
                        
                        guide=guide_colorbar(reverse=TRUE, order=1)) +
  
  scale_size_continuous(range=c(2, 10)) +
  
  theme_dose(6) +
  
  xlab("Gene ratio") +
  xlim(0,0.55)+
  ylab(NULL) +
  
  
  
  ggtitle("T30sec Pathway enrichment")

#other colors
#=c("#f7ca64", "#46bac2", "#7e62a3"),
#=c("#9D4EDD","#5A189A","#240046")
#=c("#5ec962", "#21918c","#3b528b" )
##############################################################
##############################################################
##############################################################



