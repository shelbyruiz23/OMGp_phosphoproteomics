#Grubisha omgp data analysis
#PhosR data setup
#################
#08/08/2023
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

#load(paste0("PhosR/WorkData/PhosR_setup_final.RData"))
load(paste0("PhosR/WorkData/PhosR_begin_final.RData"))

###################
# Calculate the minimum non-NA values required for each condition
min_values <- 3

#store rownames
df$rowname=meta.pep$PhosName

# Filter rows that meet the condition
filtered.df <- df %>%
  rowwise() %>%
  filter(
    sum(!is.na(c_across(starts_with("X1_")))) >= min_values &&
      sum(!is.na(c_across(starts_with("X30_")))) >= min_values &&
      sum(!is.na(c_across(starts_with("X300_")))) >= min_values &&
      sum(!is.na(c_across(starts_with("X900_")))) >= min_values
  ) %>%
  ungroup() #8088

filtered.df2=as.data.frame(filtered.df[,-25])
rownames(filtered.df2)=filtered.df$rowname

#fix meta.pep file
meta.pep=meta.pep %>%
  filter(PhosName %in% rownames(filtered.df2))

#save file
#write.csv(filtered.df2,"PhosR/Results/files/phosphosite_df_raw.csv")
#####
###
##
#

GeneSymbol=meta.pep$Gene
Residue=meta.pep$Rep_character
Site=meta.pep$Rep_numeric
Sequence=meta.pep$AA.centered



ppe <- PhosR::PhosphoExperiment(assays = list(Quantification = as.matrix(filtered.df2)), 
                                Site = Site, 
                                GeneSymbol = GeneSymbol, 
                                Residue = Residue, Sequence = Sequence)
#####
###
##
#



#####
df=S4Vectors::DataFrame(
  Time=meta.data$Time,
  Replicate=meta.data$Replicate,
  Plate=meta.data$Plate,
  Treatment=meta.data$Treatment,
  Replicate2=meta.data$Replicate2
)

rownames(df)=colnames(ppe)
SummarizedExperiment::colData(ppe)=df

ppe
#dim(ppe)
###################################################
###################################################
###################################################


#perform log2 transformation of the data
logmat <- log2(SummarizedExperiment::assay(ppe, "Quantification"))

# mark any missing values as NA
logmat[is.infinite(logmat)] <- NA
SummarizedExperiment::assay(ppe, "Quantification") <- logmat


#Perform filtering
#We first extract the grouping information 
grps <- paste0(colData(ppe)$Time)

#We filter for sites with at least 50% quantification rate (q ≥ 0.5) in one or more conditions
ppe <- selectGrps(ppe, grps, 0.5, n=1)


#Check the filtering results
ppe
dim(ppe) #8001 x 24


#Perform imputation of missing values
#PhosR enables site- and condition-specific imputation. Here, for each phosphosite in each condition, we impute its missing values in that condition (if any) using site- 
#and condition-specific imputation if the quantification rate within that condition is equal to or greater than a desired percentage (such as ≥ 50% in the example below).
set.seed(123)

ppe <- scImpute(ppe, 0.5, grps)

ppe

#Lastly, we can impute the remaining sites using tail-based imputation
ppe <- tImpute(ppe, assay = "imputed")


#save imputed df
#write.csv(ppe@assays@data$imputed,"PhosR/Results/files/phosphosite_df_imputed.csv")


boxplot(assay(ppe,"Quantification"))
boxplot(assay(ppe,"imputed"))
#############################################
#############################################
#############################################


#############################################
#############################################
#############################################
#We use limma package for calling for differentially phosphorylated sites 

#factorize covariates
Plate=factor(colData(ppe)$Plate, levels=unique(colData(ppe)$Plate))
Time=factor(colData(ppe)$Time,levels=unique(colData(ppe)$Time))
Replicate=factor(colData(ppe)$Replicate,levels=unique(colData(ppe)$Replicate))
Replicate2=factor(colData(ppe)$Replicate2,levels=unique(colData(ppe)$Replicate2))


###########
#limma remove batch effect 
design=model.matrix( ~ Time + 0)
batch_removal=removeBatchEffect(ppe@assays@data$imputed,
                                batch = Replicate,
                                batch2 = Plate,
                                design = design )

assay(ppe, "normalised") <- batch_removal
#write.csv(ppe@assays@data$normalised ,"PhosR/Results/files/residuals.csv")
######


####
######
#tmp
#can use instead of remove batch effect
#returns scaled residuals
design <- model.matrix(~ Plate + Replicate + 0)
fit <- lmFit(ppe@assays@data$imputed, design)

residuals.test=residuals(fit,ppe@assays@data$imputed)

#assay(ppe, "normalised") <- residuals.test
#write.csv(residuals.test,"PhosR/Results/files/scaled_residuals.csv")

###

#####
#######

#use limma package for calling for differentially phosphorylated sites 

design <- model.matrix(~ Time + 0)

# fit linear model for each phosphosite
fit <- lmFit(ppe@assays@data$normalised, design)


#define contrasts
cont.matrix <- makeContrasts(
  TimeT30 - TimeT1,
  TimeT300 - TimeT1,
  TimeT900 - TimeT1,
  levels=design)


fit2 <- contrasts.fit(fit, cont.matrix)

fit3 <- eBayes(fit2)


#plotSA(fit)

#save limma results
f=as.data.frame(fit3)
rownames(f)=rownames(ppe)
f=rownames_to_column(f)

colnames(f)[1]="PhosName"

#join f with other meta.pep information
limma.results=left_join(f,meta.pep,by="PhosName")

#save output files
#write.csv(f,"PhosR/Results/files/limmaresults.short.csv")
#write.csv(limma.results,"PhosR/Results/files/limmaresults.long.csv")


#save residual output as separate dataframe
residuals_with_preserved_time=as.data.frame(ppe@assays@data$normalised)
residuals_with_preserved_time=rownames_to_column(residuals_with_preserved_time)
colnames(residuals_with_preserved_time)[1]="PhosName"

#join residuals_with_preserved_time with other meta.pep information
residuals=left_join(residuals_with_preserved_time,meta.pep,by="PhosName")

#save output files
#write.csv(residuals_with_preserved_time,"PhosR/Results/files/residuals.short.csv")
#write.csv(residuals,"PhosR/Results/files/residuals.long.csv")

#############################################
#############################################
#############################################
#extract top-differentially regulated phosphosites from each condition
table.T30=topTable(eBayes(fit3),number = Inf,coef=1)
table.T300=topTable(eBayes(fit3),number = Inf,coef=2)
table.T900=topTable(eBayes(fit3),number = Inf,coef=3)


#write.csv(table.T30,"PhosR/Results/files/topTable.T30.csv")
#write.csv(table.T300,"PhosR/Results/files/topTable.T300.csv")
#write.csv(table.T900,"PhosR/Results/files/topTable.T900.csv")

DE2.RUV <- c(sum(table.T30[,"P.Value"] < 0.05),
             sum(table.T300[,"P.Value"] < 0.05),
             sum(table.T900[,"P.Value"] < 0.05))

#############################################
#############################################
#############################################
#check the grouping of the data throughout stepwise iterations

plotQC(SummarizedExperiment::assay(ppe,"Quantification"), 
       panel = "dendrogram", 
       grps = SummarizedExperiment::colData(ppe)$Time, 
       labels = colnames(ppe)) + 
  ggplot2::ggtitle("before batch correction")

plotQC(SummarizedExperiment::assay(ppe,"imputed"), 
       panel = "dendrogram", 
       grps = SummarizedExperiment::colData(ppe)$Time, 
       labels = colnames(ppe)) + 
  ggplot2::ggtitle("before batch correction")


plotQC(SummarizedExperiment::assay(ppe,"normalised"), 
       panel = "dendrogram", 
       grps = SummarizedExperiment::colData(ppe)$Time, 
       labels = colnames(ppe)) + 
  ggplot2::ggtitle("after batch correction")

###
plotQC(SummarizedExperiment::assay(ppe,"Quantification"), 
       panel = "pca", 
       grps = SummarizedExperiment::colData(ppe)$Time, 
       labels = colnames(ppe)) + 
  ggplot2::ggtitle("before batch correction")

plotQC(SummarizedExperiment::assay(ppe,"imputed"), 
       panel = "pca", 
       grps = SummarizedExperiment::colData(ppe)$Time, 
       labels = colnames(ppe)) + 
  ggplot2::ggtitle("before batch correction")


col=c("grey","#35b779","#31688e","#440154")

plotQC(SummarizedExperiment::assay(ppe,"normalised"), 
       panel = "pca", 
       grps = SummarizedExperiment::colData(ppe)$Time, 
       labels = colnames(ppe)) + 
  ggplot2::ggtitle("after batch correction") +
  theme_classic()+
  scale_color_manual(values=col)+
  theme(legend.position = "top", legend.justification = "right")

plotQC(SummarizedExperiment::assay(ppe,"Quantification"), 
       panel = "pca", 
       grps = SummarizedExperiment::colData(ppe)$Time, 
       labels = colnames(ppe)) + 
  ggplot2::ggtitle("before batch correction") +
  theme_classic()+
  scale_color_manual(values=col)+
  theme(legend.position = "top", legend.justification = "right")

#############################################
#############################################
#############################################
######################################
#plot p value distribution
tmp=limma.results[,c(16,17,18)]
tmp=gather(tmp,time,p,p.value.TimeT30...TimeT1:p.value.TimeT900...TimeT1)


#col=c("#35b779","#31688e","#440154")
#col=c("#2CB2D3","#5C5BBC","#8C04A4")
col=c("#9D4EDD","#5A189A","#240046")
#col=c("#1780A1","#5C4D7D","#892B64")

tmp %>%
  ggplot( aes(x=p,fill=time))+
  geom_histogram( color="black",alpha=0.7,position = 'identity')+
  theme_classic()+
  scale_fill_manual(values=col)+
  theme(legend.position = "top", legend.justification = "right")


#############################################
#############################################
#plot barplot of no. DE phosphopeptides 

#plot results by individual time point
#P.Value
tmp=data.frame(
  X=c("total", "DE 30 sec","DE 5 min", "DE 15 min"),
  Y=c(8001,1824,1685,880)
)



test=ggbarplot(
  tmp, x='X', y='Y',
  fill="black",
  stat="identity",
  width=0.9,
  position=position_dodge(0.9),
  orientation="vertical",
  order=c("total", "DE 30 sec","DE 5 min", "DE 15 min"),
  title="P.Value")+
  ylab("No. phosphopeptides")+
  xlab("")+
  geom_text(aes(label=Y),vjust=-0.2)+
  theme_classic()+
  theme(axis.text.x = element_text(angle=45,hjust=1),legend.position = c(0.6, 0.8))

plot(test)


#plot results by individual time point
#adj.P.Val
tmp=data.frame(
  X=c("total", "DE 30 sec","DE 5 min", "DE 15 min"),
  Y=c(8001,635,594,100)
)

test=ggbarplot(
  tmp, x='X', y='Y',
  fill='black',
  stat="identity",
  width=0.9,
  position=position_dodge(0.9),
  orientation="vertical",
  order=c("total", "DE 30 sec","DE 5 min", "DE 15 min"),
  title="adj.P.Value")+
  ylab("No. phosphopeptides")+
  xlab("")+
  geom_text(aes(label=Y),vjust=-0.2)+
  theme_classic()+
  theme(axis.text.x = element_text(angle=45,hjust=1),legend.position = c(0.6, 0.8))

plot(test)

#############################################
#############################################


###################################################
###################################################
###################################################
#make venn diagram of resulting DE phosphopeptides

#make one df with all values in it
tmp.T30=rownames_to_column(table.T30)
tmp.T30=tmp.T30[,c(1,5,6)]
names(tmp.T30)[1]="PhosName"
names(tmp.T30)[2]="P.Value.T30"
names(tmp.T30)[3]="adj.P.Val.T30"

tmp.T300=rownames_to_column(table.T300)
tmp.T300=tmp.T300[,c(1,5,6)]
names(tmp.T300)[1]="PhosName"
names(tmp.T300)[2]="P.Value.T300"
names(tmp.T300)[3]="adj.P.Val.T300"

tmp.T900=rownames_to_column(table.T900)
tmp.T900=tmp.T900[,c(1,5,6)]
names(tmp.T900)[1]="PhosName"
names(tmp.T900)[2]="P.Value.T900"
names(tmp.T900)[3]="adj.P.Val.T900"

tmp=merge(tmp.T30,tmp.T300,by= "PhosName" )
tmp=merge(tmp,tmp.T900,by= "PhosName" )

#identify overlaps
#X30=filter(tmp,P.Value.T30 < 0.05) #1824
#X30=X30[,1] 

#X300=filter(tmp,P.Value.T300 < 0.05) #1685
#X300=X300[,1]

#X900=filter(tmp,P.Value.T900 < 0.05) #880
#X900=X900[,1]


#identify overlaps
X30=filter(tmp,adj.P.Val.T30 < 0.05) #635
X30=X30[,1] 

X300=filter(tmp,adj.P.Val.T300 < 0.05) #594
X300=X300[,1]

X900=filter(tmp,adj.P.Val.T900 < 0.05) #100
X900=X900[,1]

#find phosphopeptides that are DE at every time point
overlap=intersect(X30, intersect(X300,X900))#227, 22

#remove overlapping peptides from every timepoint
X30=setdiff(X30,overlap)
X300=setdiff(X300,overlap)
X900=setdiff(X900,overlap)

X30.X300=intersect(X30,X300)#582, 214
X30.X900=intersect(X30,X900)#145, 12
X300.X900=intersect(X300,X900)#153, 20

#plotting weighted venn diagrams
#install.packages("venneuler")     # Install & load venneuler package
library("venneuler")

#A=1824-227-582-145
#B=1685-227-582-153
#C=880-227-145-153

A=635-22-214-12
B=594-22-214-20
C=100-22-12-20

plot(venneuler(c(A = A,          # Draw pairwise venn diagram
                 B = B,
                 C = C,
                 "A&B" = 582,
                 "A&C" = 145,
                 "B&C" = 153,
                 "A&B&C" = 227)))

plot(venneuler(c(A = A,          # Draw pairwise venn diagram
                 B = B,
                 C = C,
                 "A&B" = 214,
                 "A&C" = 12,
                 "B&C" = 20,
                 "A&B&C" = 22)))
###################################################
###################################################
###################################################
#plot a Volcano plot for DE phosphopeptides at difference timepoints

#make one df with all values in it
tmp.T30=rownames_to_column(table.T30)
tmp.T30=tmp.T30[,c(1,2,5,6)]
names(tmp.T30)[1]="PhosName"
names(tmp.T30)[2]="logFC.T30"
names(tmp.T30)[3]="P.Value.T30"
names(tmp.T30)[4]="adj.P.Val.T30"



#add better labels than PhosName
tmp_split <- tmp.T30 %>%
  separate(PhosName, into = c("Protein", "Gene", "Position", "AA.Centered"), sep = "~")

# Create the Gene and Position columns
tmp <- tmp_split %>%
  mutate(Position = paste(Gene, "[", Position, "]", sep = ""),
         Gene = gsub("^.*~", "", Gene),
         Rest = NULL)

#(2) add labeling variables

#highlight SZ risk genes
tmp$Gene.label=ifelse(tmp$Protein == "F1M0Z1"|
                        tmp$Protein == "O54898"|
                        tmp$Protein == "Q00959"|
                        tmp$Protein == "P22199"|
                        tmp$Protein == "P15865"|
                        tmp$Protein == "O70196"|
                        tmp$Protein == "O88382"|
                        tmp$Protein == "Q5YLM1"|
                        tmp$Protein == "Q02563"|
                        tmp$Protein == "P81795","YES",
                      "NO")

tmp$Group4=ifelse(tmp$adj.P.Val.T30 < 0.05 & tmp$Gene.label == "YES", "significant.label",
                  ifelse(tmp$adj.P.Val.T30 > 0.05 & tmp$Gene.label == "YES", "nonsignificant.label",
                         ifelse(tmp$adj.P.Val.T30 < 0.05 & tmp$Gene.label == "NO", "significant",
                                "nonsignificant")))



#######
#(3) volcano plot
#library(viridis)

cols=c("significant.label"="#0ad5ea",
       "nonsignificant.label"="#0ad5ea",
       "significant" = "#9D4EDD",
       "nonsignificant"="grey")

tmp$Group4=factor(tmp$Group4,levels=c("nonsignificant","significant","nonsignificant.label","significant.label"))
tmp=arrange(tmp,Group4)

ggplot(tmp, 
       aes(x=logFC.T30,
           y=-log10(P.Value.T30)),
       fill=Group4)+
  geom_point(shape=21,
             size=4,
             color="black",aes(fill=Group4))+
  scale_fill_manual(values=cols)+
  
  xlim(-1.75,1.75)+
  ylim(0,14)+
  theme_classic()+
  
  theme(legend.position = "none")+
  geom_text_repel(aes(label=ifelse(Group4 == "significant.label", as.character(c(Position)),"")),
                  max.overlaps=6000,
                  min.segment.length = 1
  )



######################################################## 
######################################################## 
######################################################## 
######################################################## 
######################################################## 
#####
###
##
#
#save(filtered.df2,
#     meta.pep,
 #    meta.data,
  #   limma.results,
   #  residuals.test,
   #  table.T30,
    # table.T300,
    # table.T900,
    # ppe,
    # file=paste0("PhosR/WorkData/PhosR_begin_final.RData"))





