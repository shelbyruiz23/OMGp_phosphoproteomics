#Grubisha omgp data analysis
#PhosR data analysis
#################
#08/15/2023
#################

rm(list=ls()); gc()
options(stringsAsFactors = F)

#setwd
setwd("C:/Users/Shelby_macdonaldlab/OneDrive - University of Pittsburgh/MacDonaldLab_Data/Collaborations/Grubisha_lab/Phospho_OMGp/")

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
#load files
load(paste0("PhosR/WorkData/PhosR_begin_final.RData"))

#read in normalized df
#average and output for RNfuzzyApp clustering
df=read.csv("PhosR/Results/files/residuals.csv")


#prepare for clustering analysis
#(1) filter all proteins for adj P value < 0.1 or top 25%
all=table.T30
all=rownames_to_column(all)

num_rows=round(nrow(all)*0.25) #2001 (round up)

top25=all %>%
  slice_min(order_by = adj.P.Val, n=num_rows)
all_s=top25

#(1.2) save selected phosphosites for extracting residuals of those phosphosites
list=all_s$rowname

#(1.3) extract rows from residual dataframe
df2=as.data.frame(df[df$X %in% list, ])

#average rows
# Identify the columns for each group
group_columns <- list(
  X1 = grep("^X1_", colnames(df2)),
  X30 = grep("^X30_", colnames(df2)),
  X300 = grep("^X300_", colnames(df2)),
  X900 = grep("^X900_", colnames(df2))
)

# Calculate row-wise averages for each group
df2$X1_avg <- rowMeans(df2[, group_columns$X1])
df2$X30_avg <- rowMeans(df2[, group_columns$X30])
df2$X300_avg <- rowMeans(df2[, group_columns$X300])
df2$X900_avg <- rowMeans(df2[, group_columns$X900])

df3=df2[,c(26:29)]
names(df3)[1]="X1"
names(df3)[2]="X30"
names(df3)[3]="X300"
names(df3)[4]="X900"

rownames(df3)=df2$X


#clustering with MFuzz
library(Mfuzz)


df3=as.matrix(df3)
eset=new("ExpressionSet", exprs=df3)

df4=Mfuzz::standardise(eset)

m1=mestimate(df4) #2.56

cl <- mfuzz(df4,c=6,m=m1)



mfuzz.plot2(df4,
            cl=cl,
            mfrow=c(2,3),
            min.mem = 0.3,
            
            time.points = c(0,100,370,970),
            time.labels = c("0","30 sec", "5 min", "15 min"),
                            xlab="",
            ylab="Abundance Change"
                            )
          


View(cl$size) #305, 336, 287, 311, 331, 431
cl.cluster=as.data.frame(cl$cluster)
cl.membership=as.data.frame(cl$membership)
cluster.info=cbind(cl.cluster,cl.membership)
names(cluster.info)[1]="Cluster"

#write.csv(cluster.info,"PhosR/Results/files/cluster.membership.csv")

mfuzzColorBar2.0(horizontal = FALSE)





###edit MFuzz color bar function to work for non "fancy" condition
mfuzzColorBar2.0= function (col, horizontal = FALSE, ...) 
{
  require(marray) || stop("Library marray is required")
  if (missing(col)) {
    col <- c("#FF8F00", "#FFA700", "#FFBF00", "#FFD700", 
             "#FFEF00", "#F7FF00", "#DFFF00", "#C7FF00", "#AFFF00", 
             "#97FF00", "#80FF00", "#68FF00", "#50FF00", "#38FF00", 
             "#20FF00", "#08FF00", "#00FF10", "#00FF28", "#00FF40", 
             "#00FF58", "#00FF70", "#00FF87", "#00FF9F", "#00FFB7", 
             "#00FFCF", "#00FFE7", "#00FFFF", "#00E7FF", "#00CFFF", 
             "#00B7FF", "#009FFF", "#0087FF", "#0070FF", "#0058FF", 
             "#0040FF", "#0028FF", "#0010FF", "#0800FF", "#2000FF", 
             "#3800FF", "#5000FF", "#6800FF", "#8000FF", "#9700FF", 
             "#AF00FF", "#C700FF", "#DF00FF", "#F700FF", "#FF00EF", 
             "#FF00D7", "#FF00BF", "#FF00A7", "#FF008F", "#FF0078", 
             "#FF0060", "#FF0048", "#FF0030", "#FF0018")
  }

  par(mar = c(5, 2, 4, 3) + 0.1)
  maColorBar(seq(0, 1, 0.01), col = col, horizontal = FALSE, 
             k = 11, ...)
}




