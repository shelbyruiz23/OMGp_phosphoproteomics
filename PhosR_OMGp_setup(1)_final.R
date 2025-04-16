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


###################################################
#load and format files
###################################################
###################################################

#load and clean raw data file
raw=read.csv("PhosR/Input/Grubisha_omgp_phospho_raw.csv",header=T) #13,124

#remove any peptides with no Quan values
raw=raw %>%
  rowwise() %>%
  filter(!all(is.na(c_across(Abundances..Grouped...F1..126:Abundances..Grouped...F2..132N)))) # 13,009

#clean up raw data file
#keep only relevant columns, replace spaces with periods, make all sequences unique, reorder columns
raw.df=raw[,c(3,4,5,6,10,11,12,13,16:39)]
raw.df$Modifications.in.Master.Proteins=gsub(" ",".",raw.df$Modifications.in.Master.Proteins)

#remove peptides that are not unique to a single protein (isoforms considered to be the same protein)
#function to check if value prior to ; is equivalent (unique to the same protein)
is_equivalent <- function(x, y) {
  x_key <- sub("-.*", "", x)
  y_key <- sub("-.*", "", y)
  x_key == y_key
}

#function to extract the part before the hyphen '-'
extract_before_hyphen <- function(x) {
  sub("-.*", "", x)
}

#function to extract the part before the underscore '_'
extract_before_underscore <- function(x) {
  sub("_.*", "", x)
}

#remove rows where peptides map to multiple proteins
raw.df <- raw.df %>% 
  filter(
    !grepl(";", Master.Protein.Accessions) | sapply(Master.Protein.Accessions, function(x) any(sapply(strsplit(x, "; ")[[1]], function(y) is_equivalent(x, y))))
  ) #13009 -> 12790

#create a new protein column that only keeps the beginning uniprot id, regardless of isoform - group
raw.df$Protein <- sapply(raw.df$Master.Protein.Accessions, extract_before_hyphen)


#isolate peptide segment length,start and end into new columns
raw.df$Peptide.Segment <- gsub(".*?(\\[[0-9-]+\\]).*", "\\1", raw.df$Positions.in.Master.Proteins)
raw.df$Peptide.Start <- as.numeric(str_match(raw.df$Peptide.Segment, "\\[([0-9]+)")[, 2])
raw.df$Peptide.End <- as.numeric(str_match(raw.df$Peptide.Segment, "-([0-9]+)")[, 2])

#extract phosphoinformation into new column
raw.df$Mods <- gsub(".*(\\d+xPhospho \\[.*?\\]).*", "\\1", raw.df$Modifications)

#extract phospho information into meaningful columns
#make column value for how many phospho modifications the observed peptide has
raw.df$Phospho <- gsub("\\s*\\[.*", "", raw.df$Mods)

# Create empty columns
raw.df$P1x <- NA
raw.df$P2x <- NA
raw.df$P3x <- NA


# Extract values for P1x, P2x, and P3x columns
for (i in seq_along(raw.df$Mods)) {
  matches <- gregexpr("\\[(.*?)\\]", raw.df$Mods[i])[[1]]
  if (matches[[1]] != -1) {
    positions <- unlist(strsplit(regmatches(raw.df$Mods[i], matches)[[1]], "; "))
    raw.df$P1x[i] <- gsub("[\\[\\]]", "", positions[1])
    if (length(positions) > 1) {
      raw.df$P2x[i] <- gsub("[\\[\\]]", "", positions[2])
      if (length(positions) > 2) {
        raw.df$P3x[i] <- gsub("[\\[\\]]", "", positions[3])
        if (grepl("/", raw.df$P3x[i])) {
          raw.df$P1x[i] <- gsub("^.*?/", "", raw.df$P3x[i])
          raw.df$P2x[i] <- gsub("^.*?;", "", raw.df$P3x[i])
          raw.df$P3x[i] <- gsub(";.*$", "", raw.df$P3x[i])
        }
      }
    }
  }
  
  # Handle special cases for 2xPhospho and 3xPhospho
  if (raw.df$Phospho[i] == "2xPhospho" && is.na(raw.df$P2x[i])) {
    raw.df$P2x[i] <- raw.df$P1x[i]
  }
  if (raw.df$Phospho[i] == "3xPhospho") {
    if (is.na(raw.df$P2x[i])) {
      raw.df$P2x[i] <- raw.df$P1x[i]
    }
    if (length(gregexpr(";", raw.df$Mods[i])[[1]]) == 1) {
      positions <- unlist(strsplit(regmatches(raw.df$Mods[i], matches)[[1]], "; "))
      raw.df$P2x[i] <- gsub("[\\[\\]]", "", positions[2])
      raw.df$P3x[i] <- gsub("[\\[\\]]", "", positions[2])
    }
    if (!grepl(";", raw.df$Mods[i]) && is.na(raw.df$P2x[i])) {
      raw.df$P2x[i] <- raw.df$P1x[i]
      raw.df$P3x[i] <- raw.df$P1x[i]
    }
  }
}

#now remove brackets from values in P1x, P2x, and P3x
raw.df$P1x=gsub("\\[|\\]", "", raw.df$P1x)
raw.df$P2x=gsub("\\[|\\]", "", raw.df$P2x)
raw.df$P3x=gsub("\\[|\\]", "", raw.df$P3x)


# Separate characters and numeric values
raw.df$P1x_character <- gsub("\\d+", "", raw.df$P1x)
raw.df$P1x_numeric <- as.numeric(gsub("\\D+", "", raw.df$P1x))

raw.df$P2x_character <- gsub("\\d+", "", raw.df$P2x)
raw.df$P2x_numeric <- as.numeric(gsub("\\D+", "", raw.df$P2x))

raw.df$P3x_character <- gsub("\\d+", "", raw.df$P3x)
raw.df$P3x_numeric <- as.numeric(gsub("\\D+", "", raw.df$P3x))

#add correct site information to phoshosite columns
raw.df$P1x_site=(raw.df$Peptide.Start + (raw.df$P1x_numeric - 1))
raw.df$P2x_site=(raw.df$Peptide.Start + (raw.df$P2x_numeric - 1))
raw.df$P3x_site=(raw.df$Peptide.Start + (raw.df$P3x_numeric - 1))

#paste back residue with site information (relative to master protein)
raw.df$P1x_Protsite=ifelse(!is.na(raw.df$P1x_character) & !is.na(raw.df$P1x_site),
                           paste0("[", paste(raw.df$P1x_character, raw.df$P1x_site, sep = ""), "]"),
                           ifelse(!is.na(raw.df$P1x_character), paste0("[", paste(raw.df$P1x_character), "]")))

raw.df$P2x_Protsite=ifelse(!is.na(raw.df$P2x_character) & !is.na(raw.df$P2x_site),
                           paste0("[", paste(raw.df$P2x_character, raw.df$P2x_site, sep = ""), "]"),
                           ifelse(!is.na(raw.df$P2x_character), paste0("[", paste(raw.df$P2x_character), "]"),NA))

raw.df$P3x_Protsite=ifelse(!is.na(raw.df$P3x_character) & !is.na(raw.df$P3x_site),
                           paste0("[", paste(raw.df$P3x_character, raw.df$P3x_site, sep = ""), "]"),
                           ifelse(!is.na(raw.df$P3x_character), paste0("[", paste(raw.df$P3x_character), "]"),NA))



###################################################
###################################################
###################################################
###################################################
#make all phosphopeptides unique by adding their aa sequence and modifications
raw.df$Sequence=paste(raw.df$Annotated.Sequence,raw.df$Protein,raw.df$Phospho,raw.df$P1x_Protsite,raw.df$P2x_Protsite,raw.df$P3x_Protsite, sep="_")
#remove NAs from ending
raw.df$Sequence <- gsub("_NA(_NA)?$", "", raw.df$Sequence)

#check that all sequences are unique, (numrows(t)=numrows(raw.df))
t=unique(raw.df$Sequence) #12790 != 12384
#since not all peptides are unique, choose one representative peptide
#duplicate peptides with the same sequence and same phospho sites may occur due to one peptide being observed with oxidation modification and other not
#choose peptides that have the least amount of NA values 
#in cases of a tie, highest peptide quan values rowsum
#one single peptide is chosen based on having the lowest Qvality.pep value

#isolate any nonduplicated sequences from the rest of the dataframe
raw.df2=raw.df[!duplicated(raw.df$Sequence) & !duplicated(raw.df$Sequence, fromLast = TRUE), ] #11,993

#isolate any duplicate sequences from the rest of the dataframe and choose one representative peptide
tmp=raw.df[duplicated(raw.df$Sequence) | duplicated(raw.df$Sequence, fromLast = TRUE), ] # 797


#now reduce duplicates of the tmp subrows to one representative peptide
#define plex values
plex1_columns=colnames(tmp[,c(9:20)])
plex2_columns=colnames(tmp[,c(21:32)])

# Your dataframe 
df <- tmp

# Step 1: Loop along each set of duplicated "Sequence" values
unique_sequences <- unique(df$Sequence)
final_rows <- vector("numeric", length(unique_sequences))

# Initialize an empty data frame to store the final results
final_df <- data.frame()

for (i in seq_along(unique_sequences)) {
  sequence <- unique_sequences[i]
  subset_df <- df[df$Sequence == sequence, ]
  
  # Step 2: Determine which row to keep based on NA values in Plex1 and Plex2
  na_counts <- rowSums(is.na(subset_df[, c(plex1_columns, plex2_columns)]))
  
  # Step 3: Keep rows with no NA values, if available
  no_na_rows <- subset_df[na_counts == 0, ]
  if (nrow(no_na_rows) > 0) {
    selected_row <- no_na_rows[which.min(no_na_rows$Qvality.PEP), ]
  } else {
    # Step 4: Keep the row with least NA values and largest total row sum
    least_na_rows <- subset_df[na_counts == min(na_counts), ]
    total_row_sums <- rowSums(least_na_rows[, c(plex1_columns, plex2_columns)], na.rm = TRUE)
    top_two_indices <- order(-total_row_sums)[1:2]
    selected_row <- least_na_rows[top_two_indices[1], ]
    
  }
  
  final_df <- rbind(final_df, selected_row)
}

# Step 5: Keep unique "Sequence" values based on "Qvality.PEP"
final_df <- final_df %>%
  group_by(Sequence) %>%
  filter(Qvality.PEP == min(Qvality.PEP) | Sequence %in% unique_sequences) %>%
  ungroup() # 391


#now, add rows from tmp back to original dataframe
raw.df=rbind(raw.df2,final_df)

#check that all peptides are unique
t=unique(raw.df$Sequence) #12384 = 12384

#remove unneeded dataframes
rm(tmp)
rm(raw.df2)
rm(final_df)
rm(df)
rm(least_na_rows)
rm(no_na_rows)
rm(selected_row)
rm(subset_df)

###################################################
###################################################
###################################################
###################################################

#add gene names and protein information to raw data file
#load rat uniprot export (export dates may vary protein/gene mapping slightly)
rat.db=read_excel("PhosR/Input/uniprotkb_rat_export2023_07_31.xlsx")
names(rat.db)[1]="Protein"
names(rat.db)[9]="AA.Sequence"
rat.db=rat.db[,-c(2,6,8)]

raw.df=left_join(raw.df,rat.db,by="Protein")
names(raw.df)[55]='Gene'
names(raw.df)[56]='Protein.names'
names(raw.df)[57]='Gene.names'
raw.df$Gene <- sapply(raw.df$Gene, extract_before_underscore) #12,384

#####
###
##
#

#remove any row that does not have a specific site
#start with phosphorylated peptides with no known specific sites
filtered.df <- subset(raw.df, !(Phospho == "1xPhospho" & is.na(P1x_numeric))) #9798
filtered.df <- subset(filtered.df, !(Phospho == "2xPhospho" & is.na(P1x_numeric & is.na(P2x_numeric)))) #9532
filtered.df <- subset(filtered.df, !(Phospho == "3xPhospho" & is.na(P1x_numeric & is.na(P2x_numeric) & is.na(P3x_numeric)))) #9504


#can do this multiple ways
#(1)can just take the first specific site on 2x and 3xphospho peptides

# Create new columns and choose representative values
filtered.df$Rep_numeric <- apply(filtered.df[, c("P1x_site", "P2x_site", "P3x_site")], 1, function(x) x[which(!is.na(x))[1]])
#this works nicely bc there are no multi-phosphorylated sites that dont have a specific known position first
filtered.df$Rep_character <- apply(filtered.df[, c("P1x_character", "P2x_character", "P3x_character")], 1, function(x) x[which(!is.na(x))[1]])

#remove random phosphopeptides
#remove phosphopeptides whose target start site is larger than total protein length
filtered.df$Length=as.numeric(filtered.df$Length)
filtered.df=filtered.df[filtered.df$Rep_numeric <= filtered.df$Length, ]
#remove phosphopeptides that are not mapped to a protein
filtered.df=filtered.df[!is.na(filtered.df$Protein), ] #9,487

#####
###
##
#

#####
###
##
#
#make phosphosite of interest centered
filtered.df <- filtered.df %>%
  mutate(
    Seq_before_start = substr(AA.Sequence, Rep_numeric - 7, Rep_numeric - 1),
    Seq_before_start = str_pad(Seq_before_start, width = 7, pad = "_", side = "left"),
    Seq_after_start = substr(AA.Sequence, Rep_numeric + 1, Rep_numeric + 7),
    Seq_after_start = str_pad(Seq_after_start, width = 7, pad = "_", side = "right")
  )

filtered.df$AA.centered=paste(filtered.df$Seq_before_start,filtered.df$Rep_character,filtered.df$Seq_after_start,sep="")

#remove columns not needed
filtered.df2=filtered.df[,c(-62,-63)] #9,487

#####
###
##
#

###################################################
###################################################
###################################################
###################################################

#make column for desired PhosR naming scheme
filtered.df2$Rep=paste(filtered.df2$Rep_character,filtered.df2$Rep_numeric,sep="")
filtered.df2$PhosName=paste(filtered.df2$Protein,filtered.df2$Gene,filtered.df2$Rep,filtered.df2$AA.centered,sep="~")

#since multiple peptides measured can overlap with the same site, check that each site has only one row value
t=unique(filtered.df2$PhosName) #9487 != 7757


#isolate any nonduplicated Phosnames from the rest of the dataframe
filtered.df3=filtered.df2[!duplicated(filtered.df2$PhosName) & !duplicated(filtered.df2$PhosName, fromLast = TRUE), ] #6440

#isolate any duplicate sequences from the rest of the dataframe and then additionally split by number of phosphorylations
tmp=filtered.df2[duplicated(filtered.df2$PhosName) | duplicated(filtered.df2$PhosName, fromLast = TRUE), ] #3047
tmp.1x=filter(tmp, Phospho == "1xPhospho") #2316
tmp.2x=filter(tmp, Phospho == "2xPhospho") #657
tmp.3x=filter(tmp, Phospho == "3xPhospho") #74

############################################
#start with simplest case (1xPhospho) and repeat like for nonunique sequences

#define plex values
plex1_columns=colnames(tmp.1x[,c(9:20)])
plex2_columns=colnames(tmp.1x[,c(21:32)])

# Your dataframe 
df <- tmp.1x

# Step 1: Loop along each set of duplicated "Sequence" values
unique_sequences <- unique(df$PhosName)
final_rows <- vector("numeric", length(unique_sequences))

# Initialize an empty data frame to store the final results
final_df <- data.frame()

for (i in seq_along(unique_sequences)) {
  sequence <- unique_sequences[i]
  subset_df <- df[df$PhosName == sequence, ]
  
  # Step 2: Determine which row to keep based on NA values in Plex1 and Plex2
  na_counts <- rowSums(is.na(subset_df[, c(plex1_columns, plex2_columns)]))
  
  # Step 3: Keep rows with no NA values, if available
  no_na_rows <- subset_df[na_counts == 0, ]
  if (nrow(no_na_rows) > 0) {
    selected_row <- no_na_rows[which.min(no_na_rows$Qvality.PEP), ]
  } else {
    # Step 4: Keep the row with least NA values and largest total row sum
    least_na_rows <- subset_df[na_counts == min(na_counts), ]
    total_row_sums <- rowSums(least_na_rows[, c(plex1_columns, plex2_columns)], na.rm = TRUE)
    top_two_indices <- order(-total_row_sums)[1:2]
    selected_row <- least_na_rows[top_two_indices[1], ]
    
 
  }
  
  final_df <- rbind(final_df, selected_row)
}

# Step 5: Keep unique "Sequence" values based on "Qvality.PEP"
final_df <- final_df %>%
  group_by(PhosName) %>%
  filter(Qvality.PEP == min(Qvality.PEP) | PhosName %in% unique_sequences) %>%
  ungroup()


#store final df as 1xPhospho rows to read to filter.df3 dataframe
tmp.1x=final_df #1173

############################################
t=unique(tmp.2x$PhosName) # 418


############################################
#now look at a more complex case (2xPhospho) and adjust
#separate out 2xPhospho sites that dont have duplicates with another 2xPhospho site
tmp.2x_unique=tmp.2x[!duplicated(tmp.2x$PhosName) & !duplicated(tmp.2x$PhosName, fromLast = TRUE), ] #242

#separate out 2xPhospho sites that do have duplicates with another 2xPhospho site
tmp.2x_n=tmp.2x[duplicated(tmp.2x$PhosName) | duplicated(tmp.2x$PhosName, fromLast = TRUE), ] #415

#remove any rows that dont have a known secondary site 
tmp.2x_n=tmp.2x_n[!is.na(tmp.2x_n$P2x_numeric), ] # 379

#add any new unique rows to unique rows of tmp.2x_unique
tmp.2x_unique2=tmp.2x_n[!duplicated(tmp.2x_n$PhosName) & !duplicated(tmp.2x_n$PhosName, fromLast=TRUE), ] # 20
tmp.2x_unique=rbind(tmp.2x_unique,tmp.2x_unique2) #263

#separate out 2xPhospho sites that do have duplicates with another 2xPhospho site (again)
tmp.2x_n2=tmp.2x_unique[duplicated(tmp.2x_unique$PhosName) | duplicated(tmp.2x_unique$PhosName, fromLast = TRUE), ] #0
rm(tmp.2x_n2)

#now isolate duplicated rows of tmp.2x_n
tmp.2x_n=tmp.2x_n[duplicated(tmp.2x_n$PhosName) | duplicated(tmp.2x_n$PhosName, fromLast = TRUE), ] #359


#########################################
#########################################
#make a unique phosName based on the second phospho site
# Your dataframe
df <- tmp.2x_n

# Step 1: Loop along each set of duplicated "PhosName" values
unique_phos_names <- unique(df$PhosName)
final_rows <- list()

for (phos_name in unique_phos_names) {
  subset_df <- df[df$PhosName == phos_name, ]
  
  # Step 2: For exactly 2 duplicates
  if (nrow(subset_df) == 2) {
    row1 <- subset_df[1, ]
    row2 <- subset_df[2, ]
    
    row1$Rep_numeric <- row1$P1x_site
    row1$Rep_character <- row1$P1x_character
    
    row2$Rep_numeric <- row2$P2x_site
    row2$Rep_character <- row2$P2x_character
    
    final_rows <- c(final_rows, list(row1, row2))
  }
  # Step 3: For more than 2 duplicates
  else if (nrow(subset_df) >= 3) {
    top_2_rows <- subset_df %>% 
      arrange(Qvality.PEP) %>% 
      slice_head(n = 2)
    
    row1 <- top_2_rows[1, ]
    row2 <- top_2_rows[2, ]
    
    row1$Rep_numeric <- row1$P1x_site
    row1$Rep_character <- row1$P1x_character
    
    row2$Rep_numeric <- row2$P2x_site
    row2$Rep_character <- row2$P2x_character
    
    final_rows <- c(final_rows, list(row1, row2))
  }
}

# Convert the list of rows to a dataframe
final_df <- bind_rows(final_rows)

#edit the AA.centered column
#make phosphosite of interest centered
final_df <- final_df %>%
  mutate(
    Seq_before_start = substr(AA.Sequence, Rep_numeric - 7, Rep_numeric - 1),
    Seq_before_start = str_pad(Seq_before_start, width = 7, pad = "_", side = "left"),
    Seq_after_start = substr(AA.Sequence, Rep_numeric + 1, Rep_numeric + 7),
    Seq_after_start = str_pad(Seq_after_start, width = 7, pad = "_", side = "right")
  )

final_df$AA.centered=paste(final_df$Seq_before_start,final_df$Rep_character,final_df$Seq_after_start,sep="")

#remove columns not needed
final_df=final_df[,c(-65,-66)]

#edit the Rep column
final_df$Rep=paste(final_df$Rep_character,final_df$Rep_numeric,sep="")
#edit the PhosName column
final_df$PhosName=paste(final_df$Protein,final_df$Gene,final_df$Rep,final_df$AA.centered,sep="~")

tmp.2x_n=final_df

#rowbind previously duplicated 2xPhospho values with unique 2xPhospho values
tmp.2x_unique=rbind(tmp.2x_unique,tmp.2x_n)

#check how many unique values now
t=unique(tmp.2x_unique$PhosName) #557/564

#isolate duplicated and nonduplicated values from tmp.2x_unique
tmp.2x_unique2=tmp.2x_unique[!duplicated(tmp.2x_unique$PhosName) & !duplicated(tmp.2x_unique$PhosName, fromLast = TRUE), ] #550
tmp.2x_n=tmp.2x_unique[duplicated(tmp.2x_unique$PhosName) | duplicated(tmp.2x_unique$PhosName, fromLast = TRUE), ] #14


############################################
#find the best of the remain sequences from tmp.2x_n to use

#define plex values
plex1_columns=colnames(tmp.2x_n[,c(9:20)])
plex2_columns=colnames(tmp.2x_n[,c(21:32)])

# Your dataframe 
df <- tmp.2x_n

# Step 1: Loop along each set of duplicated "Sequence" values
unique_sequences <- unique(df$PhosName)
final_rows <- vector("numeric", length(unique_sequences))

# Initialize an empty data frame to store the final results
final_df <- data.frame()

for (i in seq_along(unique_sequences)) {
  sequence <- unique_sequences[i]
  subset_df <- df[df$PhosName == sequence, ]
  
  # Step 2: Determine which row to keep based on NA values in Plex1 and Plex2
  na_counts <- rowSums(is.na(subset_df[, c(plex1_columns, plex2_columns)]))
  
  # Step 3: Keep rows with no NA values, if available
  no_na_rows <- subset_df[na_counts == 0, ]
  if (nrow(no_na_rows) > 0) {
    selected_row <- no_na_rows[which.min(no_na_rows$Qvality.PEP), ]
  } else {
    # Step 4: Keep the row with least NA values and largest total row sum
    least_na_rows <- subset_df[na_counts == min(na_counts), ]
    total_row_sums <- rowSums(least_na_rows[, c(plex1_columns, plex2_columns)], na.rm = TRUE)
    top_two_indices <- order(-total_row_sums)[1:2]
    selected_row <- least_na_rows[top_two_indices[1], ]

  }
  
  final_df <- rbind(final_df, selected_row)
}

# Step 5: Keep unique "Sequence" values based on "Qvality.PEP"
final_df <- final_df %>%
  group_by(PhosName) %>%
  filter(Qvality.PEP == min(Qvality.PEP) | PhosName %in% unique_sequences) %>%
  ungroup()


#store final df as 1xPhospho rows to readd to filter.df3 dataframe
tmp.2x_n=final_df #7

############################################
#rowbind previously duplicated 2xPhospho values with unique 2xPhospho values
tmp.2x_unique=rbind(tmp.2x_unique2,tmp.2x_n) #557

#########################################
#########################################
#now look at a more complex case (3xPhospho) and adjust
#separate out 3xPhospho sites that dont have duplicates with another 3xPhospho site
tmp.3x_unique=tmp.3x[!duplicated(tmp.3x$PhosName) & !duplicated(tmp.3x$PhosName, fromLast = TRUE), ] #44

#separate out 3xPhospho sites that do have duplicates with another 3xPhospho site
tmp.3x_n=tmp.3x[duplicated(tmp.3x$PhosName) | duplicated(tmp.3x$PhosName, fromLast = TRUE), ] #30

############################################

# Your dataframe
df <- tmp.3x_n

# Step 1: Loop along each set of duplicated "PhosName" values
unique_phos_names <- unique(df$PhosName)
final_rows <- list()

for (phos_name in unique_phos_names) {
  subset_df <- df[df$PhosName == phos_name, ]
  
  num_duplicates <- nrow(subset_df)
  
  if (num_duplicates <= 3) {
    # Step 2: For less than or equal to 3 duplicates
    sorted_df <- subset_df %>%
      arrange(is.na(P2x_site), is.na(P3x_site), Qvality.PEP) %>%
      slice_head(n = num_duplicates)
    
    for (i in 1:num_duplicates) {
      selected_row <- sorted_df[i, ]
      selected_row$Rep_numeric <- selected_row[[paste0("P", i, "x_site")]]
      selected_row$Rep_character <- selected_row[[paste0("P", i, "x_character")]]
      final_rows <- c(final_rows, list(selected_row))
    }
  } else {
    # Step 3: For more than or equal to 4 duplicates
    top_3_rows <- subset_df %>%
      arrange(Qvality.PEP) %>%
      slice_head(n = 3)
    
    sorted_df <- top_3_rows %>%
      arrange(is.na(P2x_site), is.na(P3x_site))
    
    for (i in 1:3) {
      selected_row <- sorted_df[i, ]
      selected_row$Rep_numeric <- selected_row[[paste0("P", i, "x_site")]]
      selected_row$Rep_character <- selected_row[[paste0("P", i, "x_character")]]
      final_rows <- c(final_rows, list(selected_row))
    }
  }
}

# Convert the list of rows to a dataframe
final_df <- bind_rows(final_rows)

#remove rows that have NA values in Rep_numeric column
final_df <- final_df[complete.cases(final_df$Rep_numeric), ]

#edit the AA.centered column
#make phosphosite of interest centered
final_df <- final_df %>%
  mutate(
    Seq_before_start = substr(AA.Sequence, Rep_numeric - 7, Rep_numeric - 1),
    Seq_before_start = str_pad(Seq_before_start, width = 7, pad = "_", side = "left"),
    Seq_after_start = substr(AA.Sequence, Rep_numeric + 1, Rep_numeric + 7),
    Seq_after_start = str_pad(Seq_after_start, width = 7, pad = "_", side = "right")
  )

final_df$AA.centered=paste(final_df$Seq_before_start,final_df$Rep_character,final_df$Seq_after_start,sep="")

#remove columns not needed
final_df=final_df[,c(-65,-66)]

#edit the Rep column
final_df$Rep=paste(final_df$Rep_character,final_df$Rep_numeric,sep="")
#edit the PhosName column
final_df$PhosName=paste(final_df$Protein,final_df$Gene,final_df$Rep,final_df$AA.centered,sep="~")

tmp.3x_n=final_df

#rowbind previously duplicated 3xPhospho values with unique 3xPhospho values
tmp.3x_unique=rbind(tmp.3x_unique,tmp.3x_n) #70


############################################
#bin back together filtered.df3,tmp.1x,tmp.2x_unique,tmp.3x_unique and check for new duplicates
filtered.df4=rbind(filtered.df3,tmp.1x,tmp.2x_unique,tmp.3x_unique) #8240
t=unique(filtered.df4$PhosName) #7849

#isolate any nonduplicated Phosnames from the rest of the dataframe
filtered.df5=filtered.df4[!duplicated(filtered.df4$PhosName) & !duplicated(filtered.df4$PhosName, fromLast = TRUE), ] #7464

#isolate any duplicate sequences from the rest of the dataframe and then additionally split by number of phosphorylations
tmp=filtered.df4[duplicated(filtered.df4$PhosName) | duplicated(filtered.df4$PhosName, fromLast = TRUE), ] #776



#########################################
#########################################
###
df=tmp

# Step 1: Loop along each set of duplicated "PhosName" values
unique_phos_names <- unique(df$PhosName)

# Initialize an empty dataframe to store the results
final_df <- data.frame()

for (phos_name in unique_phos_names) {
  subset_df <- df[df$PhosName == phos_name, ]
  
  # Step 2: Determine the unique values
  unique_values <- unique(na.omit(c(subset_df$P1x_site, subset_df$P2x_site, subset_df$P3x_site)))
  
  # Step 3: Arrange rows within the duplicated set
  sorted_df <- subset_df %>%
    arrange(
      desc(is.na(P2x_site) | is.na(P3x_site)),
      Qvality.PEP
    )
  
  # Initialize variables to keep track of retained values and rows
  retained_values <- c()
  retained_rows <- data.frame()
  
  # Step 4: Identify rows to retain
  for (idx in 1:nrow(sorted_df)) {
    current_row <- sorted_df[idx, ]
    row_values <- c(current_row$P1x_site, current_row$P2x_site, current_row$P3x_site)
    
    # Find the first unique numeric value in the row
    unique_col_value <- NULL
    unique_col_character <- NULL
    for (col_value in row_values) {
      if (!is.na(col_value) && !(col_value %in% retained_values)) {
        unique_col_value <- col_value
        retained_values <- c(retained_values, col_value)
        
        if (col_value == current_row$P1x_site) {
          unique_col_character <- current_row$P1x_character
        } else if (col_value == current_row$P2x_site) {
          unique_col_character <- current_row$P2x_character
        } else if (col_value == current_row$P3x_site) {
          unique_col_character <- current_row$P3x_character
        }
        
        break
      }
    }
    
    # If a unique value is found, retain the row
    if (!is.null(unique_col_value)) {
      current_row$Rep_numeric <- unique_col_value
      current_row$Rep_character <- unique_col_character
      retained_rows <- rbind(retained_rows, current_row)
    }
    
    if (length(retained_values) == length(unique_values)) {
      break
    }
  }
  
  # Add retained rows to the final dataframe
  final_df <- rbind(final_df, retained_rows)
}

#edit the AA.centered column
#make phosphosite of interest centered
final_df <- final_df %>%
  mutate(
    Seq_before_start = substr(AA.Sequence, Rep_numeric - 7, Rep_numeric - 1),
    Seq_before_start = str_pad(Seq_before_start, width = 7, pad = "_", side = "left"),
    Seq_after_start = substr(AA.Sequence, Rep_numeric + 1, Rep_numeric + 7),
    Seq_after_start = str_pad(Seq_after_start, width = 7, pad = "_", side = "right")
  )

final_df$AA.centered=paste(final_df$Seq_before_start,final_df$Rep_character,final_df$Seq_after_start,sep="")

#remove columns not needed
final_df=final_df[,c(-65,-66)]

#edit the Rep column
final_df$Rep=paste(final_df$Rep_character,final_df$Rep_numeric,sep="")
#edit the PhosName column
final_df$PhosName=paste(final_df$Protein,final_df$Gene,final_df$Rep,final_df$AA.centered,sep="~")

tmp=final_df

#rowbind previously duplicated tmp values with all unique values
filtered.df6=rbind(filtered.df5,tmp) #8133

#check that all PhosNames are unique
t=unique(filtered.df6$PhosName)#8088


##########repeat as before
#isolate any nonduplicated Phosnames from the rest of the dataframe
filtered.df7=filtered.df6[!duplicated(filtered.df6$PhosName) & !duplicated(filtered.df6$PhosName, fromLast = TRUE), ] #8045

#isolate any duplicate sequences from the rest of the dataframe and then additionally split by number of phosphorylations
tmp=filtered.df6[duplicated(filtered.df6$PhosName) | duplicated(filtered.df6$PhosName, fromLast = TRUE), ] #88


#########################################
#########################################
###

#at this point choose best phosphosite with lowest Qvality.PEP value
tmp <- tmp %>%
  group_by(PhosName) %>%
  slice_min(Qvality.PEP) %>%
  ungroup() #43

#rowbind previously duplicated tmp values with all unique values
filtered.df8=rbind(filtered.df7,tmp) #8088

#check that all PhosNames are unique
t=unique(filtered.df8$PhosName)#8088


#DONE with phosphosite reduction

###################################################
###################################################
###################################################
###################################################


################################
################################
################################

#make meta.pep file
meta.pep=filtered.df8[,c(64,54,1,2,6,7,8,33:53,55:63)]



#keep only values for phosphosite abundance measurments
#remove irrelevant peptide information
df=filtered.df8[,c(9:32)]
rownames(df)=filtered.df8$PhosName

#####
###
##
#


#read in metafile with sample information
meta.data=read.csv("PhosR/Input/Grubisha_omgp_phospho_meta.csv")

#reformat file name column in metadata to match column titles in df 
meta.data$File=gsub(" |:|,|\\(|\\)",".",meta.data$File)

#turn meta data into characters
#add T to Time column
meta.data$Time=sub("^","T",meta.data$Time)
meta.data$Plate=sub("^","P",meta.data$Plate)

#replace column headers with meaningful names
col_names=match(names(df),meta.data$File)
new_names=meta.data$Name[na.omit(col_names)]
names(df)=new_names

#reformat meta.data for limma
rownames(meta.data)=meta.data$Name
meta.data=meta.data[,c(2,3,4,6,7,8)]



#####
###
##
#
#save(df,
 #    meta.pep,
  #   meta.data,
   #  rat.db,
    # file=paste0("PhosR/WorkData/PhosR_setup_final.RData"))
