### Genotyping and analyzing round goby d-loop tissues 
# This script accomplishes the following:
# PART 1: Detecting heteroplasmy and genotyping goby individuals
# PART 2: D-loop genetic analysis

setwd("/Users/kbja10/Github/eDNA_abundance_estimation")

# Clear work environment and load packages
rm(list = ls())
library(plyr)
library(dplyr)
library(tidyr)
library (ggplot2)
library(pheatmap)
library(RColorBrewer)

#########################################################################
####### Detecting heteroplasmy and genotyping goby individuals ##########
#########################################################################

# Load DADA2 haplotyp-by-sample table 
dloop_table <- read.csv("/Users/kbja10/Desktop/d_loop_FR.seqtab.mod.csv")
dloop_table <- dloop_table[grep("_t.", colnames(dloop_table))] # subset to tissue samples (96)
colnames(dloop_table) <- gsub(".*__", "", colnames(dloop_table)) # remove duplicate sample names
dloop_table <- dloop_table[,colSums(dloop_table)>50] # drop 5 samples
dloop_table[dloop_table<10] <- 0 # don't count haplotypes w/ counts < 10
dloop_table <- dloop_table[rowSums(dloop_table)>0,] # 146 haplotypes remaining
heatmap(t(as.matrix(dloop_table)))
pheatmap(t(as.matrix(dloop_table)), treeheight_row = 0, treeheight_col = 0)

