### This code plots the qPCR results from the round goby abundance experiment
### Last updated 5.5.2022 by Kara Andres
### Part I: COI qPCR results -- 250 mL, 2L, 100L
### Part II: msat qPCR results -- 2L, 100L
### Part III: COI vs. msat -- 2L
### Part IV: contributor estimation vs. msat qPCR results

# Clear work environment and load packages
rm(list = ls())
library(plyr)
library(dplyr)
library(tidyr)
library(forcats)
library (ggplot2)
library(RColorBrewer)
library(ggpmisc)

################################################
### Part I: COI qPCR results -- 250 mL, 100L, 2L
################################################
# load COI datasets
metadata <- read.csv("/Users/kbja10/Documents/Cornell/Research/Round goby/abundance_eDNA/goby_eDNA_sample_sheet.csv")
qPCR_1 <- read.csv("/Users/kbja10/Documents/Cornell/Research/Round goby/abundance_eDNA/qPCR/COI/round_goby_qPCR_4.11.22.csv")
qPCR_2 <- read.csv("/Users/kbja10/Documents/Cornell/Research/Round goby/abundance_eDNA/qPCR/COI/round_goby_qPCR_4.15.22_1.csv")
qPCR_3 <- read.csv("/Users/kbja10/Documents/Cornell/Research/Round goby/abundance_eDNA/qPCR/COI/round_goby_qPCR_4.15.22_2.csv")
qPCR_4 <- read.csv("/Users/kbja10/Documents/Cornell/Research/Round goby/abundance_eDNA/qPCR/COI/round_goby_qPCR_4.20.22.csv")
COI_qPCR_results <- rbind(qPCR_1, qPCR_2, qPCR_3, qPCR_4)
COI_qPCR_results <- merge(COI_qPCR_results, metadata, by.x="Sample", by.y="Sample_ID") # add metadata

COI_qPCR_results <- COI_qPCR_results %>%
  filter(!grepl("STD",Sample)) %>%
  filter(!grepl("BL",Sample)) %>%
  mutate(Quantity = replace_na(Quantity, 0)) %>% # undetected == 0 copies
  mutate(Volume = case_when( # Add column for sample volume
    endsWith(Sample, "BL") ~ "Blank", 
    endsWith(Sample, ".100L") ~ "100L",
    endsWith(Sample, ".2L") ~ "2L",
    TRUE ~ "250mL")) %>%
  mutate(Quantity = case_when(
    Volume == "100L" ~ Quantity*4, # 100L samples eluted in 4x amount of buffer 
    TRUE ~ Quantity)) %>%
  mutate(Site = gsub("\\..*","", Sample)) %>%
  # mutate(Site = gsub("_B","", Site)) %>% # HAL_B & HHB_B
  arrange(Latitude)

# Plot copy number per sample
ggplot(COI_qPCR_results, aes(x=Sample, y=Quantity, fill=Volume)) +
  geom_boxplot(alpha=0.7) +
  scale_fill_brewer(palette="Set1") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 250 mL
p <- ggplot(COI_qPCR_results[COI_qPCR_results$Volume=="250mL",], aes(x=Sample, y=Quantity, fill=Site)) +
  geom_boxplot(alpha=0.7) +
  scale_fill_manual(values=colorRampPalette(brewer.pal(9, "Set1"))(14)) +
  ylab("DNA Quantity (COI copies/uL)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p
# ggsave(filename = paste("/Users/kbja10/Documents/Cornell/Research/Round goby/abundance_eDNA/Figures/qPCR_COI_250mL.pdf"), plot=p, dpi=300, width=8, height=8, units="in")

# 100 L
p <- ggplot(COI_qPCR_results[COI_qPCR_results$Volume=="100L"&COI_qPCR_results$Quantity>0,], aes(x=Sample, y=Quantity, fill=Site)) +
  geom_boxplot(alpha=0.7) +
  scale_fill_manual(values=colorRampPalette(brewer.pal(9, "Set1"))(14)) +
  ylab("DNA Quantity (COI copies/uL)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p
# ggsave(filename = paste("/Users/kbja10/Documents/Cornell/Research/Round goby/abundance_eDNA/Figures/qPCR_COI_100L.pdf"), plot=p, dpi=300, width=8, height=8, units="in")

# 2L
qPCR_5 <- read.csv("/Users/kbja10/Documents/Cornell/Research/Round goby/abundance_eDNA/qPCR/COI/round_goby_COI_qPCR_2L_4.23.21_1.csv")
qPCR_6 <- read.csv("/Users/kbja10/Documents/Cornell/Research/Round goby/abundance_eDNA/qPCR/COI/round_goby_COI_qPCR_2L_4.23.21_2.csv")
COI_qPCR_results_2L <- rbind(qPCR_5, qPCR_6)

COI_qPCR_results_2L <- COI_qPCR_results_2L %>%
  filter(!grepl("STD",Sample)) %>%
  filter(!grepl("BL",Sample)) %>%
  mutate(Quantity = replace_na(Quantity, 0)) %>% # undetected == 0 copies
  mutate(Site = gsub("\\..*","", Sample)) %>%
  mutate(Site = gsub("_B","", Site)) %>%
  arrange(Sample)

p <- ggplot(COI_qPCR_results_2L, aes(x=Sample, y=Quantity, fill=Site)) +
  geom_boxplot(alpha=0.7) +
  scale_fill_manual(values=colorRampPalette(brewer.pal(9, "Set1"))(24)) +
  ylab("DNA Quantity (COI copies/uL)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p
# ggsave(filename = paste("/Users/kbja10/Documents/Cornell/Research/Round goby/abundance_eDNA/Figures/qPCR_COI_2L.pdf"), plot=p, dpi=300, width=8, height=8, units="in")


###########################################
### Part II: msat qPCR results -- 2L, 100L
###########################################

# load datasets
qPCR_5 <- read.csv("/Users/kbja10/Documents/Cornell/Research/Round goby/abundance_eDNA/qPCR/microsatellite/round_goby_msat_qPCR_4.24.22.csv")
qPCR_6 <- read.csv("/Users/kbja10/Documents/Cornell/Research/Round goby/abundance_eDNA/qPCR/microsatellite/round_goby_msat_qPCR_5.1.22.csv")
qPCR_7 <- read.csv("/Users/kbja10/Documents/Cornell/Research/Round goby/abundance_eDNA/qPCR/microsatellite/round_goby_msat_qPCR_5.2.22.csv")
msat_qPCR_results <- rbind(qPCR_5, qPCR_6, qPCR_7)

msat_qPCR_results <- msat_qPCR_results %>%
  filter(!grepl("STD",Sample)) %>%
  filter(!grepl("BL",Sample)) %>%
  mutate(Quantity = replace_na(Quantity, 0)) %>% # undetected == 0 copies
  mutate(Volume = case_when( # Add column for sample volume
    endsWith(Sample, "BL") ~ "Blank", 
    endsWith(Sample, "100L") ~ "100L",
    endsWith(Sample, "2L") ~ "2L",
    TRUE ~ "250mL")) %>%
  mutate(Quantity = case_when(
    Volume == "100L" ~ Quantity*4, # 100L samples eluted in 4x amount of buffer 
    TRUE ~ Quantity)) %>%
  mutate(Site = gsub("\\..*","", Sample)) %>%
  mutate(Site = gsub("_B","", Site)) %>%
  arrange(Sample)

# Plot copy number per sample
p <- ggplot(msat_qPCR_results[msat_qPCR_results$Volume=="2L",], aes(x=Sample, y=Quantity, fill=Site)) +
  geom_boxplot(alpha=0.7) +
  scale_fill_manual(values=colorRampPalette(brewer.pal(9, "Set1"))(24)) +
  ylab("DNA Quantity (microsatellite copies/uL)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p
# ggsave(filename = paste("/Users/kbja10/Documents/Cornell/Research/Round goby/abundance_eDNA/Figures/qPCR_msat_2L.pdf"), plot=p, dpi=300, width=8, height=8, units="in")

p <- ggplot(msat_qPCR_results[msat_qPCR_results$Volume=="100L"&msat_qPCR_results$Quantity>0,], aes(x=Sample, y=Quantity, fill=Site)) +
  geom_boxplot(alpha=0.7) +
  scale_fill_manual(values=colorRampPalette(brewer.pal(9, "Set1"))(14)) +
  ylab("DNA Quantity (microsatellite copies/uL)") + ylim(0,20) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p
# ggsave(filename = paste("/Users/kbja10/Documents/Cornell/Research/Round goby/abundance_eDNA/Figures/qPCR_msat_100L.pdf"), plot=p, dpi=300, width=8, height=8, units="in")


########################################
### Part III: COI vs. msat results -- 2L
########################################

msat_qPCR_results_2L_means <- as.data.frame(msat_qPCR_results[msat_qPCR_results$Volume=="2L",] %>%
                                     group_by(Sample) %>%
                                     summarise(Quantity_msat = mean(Quantity)) %>%
                                     mutate(Site = gsub("\\..*","", Sample)) %>%
                                     arrange(Sample))
COI_qPCR_results_2L_means <- as.data.frame(COI_qPCR_results_2L %>%
                                    group_by(Sample) %>%
                                    summarise(Quantity_COI = mean(Quantity)) %>%
                                    mutate(Site = gsub("\\..*","", Sample)) %>%
                                    arrange(Sample))

qPCR_results_2L <- merge(qPCR_results_msat_means, COI_qPCR_results_2L_means, by="Sample")

p <- ggplot(qPCR_results_2L, aes(x=Quantity_msat, y=Quantity_COI, color=Site.x)) +
  geom_point(alpha=0.7, size=3) +
  scale_color_manual(values=colorRampPalette(brewer.pal(9, "Set1"))(24)) +
  ylab("Mitochondrial DNA concentration (copies/uL)") + 
  xlab("Nuclear DNA concentration (copies/uL)") + 
  theme_bw()
p
# ggsave(filename = paste("/Users/kbja10/Documents/Cornell/Research/Round goby/abundance_eDNA/Figures/msat_vs_COI_qPCR_2L.pdf"), plot=p, dpi=300, width=8, height=8, units="in")

qPCR_results_2L$Ratio <- qPCR_results_2L$Quantity_COI/qPCR_results_2L$Quantity_msat
qPCR_results_2L[is.infinite(qPCR_results_2L$Ratio),]$Ratio <- NA

p <- ggplot(qPCR_results_2L, aes(x=Site.x, y=Ratio, fill=Site.x)) +
  geom_boxplot(alpha=0.7) +
  scale_fill_manual(values=colorRampPalette(brewer.pal(9, "Set1"))(24)) +
  ylab("Mitochondrial:nuclear eDNA") + ylim(0,max(qPCR_results_2L$Ratio)) +
  xlab("Site") + 
  geom_hline(alpha=0.7, yintercept=mean(qPCR_results_2L$Ratio, na.rm=TRUE), size=1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "none")
p
# ggsave(filename = paste("/Users/kbja10/Documents/Cornell/Research/Round goby/abundance_eDNA/Figures/msat_COI_ratio_qPCR_2L.pdf"), plot=p, dpi=300, width=8, height=8, units="in")

################################################
### Part I: COI qPCR vs. AUV abundance
################################################

AUV_abundance <- read.csv("/Users/kbja10/Documents/Cornell/Research/Round goby/abundance_eDNA/Datasets/AUV_NumDens_062022_final.csv")
AUV_abundance <- read.csv("/Users/kbja10/Documents/Cornell/Research/Round goby/abundance_eDNA/Datasets/AUV_NumDens_transect.csv")
COI_qPCR_results_250mL <- COI_qPCR_results[COI_qPCR_results$Volume=="250mL",]
COI_qPCR_results_250mL_mean <- 
  as.data.frame(COI_qPCR_results_250mL %>%
  group_by(Site) %>%
  summarise(COI_mean=mean(Quantity, na.rm=TRUE)))
qPCR_adundance <- merge(COI_qPCR_results_250mL_mean, AUV_abundance, by.x="Site", by.y="site_ID")  

ggplot(qPCR_adundance, aes(x=COI_mean, y=transect.density.fish.m2.)) +
  stat_poly_line() +
  stat_poly_eq() +
  geom_point(size=3) +
  ylab("AUV abundance (density/m^2)") + xlab("DNA Quantity (COI copies/uL)") +
  theme_bw()

