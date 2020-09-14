#!/usr/bin/Rscript

library("ggplot2")
library("readr")


##Gene lengths as violin plot

len = read.csv("contig_lengths.csv", header = TRUE)
len_df = as.data.frame(len)
df = stack(len_df)
pdf("contig_lengths_violin.pdf")

ggplot(df, aes(x = ind, y = values, fill = ind))+
  geom_violin(scale = "width") +
  labs(x = "CDS", y = "Sequence length") +
  guides(fill = FALSE) +
  theme_minimal()
  
dev.off()


##CDS count per gene

df$ftr = 1
df$ftr[is.na(df$values)] = 0
pdf("cds_count_separate.pdf")

ggplot(df, aes(x = reorder(ind, ftr),y = ftr)) +
  geom_col(color = "#00bfc4", fill = "#00bfc4") +
  labs(x = "CDS", y = "") +
  theme_minimal()

dev.off()


##Sequence lenght as histogram and features as bar plots

ftr = read.csv("features.csv", header = TRUE)
pdf("contig_lengths_hist.pdf")

ggplot(ftr, aes(ftr$length))+
  geom_histogram(binwidth = 130, col = "black", fill = "white") +
  labs(x = "Sequence length", y = "Frequency") +
  scale_y_continuous(expand = c(0,0),limits = c(0,20)) +
  scale_x_continuous(expand = c(0,0),limits = c(0,20000)) +
  geom_vline(xintercept = c(14000,17000), color = "orange") +
  theme_minimal()

dev.off()

pdf("gene_count.pdf")

ggplot(ftr, aes(ftr$gene, fill = factor(ifelse(gene == 13, "Highlightes", "Normal"))))+
  geom_bar() +
  labs(x = "Number of genes", y = "Count") +
  guides(fill = FALSE) +
  theme_minimal()

dev.off()

pdf("cds_count.pdf")

ggplot(ftr, aes(ftr$cds,fill = factor(ifelse(cds == 13, "Highlightes", "Normal"))))+
  geom_bar() +
  labs(x = "Number of CDSs", y = "Count") +
  guides(fill = FALSE) +
  theme_minimal()

dev.off()

pdf("trna_count.pdf")

ggplot(ftr, aes(ftr$trna, fill = factor(ifelse(trna == 22, "Highlightes", "Normal"))))+
  geom_bar() +
  labs(x = "Number of tRNAs", y = "Count") +
  guides(fill = FALSE) +
  theme_minimal()

dev.off()

pdf("rrna_count.pdf")

ggplot(ftr, aes(ftr$rrna, fill = factor(ifelse(rrna == 2, "Highlightes", "Normal"))))+
  geom_bar() +
  labs(x = "Number of rRNAs", y = "Count") +
  guides(fill = FALSE) +
  theme_minimal()
dev.off()


##Family level distribution as bar plot

family = read.csv("family_count.csv", header = TRUE)
pdf("family_count.pdf")
ggplot(family, aes(x = Family, y = Count))+
  geom_bar(stat = "identity", color = "#00bfc4", fill = "#00bfc4") +
  labs(x = "Family", y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()


##Country distribution as bar plot

country = read.csv("country_count.csv", header=TRUE)       
pdf("country_count.pdf")

ggplot(country, aes(x = Country, y = Count))+
  geom_bar(stat = "identity", color = "#00bfc4", fill = "#00bfc4") +
  labs(x = "Country", y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

dev.off()

# BUGS
# BUG1: trna and rrna plots for zero columns generated fucked histograms.
