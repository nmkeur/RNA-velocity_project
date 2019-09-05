#Set seed and load libraries
set.seed(60)
library(monocle)
library(VGAM)
library(Scribe)
library(plyr)
library(inflection)
library(igraph)
library(ggplot2)
library(GGally)
library(ggraph)
library(pracma)
library(lpSolveAPI)
library(cowplot)

#load CDS and Scribe results.
load('~/Desktop/CDS_predicted_velocity_10_genes.RData')

#extract and plot Raw RDI scores
ucrdi_list_velocity <- urdi_list_velocity$max_rdi_value
a <- pheatmap::pheatmap(ucrdi_list_velocity, cluster_rows = F, cluster_cols = F, annotation_names_col = T, border_color = NA,fontsize = 12)

# Filtering of raw RDI scores
ucrdi_list_velocity[ucrdi_list_velocity <= 0.] <- 0

#Filtering using dn sparisfier and CLR Z-scores.
ucrdi_list2 <- clr(ucrdi_list_velocity)

#Reload the raw scores and filter using z scores
load('~/Desktop/CDS_predicted_monocle_10_genes.RData')
ucrdi_list_velocity <- urdi_list_velocity$max_rdi_value
ucrdi_list_velocity[ucrdi_list2 <= 0] <-0

# Plot sparsified heatmap
b<- pheatmap::pheatmap(ucrdi_list_velocity, cluster_rows = F, cluster_cols = F, annotation_names_col = T, border_color = NA,fontsize = 12)

png("~/RNA-velocity_project/figure_25/figure/figure_25.png", width = 1500, height = 1800, units='px', res = 200)
plot_grid(a[[4]],b[[4]], labels=c("A","B"), nrow = 2)
dev.off()
