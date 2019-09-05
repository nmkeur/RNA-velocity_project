#Load required libraries
library(monocle)
library(cowplot)

#Load CDS file used in analysis
CDS_predicted <- readRDS('~/Desktop/CDS_predicted_fulldataset_INF.Rds')

#Change the variable names for figures
names(pData(CDS_predicted))[names(pData(CDS_predicted)) == 'batch_name'] <- 'Batch'
names(pData(CDS_predicted))[names(pData(CDS_predicted)) == 'stim'] <- 'Stimulation'
CDS_predicted$Batch <- as.character(CDS_predicted$Batch)
names(pData(CDS_predicted))[names(pData(CDS_predicted)) == 'sample'] <- 'Sample'

sel_genes <- c("GBP2","IRF1","IRF2","JAK1","STAT1")
CDS_predicted <- CDS_predicted[CDS_predicted@featureData@data$gene_short_name %in%
                                 sel_genes,]

a <- plot_genes_jitter(CDS_predicted,
                       color_by = 'Sample',
                       grouping = 'Sample',
                       ncol=2, 
                       cell_size = 0.2, 
                       min_expr = 0.9)+ 
  guides(colour = guide_legend(override.aes = list(size=3)))


png("~/RNA-velocity_project/figure_12/figure/figure_12.png", width = 2000, height = 1800, units='px', res = 200)
plot_grid(a)
dev.off()
