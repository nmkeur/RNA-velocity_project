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


a <- plot_genes_jitter(CDS_predicted[c('IRF1','IRF2','IRF3','STAT1')],
                       color_by = 'Stimulation',
                       grouping = 'Pseudotime',
                       ncol=2, 
                       cell_size = 0.2)+ 
  guides(colour = guide_legend(override.aes = list(size=3)))

b <- plot_genes_positive_cells(CDS_predicted[c('IRF1','IRF2','IRF3','STAT1')],
                       grouping = 'Stimulation',
                       ncol=2)+ 
  guides(colour = guide_legend(override.aes = list(size=3)))

c <- plot_genes_jitter(CDS_predicted[c('HLA-A','HLA-E','IFI6','IFIT1','IFIT3','XAF1')],
                       color_by = 'Stimulation',
                       grouping = 'Pseudotime',
                       ncol=2, 
                       cell_size = 0.2)+ 
  guides(colour = guide_legend(override.aes = list(size=3)))

d <- plot_genes_positive_cells(CDS_predicted[c('HLA-A','HLA-E','IFI6','IFIT1','IFIT3','XAF1')],
                               grouping = 'Stimulation',
                               ncol=2)+ 
  guides(colour = guide_legend(override.aes = list(size=3)))


png("~/RNA-velocity_project/figure_14/figure/figure_14.png", width = 2000, height = 1800, units='px', res = 200)
plot_grid(a, b,c,d, labels=c("A","B","C","D"))
dev.off()
