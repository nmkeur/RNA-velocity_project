#Load required libraries
library(monocle)
library(Scribe)
library(cowplot)

#Load CDS file used in analysis
CDS_predicted <- readRDS('~/Desktop/CDS_CD4_INF_1800.Rds')

#Change the variable names for figures
names(pData(CDS_predicted))[names(pData(CDS_predicted)) == 'batch_name'] <- 'Batch'
names(pData(CDS_predicted))[names(pData(CDS_predicted)) == 'stim'] <- 'Stimulation'
CDS_predicted$Batch <- as.character(CDS_predicted$Batch)
names(pData(CDS_predicted))[names(pData(CDS_predicted)) == 'sample'] <- 'Sample'

plot_genes_in_pseudotime(CDS_predicted[c('IRF1','STAT1')],
                         color_by = 'Stimulation', min_expr = 0.3)

# In order to calculate the RDI scores for TF only
# we can use the wired function. however, this function
# needs to be manually edited to get it funtion.
# colnames need to be adjusted to rownames and 
# the exprssion matrix needs to be changed to the CDS
# object. Uniformatization can be enabled optional.
trace(wired, edit=TRUE)

tt1_1800 <- wired(cds = CDS_predicted,
                  TF = c("IRF1","STAT1"),
                  informative_genes =  c("IFIT3","HLA-A","HLA-E","IFIT1","IRF3","IFI6","IRF2","XAF1","IRF1","STAT1"),
                  delays = 1,
                  include_conditioning = FALSE,
                  smoothing = TRUE,
                  cluster_TFs = FALSE,
                  cluster_targets = FALSE)

plot_genes_in_pseudotime(CDS_predicted[c('IRF1','STAT1')],
                         color_by = 'Stimulation', min_expr = 0.3)

pheatmap::pheatmap(tt1_1800$RDI_res$max_rdi_value, cluster_rows = F, cluster_cols = F,
                   annotation_names_row = T, border_color = NA, fontsize_row = 6)



