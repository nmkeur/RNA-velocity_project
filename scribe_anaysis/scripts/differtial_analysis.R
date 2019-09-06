library(Seurat)
library(monocle)
library(VGAM)
library(Scribe)
library(plyr)
library(inflection)
library(igraph)
library(plotrix)
library(parmigene)
library(netbiov)
library(cowplot)

##### PLOT TEST HSMM
#Load Seurat object
CDS_predicted <- readRDS('~/Desktop/CDS_predicted_fulldataset_INF.Rds')
names(pData(CDS_predicted))[names(pData(CDS_predicted)) == 'batch_name'] <- 'Batch'
names(pData(CDS_predicted))[names(pData(CDS_predicted)) == 'stim'] <- 'Stimulation'
CDS_predicted$Batch <- as.character(CDS_predicted$Batch)
names(pData(CDS_predicted))[names(pData(CDS_predicted)) == 'sample'] <- 'Sample'



CDS_predicted <- CDS_predicted[, order(pData(CDS_predicted)$Pseudotime)]

# Perform differtial analysis based on condition
diff_test_results_stimulation <- differentialGeneTest(CDS_predicted,
                                                      fullModelFormulaStr = "~Stimulation")
# Select genes that are significant at an FDR < 10%
sig_genes_stim <- subset(diff_test_results_stimulation, qval < 0.1)
sig_genes_stim <- sig_genes_stim[,c("gene_short_name", "pval", "qval")]


# Perform differtial analysis based on State
diff_test_results_state <- differentialGeneTest(CDS_predicted,
                                                fullModelFormulaStr = "~State")
# Select genes that are significant at an FDR < 10%
sig_genes_state <- subset(diff_test_results_state, qval < 0.01)
sig_genes_state <- sig_genes_state[,c("gene_short_name", "pval", "qval")]


# Perform differtial analysis based on Pseudotime
diff_test_results_pseudotime <- differentialGeneTest(CDS_predicted,
                                                     fullModelFormulaStr = "~sm.ns(Pseudotime)")

# Select genes that are significant at an FDR < 10%
sig_genes_pseudotime <- subset(diff_test_results_pseudotime, qval < 0.01)


intersect(intersect(sig_genes_stim$gene_short_name,sig_genes_state$gene_short_name),sig_genes_pseudotime$gene_short_name)


a <- plot_pseudotime_heatmap(CDS_predicted[sig_genes_state$gene_short_name,],
                        num_clusters = 1,
                        cores = 1,
                        show_rownames = T,
                        cluster_rows = TRUE,
                        scale_max = 3,
                        return_heatmap = TRUE)

plot_genes_branched_heatmap(CDS_predicted[sig_genes_pseudotime$gene_short_name,],
                            num_clusters = 1,
                            cores = 1,
                            show_rownames = T,
                            branch_states = c(2,3),
                            branch_labels = c("State 2", "State 3"),
                            branch_colors = c('#F8766D','#00BA38','#619CFF'),
                            cluster_rows = TRUE,
                            scale_max = 3,
                            return_heatmap = FALSE)

png("~/RNA-velocity_project/suplemenatory_2/figure/suplemenatory_2.png", width = 2000, height = 1800, units='px', res = 200)
plot_grid(a[[4]], labels=c("A"))
dev.off()
