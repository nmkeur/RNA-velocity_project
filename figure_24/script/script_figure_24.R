#load libraries
library(monocle)
library(cowplot)

#load CDS and scribe results
load('~/Desktop/CDS_predicted_monocle_10_genes.RData')

# Create matrix with gene-pairs
gene_pairs_mat = matrix(c("STAT1","IRF1","STAT1","XAF1","STAT1","IFIT3","STAT1","HLA-E",
                          "IRF1","IFIT3","IRF1","IRF2","IRF1","HLA-A"), ncol = 2, byrow = T)


#plot gene pairs in pseudotime + response target visualization

png("~/RNA-velocity_project/figure_24/figure/figure_24.png", width = 1500, height = 1800, units='px', res = 200)
plot_grid(Scribe::plot_gene_pairs_in_pseudotime(CDS_predicted_monocle, gene_pairs_mat,
                                                n_col = 2,
                                                fitting_type = 'loess',
                                                relative_expr = TRUE) +
                                                labs(color='')+
                                                scale_color_manual(labels = c("Target", "Regulator"),
                                                values = c("#F8766D", "#00BFC4")),
          Scribe::plot_lagged_drevi(CDS_predicted_monocle,
                                    gene_pairs_mat,
                                    d=1,
                                    n_col = 2),
          labels = c("A", "B"))
dev.off()
