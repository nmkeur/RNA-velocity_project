#Load required libraries
library(monocle)
library(Scribe)
library(cowplot)
library(reshape2)

#Load CDS file used in analysis
CDS_predicted <- readRDS('~/Desktop/CDS_CD4_INF_3600.Rds')
CDS_predicted <- CDS_predicted[, order(pData(CDS_predicted)$Pseudotime)]

# Create dataframe using expr matrx
exp_matrix <- as.data.frame(t((as.matrix((exprs(CDS_predicted))))))
exp_matrix$ID <- seq.int(nrow(exp_matrix))

# Extract projection data
exp_matrix$Projection <- CDS_predicted$Projection
exp_matrix$Projection[exp_matrix$Projection == 0] <- 'Measured'
exp_matrix$Projection[exp_matrix$Projection == 1] <- "Projected"


#Plot dataframe
plot_grid(ggplot(exp_matrix, aes(x=ID, y=IRF1)) +
            geom_point(size=2, aes(color=Projection))+xlab('Pseudotime')+ ylab('Raw Expression')+ggtitle('IRF1')+theme(plot.title = element_text(hjust = 0.5)),
          ggplot(exp_matrix, aes(x=ID, y=STAT1)) +
            geom_point(size=2, aes(color=Projection))+ xlab('Pseudotime')+ ylab('Raw Expression')+ggtitle('STAT1')+theme(plot.title = element_text(hjust = 0.5)),
          labels = c("A","B"),
          nrow = 2)
