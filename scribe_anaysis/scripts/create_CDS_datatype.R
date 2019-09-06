library(Scribe)
library(monocle)

# This script can be used to create Single cell data objects used by Scribe
# Expression data  matrix is required rows are genes columns or cells
# Information about each cell/ sample
# Information about each gene
# names in all three object need to match

#Information about each cell is required
cell_data <- read.csv('~/Desktop/test_data/sample_info.csv', sep=',', header=TRUE, row.names = 1)
#cell_data$Pseudotime <- seq.int(nrow(cell_data))
#random_samples <- sample(rownames(cell_data),size = 2000)

# Expression data  matrix is required rows are genes columns or cells
expr_0 <- read.csv('~/Desktop/test_data/Sx_sz.csv', sep=',', header=TRUE, row.names = 1)

# Gene data, 1 column is required named "gene_short_name"
gene_data <- read.csv('~/Desktop/test_data/gene_info.csv', sep=',', header=TRUE, row.names = 1)
gene_data["gene_short_name"] <- row.names(gene_data)

#Create the Cds with all the new data
cd <- new("AnnotatedDataFrame", data = cell_data)
gd <- new("AnnotatedDataFrame", data = gene_data)
CDS_predicted <- newCellDataSet(as(as.matrix(binded_expr_data), "sparseMatrix"),
                                  phenoData = cd,
                                  expressionFamily = negbinomial.size(),
                                  featureData = gd)

# Estimatate size factors and reducte dimension to infer trajectory
CDS_predicted <- estimateSizeFactors(CDS_predicted)
CDS_predicted <- reduceDimension(CDS_predicted, max_components = 2,
                                 method = 'DDRtree', norm_method = 'none',)

# Order cells based on learned trajecory - assigns states and pseudotime
CDS_predicted <- orderCells(CDS_predicted)

#plot cell trajectory
plot_cell_trajectory(CDS_predicted)
