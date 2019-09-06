
import numpy as np
import pandas as pd
import scanpy as sc

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
results_file = './write/results_pilot4.h5ad'  # the file that will store the analysis results
sc.settings.set_figure_params(dpi=300)


adata = sc.read_loom('data/raw_data/merged_notfiltered_pilot4.loom')
adata.var_names_make_unique()  # this is unnecessary if using 'gene_ids'

#AnnData object with n_obs x n_vars = 13292 x 16126
sc.pp.filter_cells(adata, min_genes=100)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.filter_cells(adata, min_genes=100)

# Select all the genes that are mitochrondrial
mito_genes = adata.var_names.str.startswith('MT-')

# for each cell compute fraction of counts in mito genes vs. all genes
# the `.A1` is only necessary as X is sparse (to transform to a dense array after summing)
adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1

# add the total counts per cell as observations-annotation to adata
# for each cell compute fraction of counts in mito genes vs. all genes
# the `.A1` is only necessary as X is sparse (to transform to a dense array after summing)
# add the total counts per cell as observations-annotation to adata
adata.obs['n_counts'] = adata.X.sum(axis=1).A1

#Create violin and scatter plot to determine thresholds used for filtering
sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'],
             jitter=False, multi_panel=True, save='violin_plot_nonfiltered.png')
sc.pl.scatter(adata, x='n_counts', y='percent_mito', save='scatter_counts-vs-mito_nonfiltered.png')
sc.pl.scatter(adata, x='n_counts', y='n_genes', save='scatter-counts-vs-ngenes_nonfiltered.png')

#This perfroms the actual filtering based on the scatter plots above
adata = adata[adata.obs['n_genes'] < 3500, :]
adata = adata[adata.obs['percent_mito'] < 0.05, :]

#Normalize the data.
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
sc.pp.log1p(adata)

#Remove genes that have no expression after we removed the cells
adata = adata[:,adata.X.sum(axis=0) > 0]

sc.pp.regress_out(adata, ['n_counts', 'percent_mito'])

adata.write_loom("data/raw_data/filtered_scaled_pilot4v01.loom")

# #Set the raw data to be the normalized data
# adata.raw = adata
#
#
#
# sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
#
# sc.pl.highly_variable_genes(adata)
#
# adata = adata[:, adata.var['highly_variable']]
#
# sc.pp.regress_out(adata, ['n_counts', 'percent_mito'])
#
# sc.pp.scale(adata, max_value=10)
#
# #Principle component analysis
#
# sc.tl.pca(adata, svd_solver='arpack')
# sc.pl.pca(adata, color='CST3')
#
# sc.pl.pca_variance_ratio(adata, log=True)
#
#
# adata.write(results_file)
#
#
#  #Computes the neighbouring graph
#
#  sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
#
# #Embedding the neighborhood graph
#
# sc.tl.umap(adata)
# sc.pl.umap(adata, color=['CST3', 'NKG7', 'PPBP'])
#
# sc.pl.umap(adata, color=['TNFSF10', 'NKG7', 'PPBP'], use_raw=False)
#
# #Calulate neighbors
# sc.tl.louvain(adata)
#
# sc.pl.umap(adata, color=['louvain'])
#
#
# adata.write(results_file)
#
# sc.tl.rank_genes_groups(adata, 'louvain', method='t-test')
# sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
# sc.settings.verbosity = 2  # reduce the verbosity
#
# sc.tl.rank_genes_groups(adata, 'louvain', method='wilcoxon')
# sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
#
# adata.write(results_file)
#
#
# sc.tl.rank_genes_groups(adata, 'louvain', method='logreg')
# sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
#
# marker_genes = ['IL7R', 'CD79A', 'MS4A1', 'CD8A', 'CD8B', 'LYZ', 'CD14',
#                 'LGALS3', 'S100A8', 'GNLY', 'NKG7', 'KLRB1',
#                 'FCGR3A', 'MS4A7', 'FCER1A', 'CST3', 'PPBP']
#
# adata = sc.read(results_file)
#
# pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(5)
#
#
# result = adata.uns['rank_genes_groups']
# groups = result['names'].dtype.names
# pd.DataFrame(
#     {group + '_' + key[:1]: result[key][group]
#     for group in groups for key in ['names', 'pvals']}).head(5)
#
# sc.tl.rank_genes_groups(adata, 'louvain', groups=['0'], reference='1', method='wilcoxon')
# sc.pl.rank_genes_groups(adata, groups=['0'], n_genes=20)
#
#
#
# sc.pl.rank_genes_groups_violin(adata, groups='0', n_genes=8)
#
#
# adata = sc.read(results_file)
#
# sc.pl.rank_genes_groups_violin(adata, groups='0', n_genes=8)
#
# sc.pl.violin(adata, ['CST3', 'NKG7', 'PPBP'], groupby='louvain')
#
#
# new_cluster_names = [
#     'CD4 T', 'CD14+ Monocytes',
#     'B', 'CD8 T',
#     'NK', 'FCGR3A+ Monocytes',
#     'Dendritic', 'Megakaryocytes']
# adata.rename_categories('louvain', new_cluster_names)
#
#
# sc.pl.umap(adata, color='louvain', legend_loc='on data', title='', frameon=False, save='.pdf')
#
# ax = sc.pl.dotplot(adata, marker_genes, groupby='louvain')
#
# ax = sc.pl.stacked_violin(adata, marker_genes, groupby='louvain', rotation=90)
#
#
# adata.write(results_file, compression='gzip')  # `compression='gzip'` saves disk space, but slows down writing and subsequent reading

#Export
#adata.X = None
#adata.write('./write/pbmc3k_withoutX.h5ad', compression='gzip')
