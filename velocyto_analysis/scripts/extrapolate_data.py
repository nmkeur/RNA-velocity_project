# import libraries
import velocyto as vcy
import pandas as pd

# path to loom file
vlm = vcy.VelocytoLoom("data/raw_data/merged_filtered_pilot4.loom")



# Normalize spliced and unspliced counts seperatly
vlm._normalize_S(relative_size=vlm.S.sum(0),
             target_size=vlm.S.sum(0).mean())
vlm._normalize_U(relative_size=vlm.U.sum(0),
             target_size=vlm.U.sum(0).mean())

#pd.DataFrame(vlm.S_norm)
#pd.DataFrame(vlm.Sx_sz_t)

# Perform pca and knn smooting of the data
vlm.perform_PCA()
vlm.knn_imputation(n_pca_dims=25, k=25, balanced=True, b_sight=150, b_maxl=150, n_jobs=16)

# Use linear regression to fit gamma's
vlm.fit_gammas()

# Predict unspliced and calculate velocity and extrapolate
vlm.predict_U()
vlm.calculate_velocity()
vlm.calculate_shift(assumption="constant_velocity")
vlm.extrapolate_cell_at_t(delta_t=1.)

# pandas is used to create dataframe and write the matrices to csv
pd.DataFrame(vlm.Ux_sz, index=vlm.ra['var_names'], columns=vlm.ca['obs_names']).to_csv('ux_sz.csv', sep=',')
pd.DataFrame(vlm.Upred, index=vlm.ra['var_names'], columns=vlm.ca['obs_names']).to_csv('ux_sz_pred.csv', sep=',')
pd.DataFrame(vlm.velocity, index=vlm.ra['var_names'], columns=vlm.ca['obs_names']).to_csv('ux_sz_velocity.csv', sep=',')

pd.DataFrame(vlm.S_norm, index=vlm.ra['var_names'], columns=vlm.ca['obs_names']).to_csv('Sx_sz.csv', sep=',')
pd.DataFrame(vlm.Sx_sz_t, index=vlm.ra['var_names'], columns=vlm.ca['obs_names']).to_csv('Sx_sz_pred.csv', sep=',')

#Filtering is done using boealean numpy arrays
vlm.filter_cells(cell_array_type)
vlm.S_norm = vlm.S_norm[:,cell_array_type]
vlm.S_sz = vlm.S_sz[:,cell_array_type]
vlm.U_norm = vlm.U_norm[:,cell_array_type]
vlm.U_sz = vlm.U_sz[:,cell_array_type]
vlm.Unorm_factor = vlm.Unorm_factor[cell_array_type]
vlm.Ucell_size = vlm.Ucell_size[cell_array_type]

vlm.filter_cells(cell_array_cond)
vlm.S_norm = vlm.S_norm[:,cell_array_cond]
vlm.S_sz = vlm.S_sz[:,cell_array_cond]
vlm.U_norm = vlm.U_norm[:,cell_array_cond]
vlm.U_sz = vlm.U_sz[:,cell_array_cond]
vlm.Unorm_factor = vlm.Unorm_factor[cell_array_cond]
vlm.Ucell_size = vlm.Ucell_size[cell_array_cond]
