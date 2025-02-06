"""
Vignette cospar for clonal trajectory.
"""

import os
import scanpy as sc
import cospar as cs
from plotting_utils._plotting_base import *
matplotlib.use('macOSX')

path = '/Users/IEO5505/Desktop/BC_chemo_reproducibility/downstream/MDA/cospar'
cs.logging.print_version()
cs.settings.data_path = (path)
cs.settings.figure_path = (path)

# fig, axs = plt.subplots(1,2,figsize=(8,3.5))
# sc.pl.embedding(adata_orig, basis='X_emb', color='state_info', ax=axs[0], show=False, frameon=False)
# sc.pl.embedding(adata_orig, basis='X_emb', color='time_info', ax=axs[1], show=False, frameon=False)
# fig.tight_layout()
# plt.show()

# cs.pl.barcode_heatmap(adata_orig)
# plt.show()

# cs.pl.clones_on_manifold(adata_orig, selected_clone_list=[1], color_list=adata_orig.uns['time_info_colors'])
# plt.show()

# cs.tl.fate_coupling(adata_orig, source="X_clone",normalize=False)
# cs.pl.fate_coupling(adata_orig, source="X_clone")
# plt.show()

# cs.tl.clonal_fate_bias(adata_orig, selected_fate="Fate_A", alternative="two-sided")
# cs.pl.clonal_fate_bias(adata_orig)
# plt.show()

# adata = cs.tmap.infer_Tmap_from_multitime_clones(
#     adata_orig,
#     clonal_time_points=["1", "2"],
#     smooth_array=[20, 15, 10],
#     CoSpar_KNN=10,
#     sparsity_threshold=0.1,
# )

# adata.uns['intraclone_transition_map'].A.sum()
# 
# 
# test_1 = (adata.obsm['X_clone'].A[:,0]==1) & (adata.obs['time_info'] == '1')
# test_2 = (adata.obsm['X_clone'].A[:,0]==1) & (adata.obs['time_info'] == '2')
# 
# 
# adata.uns['Tmap_cell_id_t1']
# adata.uns['intraclone_transition_map'].A[test_1, test_2]


import os
import time

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as ssp

from cospar.tmap import _tmap_core as tmap_core
from cospar.tmap import _utils as tmap_util

from cospar import hf 
import logging as logg
import cospar.settings as settings
import cospar.tool as tl


adata_orig = cs.datasets.synthetic_bifurcation()
adata = adata_orig.copy()










later_time_point = None
clonal_time_points = ['1','2']
smooth_array = [20,15,10]
save_subset = False
extend_Tmap_space = False


# cs.tmap.infer_Tmap_from_multitime_clones
t0 = time.time()
# hf.check_available_clonal_info(adata)

# check_available_clonal_info
X_clone = adata.obsm["X_clone"]
time_info = adata.obs["time_info"]

hf.update_time_ordering(adata, mode="auto")

# record time points with clonal information
if ssp.issparse(X_clone):
    clone_N_per_cell = X_clone.sum(1).A.flatten()
else:
    clone_N_per_cell = X_clone.sum(1)

clonal_time_points = []
for xx in list(set(time_info)):
    idx = np.array(time_info) == xx
    if np.sum(clone_N_per_cell[idx]) > 0:
        clonal_time_points.append(xx)
time_ordering = adata.uns["time_ordering"]
sel_idx_temp = np.in1d(time_ordering, clonal_time_points)
clonal_time_points = time_ordering[sel_idx_temp]
adata.uns["clonal_time_points"] = clonal_time_points



clonal_time_points_0 = np.array(adata_orig.uns["clonal_time_points"])
if len(clonal_time_points_0) < 2:
    raise ValueError("There are no multi-time clones. Abort the inference.")
if clonal_time_points is None:
    clonal_time_points = clonal_time_points_0
if type(later_time_point) == list:
    later_time_point = later_time_point[0]
if later_time_point is not None:
    clonal_time_points = list(clonal_time_points) + [later_time_point]
    clonal_time_points = list(set(clonal_time_points))

hf.check_input_parameters(
    adata_orig,
    later_time_point=later_time_point,
    clonal_time_points=clonal_time_points,
    smooth_array=smooth_array,
    save_subset=save_subset,
)

# order the clonal time points
time_ordering = adata_orig.uns["time_ordering"]
sel_idx_temp = np.in1d(time_ordering, clonal_time_points)
clonal_time_points = time_ordering[sel_idx_temp]
logg.info("------Compute the full Similarity matrix if necessary------")
data_path = settings.data_path






if later_time_point is None:

    logg.info("----Infer transition map between neighboring time points-----")
    logg.info("Step 1: Select time points")

    adata = tmap_util.select_time_points(
        adata_orig,
        time_point=clonal_time_points,
        extend_Tmap_space=extend_Tmap_space,
    )







def select_time_points(
    adata_orig, time_point=["day_1", "day_2"], extend_Tmap_space=False
):
    """
    Select barcoded cells at given time points for Tmap inference.

    Select cells at given time points, and prepare the right data structure
    for running core cospar function to infer the Tmap.

    Parameters
    ----------
    adata_orig: original :class:`~anndata.AnnData` object
    time_point: `list` optional (default: ['day_1','day_2'])
        Require at least two time points, arranged in ascending order.
    extend_Tmap_space: `bool` optional (default: `False`)
        If true, the initial states for Tmap will include all states at initial time points,
        and the later states for Tmap will include all states at later time points.
        Otherwise, the initial and later state
        space of the Tmap will be restricted to cells with multi-time clonal information
        alone. The latter case speeds up the computation, which is recommended.

    Returns
    -------
    Subsampled :class:`~anndata.AnnData` object
    """

    # x_emb_orig=adata_orig.obsm['X_emb'][:,0]
    # y_emb_orig=adata_orig.obsm['X_emb'][:,1]
    time_info_orig = np.array(adata_orig.obs["time_info"])
    clone_annot_orig = adata_orig.obsm["X_clone"]
    if len(time_point) == 0:  # use all clonally labelled cell states
        time_point = np.sort(
            list(set(time_info_orig))
        )  # this automatic ordering may not work

    if len(time_point) < 2:
        logg.error("Must select more than 1 time point!")
    else:

        At = []                                     # Clonal matrices subsets 
        for j, time_0 in enumerate(time_point):
            At.append(ssp.csr_matrix(clone_annot_orig[time_info_orig == time_0]))

        ### Day t - t+1
        Clonal_cell_ID_FOR_t = []
        for j in range(len(time_point) - 1):
            idx_t = np.array((At[j] * At[j + 1].T).sum(1) > 0).flatten()
            time_index_t = time_info_orig == time_point[j]
            temp = np.nonzero(time_index_t)[0][idx_t]
            Clonal_cell_ID_FOR_t.append(
                temp
            )  # this index is in the original space, without sampling etc

            logg.hint(
                f"Clonal cell fraction (day {time_point[j]}-{time_point[j+1]}):",
                len(temp) / np.sum(time_index_t),
            )

        ### Day t+1 - t
        Clonal_cell_ID_BACK_t = []
        for j in range(len(time_point) - 1):
            idx_t = np.array((At[j + 1] * At[j].T).sum(1) > 0).flatten()
            time_index_t = time_info_orig == time_point[j + 1]
            temp = np.nonzero(time_index_t)[0][idx_t]
            Clonal_cell_ID_BACK_t.append(
                temp
            )  # this index is in the original space, without sampling etc

            logg.hint(
                f"Clonal cell fraction (day {time_point[j+1]}-{time_point[j]}):",
                len(temp) / np.sum(time_index_t),
            )

        for j in range(len(time_point) - 1):
            logg.hint(
                f"Numer of cells that are clonally related -- day {time_point[j]}: {len(Clonal_cell_ID_FOR_t[j])}  and day {time_point[j+1]}: {len(Clonal_cell_ID_BACK_t[j])}"
            )

        proportion = np.ones(len(time_point))
        # flatten the list
        flatten_clonal_cell_ID_FOR = np.array(
            [sub_item for item in Clonal_cell_ID_FOR_t for sub_item in item]
        )
        flatten_clonal_cell_ID_BACK = np.array(
            [sub_item for item in Clonal_cell_ID_BACK_t for sub_item in item]
        )
        valid_clone_N_FOR = np.sum(
            clone_annot_orig[flatten_clonal_cell_ID_FOR].A.sum(0) > 0
        )
        valid_clone_N_BACK = np.sum(
            clone_annot_orig[flatten_clonal_cell_ID_BACK].A.sum(0) > 0
        )

        logg.info(f"Number of multi-time clones post selection: {valid_clone_N_FOR}")
        # logg.info("Valid clone number 'BACK' post selection",valid_clone_N_BACK)

        ###################### select initial and later cell states

        if extend_Tmap_space:
            old_Tmap_cell_id_t1 = []
            for t_temp in time_point[:-1]:
                old_Tmap_cell_id_t1 = old_Tmap_cell_id_t1 + list(
                    np.nonzero(time_info_orig == t_temp)[0]
                )
            old_Tmap_cell_id_t1 = np.array(old_Tmap_cell_id_t1)

            ########
            old_Tmap_cell_id_t2 = []
            for t_temp in time_point[1:]:
                old_Tmap_cell_id_t2 = old_Tmap_cell_id_t2 + list(
                    np.nonzero(time_info_orig == t_temp)[0]
                )
            old_Tmap_cell_id_t2 = np.array(old_Tmap_cell_id_t2)

        else:
            old_Tmap_cell_id_t1 = flatten_clonal_cell_ID_FOR
            old_Tmap_cell_id_t2 = flatten_clonal_cell_ID_BACK

        old_clonal_cell_id_t1 = flatten_clonal_cell_ID_FOR
        old_clonal_cell_id_t2 = flatten_clonal_cell_ID_BACK
        ########################

        sp_id = np.sort(
            list(set(list(old_Tmap_cell_id_t1) + list(old_Tmap_cell_id_t2)))
        )
        sp_idx = np.zeros(clone_annot_orig.shape[0], dtype=bool)
        sp_idx[sp_id] = True

        Tmap_cell_id_t1 = hf.converting_id_from_fullSpace_to_subSpace(
            old_Tmap_cell_id_t1, sp_id
        )[0]
        clonal_cell_id_t1 = hf.converting_id_from_fullSpace_to_subSpace(
            old_clonal_cell_id_t1, sp_id
        )[0]
        clonal_cell_id_t2 = hf.converting_id_from_fullSpace_to_subSpace(
            old_clonal_cell_id_t2, sp_id
        )[0]
        Tmap_cell_id_t2 = hf.converting_id_from_fullSpace_to_subSpace(
            old_Tmap_cell_id_t2, sp_id
        )[0]

        Clonal_cell_ID_FOR_t_new = []
        for temp_id_list in Clonal_cell_ID_FOR_t:
            convert_list = hf.converting_id_from_fullSpace_to_subSpace(
                temp_id_list, sp_id
            )[0]
            Clonal_cell_ID_FOR_t_new.append(convert_list)

        Clonal_cell_ID_BACK_t_new = []
        for temp_id_list in Clonal_cell_ID_BACK_t:
            convert_list = hf.converting_id_from_fullSpace_to_subSpace(
                temp_id_list, sp_id
            )[0]
            Clonal_cell_ID_BACK_t_new.append(convert_list)

        sp_id_0 = np.sort(list(old_clonal_cell_id_t1) + list(old_clonal_cell_id_t2))
        sp_idx_0 = np.zeros(clone_annot_orig.shape[0], dtype=bool)
        sp_idx_0[sp_id_0] = True

        barcode_id = np.nonzero(clone_annot_orig[sp_idx_0].A.sum(0).flatten() > 0)[0]
        # sp_id=np.nonzero(sp_idx)[0]
        clone_annot = clone_annot_orig[sp_idx][:, barcode_id]

        adata = adata_orig[sp_idx]
        adata.obsm["X_clone"] = clone_annot
        adata.uns["clonal_cell_id_t1"] = clonal_cell_id_t1
        adata.uns["clonal_cell_id_t2"] = clonal_cell_id_t2
        adata.uns["Tmap_cell_id_t1"] = Tmap_cell_id_t1
        adata.uns["Tmap_cell_id_t2"] = Tmap_cell_id_t2
        adata.uns["multiTime_cell_id_t1"] = Clonal_cell_ID_FOR_t_new
        adata.uns["multiTime_cell_id_t2"] = Clonal_cell_ID_BACK_t_new
        adata.uns["proportion"] = np.ones(len(time_point) - 1)
        adata.uns["sp_idx"] = sp_idx

        data_des_orig = adata_orig.uns["data_des"][0]
        data_des_0 = adata_orig.uns["data_des"][-1]
        time_label = "t"
        for x in time_point:
            time_label = time_label + f"*{x}"

        data_des = (
            data_des_0
            + f"_MultiTimeClone_FullSpace{int(extend_Tmap_space)}_{time_label}"
        )
        adata.uns["data_des"] = [data_des_orig, data_des]

        if logg._settings_verbosity_greater_or_equal_than(3):
            N_cell, N_clone = clone_annot.shape
            logg.info(f"Cell number={N_cell}, Clone number={N_clone}")
            x_emb = adata.obsm["X_emb"][:, 0]
            y_emb = adata.obsm["X_emb"][:, 1]
            pl.customized_embedding(x_emb, y_emb, -x_emb)

        logg.hint(f"clonal_cell_id_t1: {len(clonal_cell_id_t1)}")
        logg.hint(f"Tmap_cell_id_t1: {len(Tmap_cell_id_t1)}")
        
        return adata
    








time_point = ["1", "2"]
extend_Tmap_space = False



























































































# cs.hf.check_available_map(adata)
# adata.uns['available_map']

# cs.tl.fate_bias(adata, selected_fates=["Fate_A", "Fate_B"], source="transition_map")
# cs.pl.fate_bias(adata, source="transition_map")
# plt.show()
# 
# cs.tl.fate_bias(adata, selected_fates=["Fate_A", "Fate_B"], source="intraclone_transition_map")
# cs.pl.fate_bias(adata, source="intraclone_transition_map")
# plt.show()

# selected_state_id_list = [2]
# map_backward = True
# cs.pl.single_cell_transition(
#     adata,
#     selected_state_id_list=selected_state_id_list,
#     source="intraclone_transition_map",
#     map_backward=map_backward,
# )
# plt.show()

# cs.tl.fate_potency(
#     adata,
#     source="transition_map",
#     map_backward=True,
#     method="norm-sum",
#     fate_count=True,
# )
# cs.pl.fate_potency(adata, source="transition_map")
# plt.show()

# adata = cs.tmap.infer_Tmap_from_clonal_info_alone(adata_orig)
# cs.tl.fate_bias(
#     adata, selected_fates=["Fate_A", "Fate_B"], source="clonal_transition_map"
# )
# cs.pl.fate_bias(adata, source="clonal_transition_map")
# plt.show()

# sc.pp.neighbors(adata, n_pcs=10, n_neighbors=15)
# sc.tl.draw_graph(adata, n_jobs=-1)
# fig, axs = plt.subplots(1,2,figsize=(8,3.5))
# sc.pl.embedding(adata_orig, basis='X_emb', color='time_info', ax=axs[1], show=False, frameon=False)
# sc.pl.embedding(adata, basis='X_draw_graph_fa', color='time_info', show=False, ax=axs[0], frameon=False)
# fig.tight_layout()
# plt.show()
