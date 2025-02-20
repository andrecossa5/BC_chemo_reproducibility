import logging as logg
import numpy as np
import cospar as cs
import scipy.sparse as ssp
from cospar.hf import *
from cospar import hf

adata_orig = cs.datasets.synthetic_bifurcation()


# v1 version, allows to set later time point
def infer_Tmap_from_multitime_clones(
    adata_orig,
    clonal_time_points=None,
    later_time_point=None,
    smooth_array=[15, 10, 5],
    CoSpar_KNN=20,
    sparsity_threshold=0.1,
    intraclone_threshold=0.05,
    normalization_mode=1,
    extend_Tmap_space=False,
    save_subset=True,
    use_full_Smatrix=True,
    trunca_threshold=[0.001, 0.01],
    compute_new=False,
    max_iter_N=5,
    epsilon_converge=0.05,
):
    """
    Infer transition map for clonal data with multiple time points.

    It prepares adata object for cells of targeted time points by
    :func:`cospar.tmap._utils.select_time_points`, generates the similarity matrix
    via :func:`cospar.tmap._utils.generate_similarity_matrix`, and iteratively calls
    the core function :func:`.refine_Tmap_through_cospar` to update
    the transition map.

    * If `later_time_point=None`, the inferred map allows transitions
      between neighboring time points. For example, if
      clonal_time_points=['day1','day2','day3'], then it computes transitions
      for pairs (day1, day2) and (day2, day3), but not (day1, day3).

    * If `later_time_point` is specified, the function produces a map
      between earlier time points and this later time point. For example, if
      `later_time_point='day3`, the map allows transitions for pairs (day1, day3)
      and (day2, day3), but not (day1,day2).

    Parameters
    ------------
    adata_orig: :class:`~anndata.AnnData` object
        Should be prepared from our anadata initialization.
    clonal_time_points: `list` of `str`, optional (default: all time points)
        List of time points to be included for analysis.
        We assume that each selected time point has clonal measurements.
    later_time_points: `list`, optional (default: None)
        If specified, the function will produce a map T between these early
        time points among `clonal_time_points` and the `later_time_point`.
        If not specified, it produces a map T between neighboring clonal time points.
    smooth_array: `list`, optional (default: [15,10,5])
        List of smooth rounds at initial runs of iteration.
        Suppose that it has a length N. For iteration n<N, the n-th entry of
        smooth_array determines the kernel exponent to build the S matrix at the n-th
        iteration. When n>N, we use the last entry of smooth_array to compute
        the S matrix. We recommend starting with more smoothing depth and gradually
        reduce the depth, as inspired by simulated annealing. Data with higher
        clonal dispersion should start with higher smoothing depth. The final depth should
        depend on the manifold itself. For fewer cells, it results in a small KNN graph,
        and a small final depth should be used. We recommend to use a number at
        the multiple of 5 for computational efficiency i.e.,
        smooth_array=[20, 15, 10, 5], or [20,15,10]
    max_iter_N: `int`, optional (default: 5)
        The maximum iterations used to compute the transition map, regardless of epsilon_converge.
    epsilon_converge: `float`, optional (default: 0.05)
        The convergence threshold for the change of map correlations between consecutive iterations.
        This convergence test is activated only when CoSpar has iterated for 3 times.
    CoSpar_KNN: `int`, optional (default: 20)
        The number of neighbors for KNN graph used for computing the
        similarity matrix.
    trunca_threshold: `list`, optional (default: [0.001,0.01])
        Threshold to reset entries of a matrix to zero. The first entry is for
        Similarity matrix; the second entry is for the Tmap.
        This is only for computational and storage efficiency.
    sparsity_threshold: `float`, optional (default: 0.1)
        The relative threshold to remove noises in the updated transition map,
        in the range [0,1].
    intraclone_threshold: `float`, optional (default: 0.05)
        The threshold to remove noises in the demultiplexed (un-smoothed) map,
        in the range [0,1]
    normalization_mode: `int`, optional (default: 1)
        Normalization method. Choice: [0,1].
        0, single-cell normalization; 1, Clone normalization. The clonal
        normalization suppresses the contribution of large
        clones, and is much more robust.
    extend_Tmap_space: `bool` optional (default: `False`)
        If true, the initial states for Tmap will include all states at initial time points,
        and the later states for Tmap will include all states at later time points.
        Otherwise, the initial and later state space of the Tmap will be
        restricted to cells with multi-time clonal information
        alone. The latter case speeds up the computation, which is recommended.
        This option is ignored when `later_time_points` is not None.
    save_subset: `bool`, optional (default: True)
        If true, save only Smatrix at smooth round [5,10,15,...];
        Otherwise, save Smatrix at each round.
    use_full_Smatrix: `bool`, optional (default: True)
        If true, extract the relevant Smatrix from the full Smatrix defined by all cells.
        This tends to be more accurate. The package is optimized around this choice.
    Compute_new: `bool`, optional (default: False)
        If True, compute Smatrix from scratch, whether it was
        computed and saved before or not. This is activated only when
        `use_full_Smatrix=False`.

    Returns
    -------
    adata: :class:`~anndata.AnnData` object
        Store results at adata.uns['transition_map']
        and adata.uns['intraclone_transition_map']. This adata is different
        from the input adata_orig due to subsampling cells.
    """

    t0 = time.time()
    hf.check_available_clonal_info(adata_orig)
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
    if (
        use_full_Smatrix
    ):  # prepare the similarity matrix with all state info, all subsequent similarity will be down-sampled from this one.

        temp_str = "0" + str(trunca_threshold[0])[2:]
        round_of_smooth = np.max(smooth_array)
        data_des = adata_orig.uns["data_des"][0]
        similarity_file_name = os.path.join(
            data_path,
            f"{data_des}_Similarity_matrix_with_all_cell_states_kNN{CoSpar_KNN}_Truncate{temp_str}",
        )
        if not (
            os.path.exists(similarity_file_name + f"_SM{round_of_smooth}.npz")
            and (not compute_new)
        ):
            similarity_matrix_full = tmap_util.generate_similarity_matrix(
                adata_orig,
                similarity_file_name,
                round_of_smooth=round_of_smooth,
                neighbor_N=CoSpar_KNN,
                truncation_threshold=trunca_threshold[0],
                save_subset=save_subset,
                compute_new_Smatrix=compute_new,
            )

    # compute transition map between neighboring time points
    if later_time_point is None:
        logg.info("----Infer transition map between neighboring time points-----")
        logg.info("Step 1: Select time points")
        adata = tmap_util.select_time_points(
            adata_orig,
            time_point=clonal_time_points,
            extend_Tmap_space=extend_Tmap_space,
        )

        logg.info("Step 2: Optimize the transition map recursively")
        tmap_core.infer_Tmap_from_multitime_clones_private(
            adata,
            smooth_array=smooth_array,
            neighbor_N=CoSpar_KNN,
            sparsity_threshold=sparsity_threshold,
            intraclone_threshold=intraclone_threshold,
            normalization_mode=normalization_mode,
            save_subset=save_subset,
            use_full_Smatrix=use_full_Smatrix,
            trunca_threshold=trunca_threshold,
            compute_new_Smatrix=compute_new,
            max_iter_N=max_iter_N,
            epsilon_converge=epsilon_converge,
        )

        if "Smatrix" in adata.uns.keys():
            adata.uns.pop("Smatrix")

        logg.info(f"-----------Total used time: {time.time()-t0} s ------------")
        return adata

    else:
        # compute transition map between initial time points and the later time point
        sel_id = np.nonzero(np.in1d(clonal_time_points, later_time_point))[0][0]
        initial_time_points = clonal_time_points[:sel_id]

        time_info_orig = np.array(adata_orig.obs["time_info"])
        sp_idx = np.zeros(adata_orig.shape[0], dtype=bool)
        all_time_points = list(initial_time_points) + [later_time_point]
        label = "t"
        for xx in all_time_points:
            id_array = np.nonzero(time_info_orig == xx)[0]
            sp_idx[id_array] = True
            label = label + "*" + str(xx)

        adata = adata_orig[sp_idx]
        data_des_orig = adata_orig.uns["data_des"][0]
        data_des_0 = adata_orig.uns["data_des"][-1]
        data_des = (
            data_des_0
            + f"_MultiTimeClone_Later_FullSpace{int(extend_Tmap_space)}_{label}"
        )
        adata.uns["data_des"] = [data_des_orig, data_des]

        time_info = np.array(adata.obs["time_info"])
        time_index_t2 = time_info == later_time_point
        time_index_t1 = ~time_index_t2

        #### used for similarity matrix generation
        Tmap_cell_id_t1 = np.nonzero(time_index_t1)[0]
        Tmap_cell_id_t2 = np.nonzero(time_index_t2)[0]
        adata.uns["Tmap_cell_id_t1"] = Tmap_cell_id_t1
        adata.uns["Tmap_cell_id_t2"] = Tmap_cell_id_t2
        adata.uns["clonal_cell_id_t1"] = Tmap_cell_id_t1
        adata.uns["clonal_cell_id_t2"] = Tmap_cell_id_t2
        adata.uns["sp_idx"] = sp_idx
        data_path = settings.data_path

        transition_map = np.zeros((len(Tmap_cell_id_t1), len(Tmap_cell_id_t2)))
        intraclone_transition_map = np.zeros(
            (len(Tmap_cell_id_t1), len(Tmap_cell_id_t2))
        )

        logg.info(
            "------Infer transition map between initial time points and the later time one------"
        )
        for yy in initial_time_points:

            logg.info(f"--------Current initial time point: {yy}--------")

            logg.info("Step 1: Select time points")
            adata_temp = tmap_util.select_time_points(
                adata_orig, time_point=[yy, later_time_point], extend_Tmap_space=True
            )  # for this to work, we need to set extend_Tmap_space=True, otherwise for different initial time points, the later Tmap_cell_id_t2 may be different

            logg.info("Step 2: Optimize the transition map recursively")
            tmap_core.infer_Tmap_from_multitime_clones_private(
                adata_temp,
                smooth_array=smooth_array,
                neighbor_N=CoSpar_KNN,
                sparsity_threshold=sparsity_threshold,
                intraclone_threshold=intraclone_threshold,
                normalization_mode=normalization_mode,
                save_subset=save_subset,
                use_full_Smatrix=use_full_Smatrix,
                trunca_threshold=trunca_threshold,
                compute_new_Smatrix=compute_new,
                max_iter_N=max_iter_N,
                epsilon_converge=epsilon_converge,
            )

            temp_id_t1 = np.nonzero(time_info == yy)[0]
            sp_id_t1 = hf.converting_id_from_fullSpace_to_subSpace(
                temp_id_t1, Tmap_cell_id_t1
            )[0]

            transition_map[sp_id_t1, :] = adata_temp.uns["transition_map"].A
            intraclone_transition_map[sp_id_t1, :] = adata_temp.uns[
                "intraclone_transition_map"
            ].A

            if "Smatrix" in adata_temp.uns.keys():
                adata_temp.uns.pop("Smatrix")

        adata.uns["transition_map"] = ssp.csr_matrix(transition_map)
        adata.uns["intraclone_transition_map"] = ssp.csr_matrix(
            intraclone_transition_map
        )

        logg.info(f"-----------Total used time: {time.time()-t0} s ------------")
        return adata



# Initial checks
clonal_time_points = ["1", "2"]
later_time_point = None

t0 = time.time()
hf.check_available_clonal_info(adata_orig)
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


smooth_array = [15,10,5]
save_subset = True
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


# Optimize tmap

