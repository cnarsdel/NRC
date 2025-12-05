# tcr_tda/pipeline.py

#!/usr/bin/env python3
"""
General TCR topology pipeline with parallelization.

- Parallel distance matrices across datasets.
- Parallel TDA across datasets.
- Optional graphs, with warnings for >2000 nodes.
- Optional bulk node-removal in parallel (per dataset).
"""

import glob
import multiprocessing
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Callable, Dict, List, Optional, Tuple

import networkx as nx
import numpy as np
import pandas as pd

from .distance import compute_distance_matrix, compute_distance_matrices_parallel
from .network import dm_to_graph, dm_to_graph_parallel_threshold, basic_network_metrics
from .tda_tools import compute_persistence, compute_persistence_batch

# OLGA (optional)
try:
    import olga
    import olga.load_model as load_model
    import olga.generation_probability as gp_mod
except ImportError:
    olga = None
    load_model = None
    gp_mod = None


def build_olga_pgen_human_TRB() -> Callable[[str], float]:
    """
    Build a callable: AA_CDR3 -> Pgen (float) using OLGA's human TRB default model.
    """
    if olga is None or load_model is None or gp_mod is None:
        raise ImportError("OLGA and its submodules are required to build Pgen.")

    model_path = os.path.join(olga.__path__[0], "default_models", "human_T_beta")
    params_file = os.path.join(model_path, "model_params.txt")
    marginals_file = os.path.join(model_path, "model_marginals.txt")
    V_anchor_file = os.path.join(model_path, "V_gene_CDR3_anchors.csv")
    J_anchor_file = os.path.join(model_path, "J_gene_CDR3_anchors.csv")

    genomic_data = load_model.GenomicDataVDJ()
    genomic_data.load_igor_genomic_data(params_file, V_anchor_file, J_anchor_file)

    generative_model = load_model.GenerativeModelVDJ()
    generative_model.load_and_process_igor_model(marginals_file)

    pgen_model = gp_mod.GenerationProbabilityVDJ(generative_model, genomic_data)
    return pgen_model.compute_aa_CDR3_pgen


def stem(path: str) -> str:
    return os.path.splitext(os.path.basename(path))[0]


def _node_removal_worker(args):
    """
    Worker to compute TDA after removing one node/index.
    Returns: (node_seq, betti_curve_after_removal, diagram_after_removal)
    """
    from .tda_tools import compute_persistence  # local import for multiprocessing

    M, idx, node_seq, hom_dim = args
    M_sub = np.delete(np.delete(M, idx, axis=0), idx, axis=1)
    diagram_sub, betti_sub = compute_persistence(M_sub, hom_dim=hom_dim)
    return node_seq, betti_sub, diagram_sub


def run_pipeline(
    data_glob: str,
    out_dir_dm: str = "networks",
    out_dir_tda: str = "Topology",
    max_rows_per_dataset: Optional[int] = None,
    hom_dim: int = 1,
    cdr3_col: str = "cdr3aa",
    normalize_distances: bool = True,
    dataset_name_filter: Optional[Callable[[str], bool]] = None,
    # distance metric
    distance_metric: Optional[Callable] = None,
    # Graph/network options
    build_graphs: bool = False,
    graph_mode: str = "threshold",   # "threshold" or "knn"
    graph_epsilon: Optional[float] = None,
    graph_k: Optional[int] = None,
    warn_large_graph: bool = True,
    warn_graph_threshold: int = 2000,
    parallel_graph_threshold: bool = False,
    # Specificity / bulk node-removal
    specificity_filename_fn: Optional[Callable[[str], Optional[str]]] = None,
    specificity_cdr3_col: str = "CDR3b",
    specificity_ag_col: str = "Ag_gene",
    specificity_treatment_col: Optional[str] = None,
    treatment_resolver: Optional[Callable[[str], str]] = None,
    # Pgen builder (optional)
    pgen_builder: Optional[Callable[[], Callable[[str], float]]] = build_olga_pgen_human_TRB,
    # Parallelization controls
    max_workers_dm: Optional[int] = None,
    max_workers_tda: Optional[int] = None,
    max_workers_node_removal: Optional[int] = None,
) -> None:
    """
    General TCR topology pipeline with parallelization.
    """
    os.makedirs(out_dir_dm, exist_ok=True)
    os.makedirs(out_dir_tda, exist_ok=True)

    # --------------------- Pgen ---------------------
    if pgen_builder is not None:
        try:
            pgen_fn = pgen_builder()
        except Exception as e:
            print("Warning: OLGA Pgen not available, setting pgen=NaN. Error:", e)
            pgen_fn = lambda s: np.nan
    else:
        pgen_fn = lambda s: np.nan

    # --------------------- Load data ---------------------
    file_paths = sorted(glob.glob(data_glob))
    if not file_paths:
        raise FileNotFoundError(f"No input files matched {data_glob!r}")

    name_by_path = {fp: stem(fp) for fp in file_paths}
    file_map: Dict[str, str] = {}
    names: List[str] = []

    for fp in file_paths:
        name = name_by_path[fp]
        if dataset_name_filter is not None and not dataset_name_filter(name):
            continue
        names.append(name)
        file_map[name] = fp

    if not names:
        raise RuntimeError("No datasets selected after applying dataset_name_filter.")

    dfs: Dict[str, pd.DataFrame] = {}
    for name in names:
        df = pd.read_csv(file_map[name], sep="\t")
        if max_rows_per_dataset is not None:
            df = df.iloc[:max_rows_per_dataset].copy()
        if cdr3_col not in df.columns:
            raise KeyError(f"{name}: missing {cdr3_col!r} column")
        df[cdr3_col] = df[cdr3_col].astype(str)
        dfs[name] = df

    print("Datasets loaded:", names)

    seqs_per_dataset: Dict[str, np.ndarray] = {
        name: df[cdr3_col].to_numpy(dtype=object) for name, df in dfs.items()
    }

    # --------------------- Distance matrices (parallel) ---------------------
    if max_workers_dm is None:
        max_workers_dm = min(len(seqs_per_dataset), max(1, multiprocessing.cpu_count() // 2))

    print(f"Computing distance matrices in parallel with {max_workers_dm} worker(s)...")
    try:
        DMS = compute_distance_matrices_parallel(
            seqs_per_dataset,
            metric=distance_metric,
            max_workers=max_workers_dm,
        )
    except Exception as e:
        print("Parallel distance computation failed, falling back to sequential:", e)
        DMS = {
            name: compute_distance_matrix(seqs, metric=distance_metric)
            for name, seqs in seqs_per_dataset.items()
        }

    DMS = {name: DMS[name] for name in names}

    # --------------------- Normalize distances ---------------------
    if normalize_distances:
        mins = [float(np.min(m)) for m in DMS.values()]
        maxs = [float(np.max(m)) for m in DMS.values()]
        global_min, global_max = min(mins), max(maxs)
        if global_max == global_min:
            NDMS: Dict[str, np.ndarray] = {
                name: np.zeros_like(m, dtype=float) for name, m in DMS.items()
            }
        else:
            denom = (global_max - global_min)
            NDMS = {name: (m - global_min) / denom for name, m in DMS.items()}
    else:
        NDMS = DMS

    NODE_NAMES: Dict[str, List[str]] = {
        name: dfs[name][cdr3_col].astype(str).tolist() for name in names
    }

    # Save all matrices + node labels in one NPZ
    save_dict = {}
    for name, mat in NDMS.items():
        save_dict[f"{name}_matrix"] = mat
    for name, nodes_list in NODE_NAMES.items():
        save_dict[f"{name}_nodes"] = np.asarray(nodes_list, dtype=object)

    dm_npz = os.path.join(out_dir_dm, "distance_matrices_and_nodes_all.npz")
    np.savez_compressed(dm_npz, **save_dict)
    print("Saved global distance matrix NPZ:", dm_npz)

    # --------------------- TDA on full matrices (parallel) ---------------------
    if max_workers_tda is None:
        max_workers_tda = min(len(NDMS), max(1, multiprocessing.cpu_count() // 2))

    print(f"Running TDA in parallel with {max_workers_tda} worker(s)...")
    tda_results = compute_persistence_batch(NDMS, hom_dim=hom_dim, max_workers=max_workers_tda)

    dgms_full: Dict[str, np.ndarray] = {}
    bts_full: Dict[str, np.ndarray] = {}
    for name in names:
        diagram, betti = tda_results[name]
        dgms_full[name] = diagram
        bts_full[name] = betti

        np.savez_compressed(
            os.path.join(out_dir_tda, f"diagrams_{name}_dimension_{hom_dim}.npz"),
            diagram=diagram,
        )
        np.savez_compressed(
            os.path.join(out_dir_tda, f"betti_curves_{name}_dimension_{hom_dim}.npz"),
            betti_curve=betti,
        )
        print(f"TDA saved for {name}")

    # --------------------- Graphs & network metrics ---------------------
    network_metrics_rows: List[Dict[str, object]] = []
    graphs_by_name: Dict[str, nx.Graph] = {}

    if build_graphs:
        print("Building graphs for each dataset...")
        for name in names:
            dm = NDMS[name]
            nodes = NODE_NAMES[name]
            if graph_mode == "threshold" and parallel_graph_threshold:
                G = dm_to_graph_parallel_threshold(
                    dm,
                    nodes=nodes,
                    epsilon=(graph_epsilon if graph_epsilon is not None else 0.5),
                    include_weights=True,
                    max_workers=max(1, multiprocessing.cpu_count() // 2),
                    warn_if_large=warn_large_graph,
                    warn_threshold=warn_graph_threshold,
                )
            else:
                G = dm_to_graph(
                    dm,
                    nodes=nodes,
                    mode=graph_mode,
                    epsilon=graph_epsilon,
                    k=graph_k,
                    warn_if_large=warn_large_graph,
                    warn_threshold=warn_graph_threshold,
                )
            graphs_by_name[name] = G

            graph_path = os.path.join(out_dir_dm, f"graph_{name}.gpickle")
            nx.write_gpickle(G, graph_path)
            print(f"Graph saved for {name} -> {graph_path}")

            metrics = basic_network_metrics(G)
            metrics["dataset"] = name
            network_metrics_rows.append(metrics)

        if network_metrics_rows:
            df_metrics = pd.DataFrame(network_metrics_rows)
            csv_path = os.path.join(out_dir_dm, "network_metrics_summary.csv")
            df_metrics.to_csv(csv_path, index=False)
            print("Network metrics summary saved:", csv_path)

    # --------------------- Bulk node-removal ---------------------
    if specificity_filename_fn is None:
        print("No specificity_filename_fn provided; skipping bulk node-removal.")
        return

    if max_workers_node_removal is None:
        max_workers_node_removal = max(1, multiprocessing.cpu_count() // 2)

    for name in names:
        print(f"\nBulk node-removal analysis for dataset {name}")

        if treatment_resolver is not None:
            treatment = treatment_resolver(name)
        else:
            treatment = "all"

        spec_path = specificity_filename_fn(name)
        if spec_path is None:
            print(f"  Skipping {name}: specificity_filename_fn returned None")
            continue
        if not os.path.exists(spec_path):
            print(f"  Warning: missing specificity file {spec_path}; skipping {name}.")
            continue

        spec_df = pd.read_csv(spec_path)

        sp = spec_df
        if (
            specificity_treatment_col is not None
            and specificity_treatment_col in spec_df.columns
            and treatment in ("pre", "post")
        ):
            sp = spec_df[
                spec_df[specificity_treatment_col].astype(str).str.contains(treatment, na=False)
            ]

        required_cols = {specificity_cdr3_col, specificity_ag_col}
        if not required_cols.issubset(sp.columns):
            print(
                f"  Warning: specificity for {name} missing columns "
                f"{required_cols - set(sp.columns)}; skipping."
            )
            continue

        nodes_list = NODE_NAMES[name]
        node_to_idx = {s: i for i, s in enumerate(nodes_list)}
        common_nodes = list(set(sp[specificity_cdr3_col].astype(str)).intersection(nodes_list))
        if not common_nodes:
            print("  No overlap between dataset nodes and specificity CDR3s; skipping.")
            continue

        M = NDMS[name]
        base_bts = bts_full[name]

        tasks = [
            (M, node_to_idx[s], s, hom_dim)
            for s in common_nodes
            if s in node_to_idx
        ]

        print(
            f"  Node-removal: {len(tasks)} overlapping nodes, "
            f"{max_workers_node_removal} worker(s)."
        )

        noderemoval_bc: Dict[str, np.ndarray] = {}
        noderemoval_dgm: Dict[str, np.ndarray] = {}
        noderemoval_rows: List[Tuple[str, float, float, str]] = []

        with ProcessPoolExecutor(max_workers=max_workers_node_removal) as exe:
            for node_seq, bts_sub, dgms_sub in exe.map(_node_removal_worker, tasks):
                noderemoval_bc[node_seq] = bts_sub
                noderemoval_dgm[node_seq] = dgms_sub

                L = min(len(base_bts), len(bts_sub))
                delta = float(np.sum(bts_sub[:L] - base_bts[:L]))

                try:
                    pg = float(pgen_fn(node_seq))
                except Exception:
                    pg = np.nan

                subset = sp.loc[sp[specificity_cdr3_col].astype(str) == node_seq, specificity_ag_col]
                ag_gene = subset.iloc[0] if len(subset) else np.nan
                noderemoval_rows.append((node_seq, delta, pg, ag_gene))

        np.savez_compressed(
            os.path.join(out_dir_tda, f"node_removal_betti_curves_{name}_dimension_{hom_dim}.npz"),
            **noderemoval_bc,
        )
        np.savez_compressed(
            os.path.join(out_dir_tda, f"node_removal_diagrams_{name}_dimension_{hom_dim}.npz"),
            **noderemoval_dgm,
        )

        noderemoval_df = pd.DataFrame(
            noderemoval_rows,
            columns=["CDR3", "delta_betti_sum", "pgen", "Ag_gene"],
        )
        out_csv = os.path.join(out_dir_tda, f"node_removal_{name}_dimension_{hom_dim}.csv")
        noderemoval_df.to_csv(out_csv, index=False)
        print("  Bulk node-removal summary saved:", out_csv)
