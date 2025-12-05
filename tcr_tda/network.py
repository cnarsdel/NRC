# tcr_tda/network.py

from typing import Any, Dict, Iterable, List, Optional, Tuple
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing
import warnings

import networkx as nx
import numpy as np


def dm_to_graph(
    dm: np.ndarray,
    nodes: Optional[Iterable[Any]] = None,
    mode: str = "threshold",
    epsilon: Optional[float] = None,
    k: Optional[int] = None,
    include_weights: bool = True,
    warn_if_large: bool = True,
    warn_threshold: int = 2000,
) -> nx.Graph:
    """
    Build a graph from a distance matrix (sequential version).

    For large graphs (> warn_threshold), a warning is emitted.
    """
    n = dm.shape[0]
    if dm.shape[1] != n:
        raise ValueError("Distance matrix must be square.")

    if nodes is None:
        nodes = list(range(n))
    else:
        nodes = list(nodes)
        if len(nodes) != n:
            raise ValueError("Length of `nodes` must match dm.shape[0].")

    if warn_if_large and n > warn_threshold:
        msg = (
            f"Graph has {n} nodes; building and analyzing it may take a while "
            f"and use substantial memory."
        )
        warnings.warn(msg)
        print("[tcr_tda.network] WARNING:", msg)

    G = nx.Graph()
    G.add_nodes_from(nodes)

    if mode == "threshold":
        if epsilon is None:
            raise ValueError("epsilon must be provided for mode='threshold'.")

        for i in range(n):
            for j in range(i + 1, n):
                d = float(dm[i, j])
                if d <= epsilon:
                    if include_weights:
                        G.add_edge(nodes[i], nodes[j], weight=d)
                    else:
                        G.add_edge(nodes[i], nodes[j])

    elif mode == "knn":
        if k is None or k <= 0:
            raise ValueError("k must be a positive integer for mode='knn'.")

        for i in range(n):
            order = np.argsort(dm[i, :])
            neighbors = [idx for idx in order if idx != i][:k]
            for j in neighbors:
                d = float(dm[i, j])
                if include_weights:
                    G.add_edge(nodes[i], nodes[j], weight=d)
                else:
                    G.add_edge(nodes[i], nodes[j])
    else:
        raise ValueError("mode must be 'threshold' or 'knn'.")

    return G


# ----------------- Parallel edge construction (threshold mode) -----------------

def _threshold_edges_for_rows(
    args: Tuple[np.ndarray, List[Any], List[int], float, bool]
) -> List[Tuple[Any, Any, float]]:
    """
    Worker: build edges for a subset of rows in threshold mode.

    Returns list of (u, v, weight) edges.
    """
    dm, nodes, row_indices, epsilon, include_weights = args
    n = dm.shape[0]
    edges: List[Tuple[Any, Any, float]] = []

    for i in row_indices:
        for j in range(i + 1, n):
            d = float(dm[i, j])
            if d <= epsilon:
                if include_weights:
                    edges.append((nodes[i], nodes[j], d))
                else:
                    edges.append((nodes[i], nodes[j], 1.0))  # dummy weight
    return edges


def dm_to_graph_parallel_threshold(
    dm: np.ndarray,
    nodes: Optional[Iterable[Any]] = None,
    epsilon: float = 0.5,
    include_weights: bool = True,
    max_workers: Optional[int] = None,
    warn_if_large: bool = True,
    warn_threshold: int = 2000,
    chunk_size: Optional[int] = None,
) -> nx.Graph:
    """
    Parallelized graph builder for threshold-mode graphs.

    Breaks the rows into chunks and builds edges in parallel. For very
    large n this can significantly speed up edge construction, at the
    cost of higher memory overhead (edges are collected per worker).
    """
    n = dm.shape[0]
    if dm.shape[1] != n:
        raise ValueError("Distance matrix must be square.")

    if nodes is None:
        nodes = list(range(n))
    else:
        nodes = list(nodes)
        if len(nodes) != n:
            raise ValueError("Length of `nodes` must match dm.shape[0].")

    if warn_if_large and n > warn_threshold:
        msg = (
            f"Graph has {n} nodes; building and analyzing it may take a while "
            f"and use substantial memory."
        )
        warnings.warn(msg)
        print("[tcr_tda.network] WARNING:", msg)

    if max_workers is None:
        max_workers = max(1, multiprocessing.cpu_count() // 2)

    if chunk_size is None:
        # heuristic: about #workers chunks, maybe more
        chunk_size = max(64, n // (max_workers * 2) or 1)

    row_indices = list(range(n))
    chunks: List[List[int]] = [
        row_indices[i : i + chunk_size] for i in range(0, n, chunk_size)
    ]

    all_edges: List[Tuple[Any, Any, float]] = []

    with ProcessPoolExecutor(max_workers=max_workers) as exe:
        futures = [
            exe.submit(
                _threshold_edges_for_rows,
                (dm, nodes, chunk, epsilon, include_weights),
            )
            for chunk in chunks
        ]
        for fut in as_completed(futures):
            all_edges.extend(fut.result())

    G = nx.Graph()
    G.add_nodes_from(nodes)
    if include_weights:
        for u, v, w in all_edges:
            G.add_edge(u, v, weight=w)
    else:
        for u, v, _ in all_edges:
            G.add_edge(u, v)

    return G


# ----------------- Basic metrics -----------------

def basic_network_metrics(G: nx.Graph) -> Dict[str, Any]:
    """
    Compute a few simple network metrics.

    Metrics:
    - n_nodes, n_edges
    - average_degree
    - density
    - n_components
    - largest_component_size
    - average_clustering
    """
    n_nodes = G.number_of_nodes()
    n_edges = G.number_of_edges()

    if n_nodes == 0:
        return {
            "n_nodes": 0,
            "n_edges": 0,
            "average_degree": 0.0,
            "density": 0.0,
            "n_components": 0,
            "largest_component_size": 0,
            "average_clustering": 0.0,
        }

    degrees = dict(G.degree())
    avg_degree = float(np.mean(list(degrees.values()))) if degrees else 0.0
    density = nx.density(G)

    components = list(nx.connected_components(G))
    n_components = len(components)
    largest_component_size = max(len(c) for c in components) if components else 0

    try:
        avg_clustering = nx.average_clustering(G)
    except ZeroDivisionError:
        avg_clustering = 0.0

    return {
        "n_nodes": n_nodes,
        "n_edges": n_edges,
        "average_degree": avg_degree,
        "density": density,
        "n_components": n_components,
        "largest_component_size": largest_component_size,
        "average_clustering": avg_clustering,
    }
