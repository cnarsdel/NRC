# tcr_tda/nra/core.py

from dataclasses import dataclass
from typing import Any, Callable, Dict, Iterable, List, Optional
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing

import networkx as nx
import numpy as np

from ..tda_tools import compute_persistence
from ..network import basic_network_metrics


@dataclass
class NodeRemovalResult:
    """
    Core result object for NRA (Node Removal Analysis).
    """
    seq: str
    index: int
    diagram_before: np.ndarray
    diagram_after: np.ndarray
    betti_before: np.ndarray
    betti_after: np.ndarray
    delta_betti_sum: float
    network_metrics_before: Optional[Dict[str, Any]] = None
    network_metrics_after: Optional[Dict[str, Any]] = None
    pgen: float = np.nan


def node_removal_analysis(
    seq_to_remove: str,
    dm: np.ndarray,
    node_labels: List[str],
    hom_dim: int = 1,
    graph: Optional[nx.Graph] = None,
    pgen_fn: Optional[Callable[[str], float]] = None,
) -> NodeRemovalResult:
    """
    NRA – Node Removal Analysis for a single node / AA sequence.
    """
    if seq_to_remove not in node_labels:
        raise ValueError(f"Sequence {seq_to_remove!r} not found in node_labels.")

    idx = node_labels.index(seq_to_remove)

    # Topology BEFORE
    diagram_before, betti_before = compute_persistence(dm, hom_dim=hom_dim)

    # Remove node in DM
    dm_after = np.delete(np.delete(dm, idx, axis=0), idx, axis=1)

    # Topology AFTER
    diagram_after, betti_after = compute_persistence(dm_after, hom_dim=hom_dim)

    # Betti change
    L = min(len(betti_before), len(betti_after))
    delta_betti_sum = float(np.sum(betti_after[:L] - betti_before[:L]))

    # Network metrics
    metrics_before = None
    metrics_after = None
    if graph is not None:
        if seq_to_remove not in graph:
            raise ValueError(
                f"Sequence {seq_to_remove!r} not found as node in the provided graph."
            )
        G_before = graph.copy()
        G_after = graph.copy()
        G_after.remove_node(seq_to_remove)

        metrics_before = basic_network_metrics(G_before)
        metrics_after = basic_network_metrics(G_after)

    # Pgen
    if pgen_fn is not None:
        try:
            pgen = float(pgen_fn(seq_to_remove))
        except Exception:
            pgen = float("nan")
    else:
        pgen = float("nan")

    return NodeRemovalResult(
        seq=seq_to_remove,
        index=idx,
        diagram_before=diagram_before,
        diagram_after=diagram_after,
        betti_before=betti_before,
        betti_after=betti_after,
        delta_betti_sum=delta_betti_sum,
        network_metrics_before=metrics_before,
        network_metrics_after=metrics_after,
        pgen=pgen,
    )


# --------------- Parallel batch NRA for many nodes -------------------

def _nra_worker(args) -> NodeRemovalResult:
    (
        seq_to_remove,
        dm,
        node_labels,
        hom_dim,
        graph,
        pgen_fn,
    ) = args
    return node_removal_analysis(
        seq_to_remove=seq_to_remove,
        dm=dm,
        node_labels=node_labels,
        hom_dim=hom_dim,
        graph=graph,
        pgen_fn=pgen_fn,
    )


def node_removal_analysis_batch(
    seqs_to_remove: Iterable[str],
    dm: np.ndarray,
    node_labels: List[str],
    hom_dim: int = 1,
    graph: Optional[nx.Graph] = None,
    pgen_fn: Optional[Callable[[str], float]] = None,
    max_workers: Optional[int] = None,
) -> Dict[str, NodeRemovalResult]:
    """
    Run NRA in parallel for many sequences.

    NOTE: For large distance matrices this can be expensive in both time
    and memory, because each worker recomputes persistence after removal.
    """
    seqs_to_remove = list(seqs_to_remove)
    if not seqs_to_remove:
        return {}

    if max_workers is None:
        max_workers = min(len(seqs_to_remove), max(1, multiprocessing.cpu_count() // 2))

    # dm and graph will be pickled to each worker; on Linux with fork this is
    # relatively cheap, but on Windows it can be memory-heavy.
    tasks = [
        (
            seq,
            dm,
            node_labels,
            hom_dim,
            graph,
            pgen_fn,
        )
        for seq in seqs_to_remove
    ]

    results: Dict[str, NodeRemovalResult] = {}
    with ProcessPoolExecutor(max_workers=max_workers) as exe:
        futures = [exe.submit(_nra_worker, t) for t in tasks]
        for fut in as_completed(futures):
            res: NodeRemovalResult = fut.result()
            results[res.seq] = res

    return results
