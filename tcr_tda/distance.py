# tcr_tda/distance.py

from typing import Callable, Dict, Optional, Tuple
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing

import numpy as np
import pwseqdist as pw
from tcrdist.rep_funcs import _pw


def compute_distance_matrix(
    seqs: np.ndarray,
    metric: Optional[Callable] = None,
) -> np.ndarray:
    """
    Compute a pairwise distance matrix for a vector of sequences.

    Parameters
    ----------
    seqs : array-like of str (shape: [n])
        Sequences (e.g. AA CDR3).
    metric : callable or None
        A pwseqdist-compatible metric, e.g. pw.metrics.nw_metric.
        If None, default = Needleman–Wunsch (pw.metrics.nw_metric).

    Returns
    -------
    dist_matrix : np.ndarray, shape (n, n)
        Symmetric pairwise distance matrix.
    """
    if metric is None:
        metric = pw.metrics.nw_metric

    dist_matrix = _pw(
        metric=metric,
        seqs1=seqs,
        ncpus=1,        # leave per-matrix parallelism off, use outer parallelization
        uniqify=False,  # keep order aligned with `seqs`
        use_numba=False,
    )
    return dist_matrix


def _dm_worker(args: Tuple[str, np.ndarray, Optional[Callable]]) -> Tuple[str, np.ndarray]:
    """
    Worker for parallel distance computation over multiple datasets.
    """
    name, seqs, metric = args
    dm = compute_distance_matrix(seqs, metric=metric)
    return name, dm


def compute_distance_matrices_parallel(
    seqs_by_dataset: Dict[str, np.ndarray],
    metric: Optional[Callable] = None,
    max_workers: Optional[int] = None,
) -> Dict[str, np.ndarray]:
    """
    Compute distance matrices in parallel for multiple datasets.

    Parameters
    ----------
    seqs_by_dataset : dict
        Mapping dataset name -> 1D array of sequences.
    metric : callable or None
        pwseqdist metric (default: Needleman–Wunsch).
    max_workers : int or None
        Max worker processes. Default: half of CPU cores or len(datasets).

    Returns
    -------
    dms : dict
        Mapping dataset name -> distance matrix.
    """
    if max_workers is None:
        max_workers = min(len(seqs_by_dataset), max(1, multiprocessing.cpu_count() // 2))

    if not seqs_by_dataset:
        return {}

    dms: Dict[str, np.ndarray] = {}
    with ProcessPoolExecutor(max_workers=max_workers) as exe:
        fut2name = {
            exe.submit(_dm_worker, (name, seqs, metric)): name
            for name, seqs in seqs_by_dataset.items()
        }
        for fut in as_completed(fut2name):
            name, dm = fut.result()
            dms[name] = dm

    return dms
