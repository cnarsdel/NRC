# tcr_tda/tda_tools.py

from typing import Dict, Optional, Tuple
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing

import matplotlib.pyplot as plt
import numpy as np
from gtda.homology import VietorisRipsPersistence
from gtda.diagrams import BettiCurve


def compute_persistence(
    dm: np.ndarray,
    hom_dim: int = 1,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute Vietoris–Rips persistence diagram and Betti curve from a precomputed dm.

    Parameters
    ----------
    dm : np.ndarray, shape (n, n)
        Precomputed distance matrix.
    hom_dim : int
        Homology dimension (0=components, 1=loops, ...).

    Returns
    -------
    diagram : np.ndarray
        Persistence diagram array for this dimension.
    betti_curve : np.ndarray
        Betti curve (1D array).
    """
    VR = VietorisRipsPersistence(metric="precomputed", homology_dimensions=[hom_dim])
    BC = BettiCurve()

    diagrams = VR.fit_transform([dm])
    betti_curves = BC.fit_transform(diagrams)

    return diagrams[0], betti_curves[0][0]


def _persistence_worker(args):
    name, dm, hom_dim = args
    dgm, bt = compute_persistence(dm, hom_dim=hom_dim)
    return name, dgm, bt


def compute_persistence_batch(
    matrices: Dict[str, np.ndarray],
    hom_dim: int = 1,
    max_workers: Optional[int] = None,
) -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
    """
    Compute persistence diagrams and Betti curves in parallel for many matrices.

    Parameters
    ----------
    matrices : dict
        Mapping name -> distance matrix.
    hom_dim : int
        Homology dimension.
    max_workers : int or None
        Number of worker processes.

    Returns
    -------
    results : dict
        Mapping name -> (diagram, betti_curve).
    """
    if not matrices:
        return {}

    if max_workers is None:
        max_workers = min(len(matrices), max(1, multiprocessing.cpu_count() // 2))

    results: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}
    with ProcessPoolExecutor(max_workers=max_workers) as exe:
        futures = [
            exe.submit(_persistence_worker, (name, dm, hom_dim))
            for name, dm in matrices.items()
        ]
        for fut in as_completed(futures):
            name, dgm, bt = fut.result()
            results[name] = (dgm, bt)
    return results


# ------------ plotting (usually fine to keep sequential) ------------

def plot_persistence_diagram(
    diagram: np.ndarray,
    ax: Optional[plt.Axes] = None,
    title: str = "Persistence diagram",
    show: bool = True,
):
    if ax is None:
        fig, ax = plt.subplots()

    births = diagram[:, 0]
    deaths = diagram[:, 1]

    max_val = max(np.max(births), np.max(deaths))
    ax.plot([0, max_val], [0, max_val], linestyle="--", linewidth=1)

    ax.scatter(births, deaths, s=20)
    ax.set_xlabel("Birth")
    ax.set_ylabel("Death")
    ax.set_title(title)

    if show:
        plt.show()


def plot_persistence_barcode(
    diagram: np.ndarray,
    ax: Optional[plt.Axes] = None,
    title: str = "Persistence barcode",
    show: bool = True,
):
    if ax is None:
        fig, ax = plt.subplots()

    births = diagram[:, 0]
    deaths = diagram[:, 1]

    y_positions = np.arange(len(births))
    for y, b, d in zip(y_positions, births, deaths):
        ax.hlines(y, b, d, linewidth=2)

    ax.set_xlabel("Filtration value")
    ax.set_yticks([])
    ax.set_title(title)

    if show:
        plt.show()


def plot_betti_curve(
    betti_curve: np.ndarray,
    x: Optional[np.ndarray] = None,
    ax: Optional[plt.Axes] = None,
    title: str = "Betti curve",
    show: bool = True,
):
    if ax is None:
        fig, ax = plt.subplots()

    if x is None:
        x = np.arange(len(betti_curve))

    ax.plot(x, betti_curve)
    ax.set_xlabel("Filtration index")
    ax.set_ylabel("Betti number")
    ax.set_title(title)

    if show:
        plt.show()
