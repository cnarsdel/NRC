# tcr_tda/__init__.py

"""
tcr_tda – TCR Topology & Node Removal Analysis (NRA)

Main exposed API:
- run_pipeline: end-to-end repertoire → distance matrices → TDA → graphs
- distance: compute_distance_matrix / compute_distance_matrices_parallel
- network: dm_to_graph, dm_to_graph_parallel_threshold, basic_network_metrics
- TDA tools: compute_persistence(+batch) and plotting helpers
- NRA: node_removal_analysis (+batch) with NodeRemovalResult
"""

from .pipeline import run_pipeline, build_olga_pgen_human_TRB

from .distance import (
    compute_distance_matrix,
    compute_distance_matrices_parallel,
)

from .network import (
    dm_to_graph,
    dm_to_graph_parallel_threshold,
    basic_network_metrics,
)

from .tda_tools import (
    compute_persistence,
    compute_persistence_batch,
    plot_persistence_diagram,
    plot_persistence_barcode,
    plot_betti_curve,
)

from .nra import (
    NodeRemovalResult,
    node_removal_analysis,
    node_removal_analysis_batch,
)

__all__ = [
    # pipeline
    "run_pipeline",
    "build_olga_pgen_human_TRB",
    # distance
    "compute_distance_matrix",
    "compute_distance_matrices_parallel",
    # network
    "dm_to_graph",
    "dm_to_graph_parallel_threshold",
    "basic_network_metrics",
    # TDA
    "compute_persistence",
    "compute_persistence_batch",
    "plot_persistence_diagram",
    "plot_persistence_barcode",
    "plot_betti_curve",
    # NRA
    "NodeRemovalResult",
    "node_removal_analysis",
    "node_removal_analysis_batch",
]
