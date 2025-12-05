# tcr_tda/nra/__init__.py

"""
NRA – Node Removal Analysis

Tools for removing individual nodes (e.g. amino-acid CDR3s) from a
TCR similarity network and tracking changes in:

- Persistent homology (diagrams + Betti curves)
- Simple network measures (degree, density, components, clustering)
- Optional Pgen (e.g. OLGA human TRB)

Exports:
- NodeRemovalResult: dataclass container for single-node NRA
- node_removal_analysis: single-node NRA
- node_removal_analysis_batch: parallel NRA for many nodes
"""

from .core import (
    NodeRemovalResult,
    node_removal_analysis,
    node_removal_analysis_batch,
)

__all__ = [
    "NodeRemovalResult",
    "node_removal_analysis",
    "node_removal_analysis_batch",
]
