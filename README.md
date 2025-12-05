# tcr_tda – TCR Topology & Node Removal Analysis (NRA)

`tcr_tda` is a small toolkit for **T cell receptor (TCR)** repertoire analysis that focuses on:

- Building **sequence similarity networks** from TCR repertoires using configurable distance metrics
- Extracting **topological signatures** with Vietoris–Rips persistence and Betti curves (via Giotto-TDA)
- Computing **simple network metrics** (degree, density, components, clustering)
- Performing **NRA – Node Removal Analysis**:
  - Remove a specific CDR3 node (AA sequence)
  - Measure how persistent homology and network measures change
  - Optionally track generation probability (Pgen, via OLGA)

The package is written with **parallelization** in mind for:
- Distance matrix computation across datasets
- TDA across datasets
- Graph construction (threshold mode)
- Bulk node removal (per dataset)

> ⚠️ For graphs with more than ~2000 nodes, building and analyzing the graph can be slow and memory-intensive. The code will print and emit a warning when this happens.

---

## Features

- **Distance matrices**
  - `compute_distance_matrix` – per-dataset pairwise distance (default: Needleman–Wunsch)
  - `compute_distance_matrices_parallel` – parallel computation across many datasets

- **Network construction**
  - `dm_to_graph` – distance matrix → graph (threshold or kNN)
  - `dm_to_graph_parallel_threshold` – parallel edge construction for threshold graphs
  - `basic_network_metrics` – simple metrics: degree, density, components, clustering

- **Topology (TDA)**
  - `compute_persistence` – Vietoris–Rips persistence diagram + Betti curve for one matrix
  - `compute_persistence_batch` – parallel TDA across many matrices
  - Plot helpers: `plot_persistence_diagram`, `plot_persistence_barcode`, `plot_betti_curve`

- **NRA – Node Removal Analysis**
  - `node_removal_analysis` – remove a single node (AA CDR3), track:
    - Persistence diagrams & Betti curves before/after
    - Simple network metrics before/after
    - Pgen (if a function is provided)
  - `node_removal_analysis_batch` – parallel NRA for many sequences
  - Result container: `NodeRemovalResult` (dataclass)

- **End-to-end pipeline**
  - `run_pipeline` – repertoire files → distance matrices → TDA → graphs → bulk node removal
  - Optional OLGA integration for Pgen: `build_olga_pgen_human_TRB`

- **CLI**
  - `tcr-tda run` – run the full pipeline on files
  - `tcr-tda nra` – single-node NRA on a saved dataset

---

## Installation

### Requirements

- Python 3.9+
- Recommended scientific stack:
  - `numpy`, `pandas`, `networkx`, `matplotlib`
  - `pwseqdist`, `tcrdist3`
  - `giotto-tda`
- Optional:
  - `olga` (for generation probabilities, Pgen)

### Install in editable mode

Clone your repo, then from the project root:

```bash
pip install -e .
