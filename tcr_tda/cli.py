# tcr_tda/cli.py

import argparse
import os

import numpy as np
import pwseqdist as pw

from .pipeline import run_pipeline, build_olga_pgen_human_TRB
from .network import dm_to_graph, dm_to_graph_parallel_threshold
from .nra import node_removal_analysis


def _resolve_metric(name: str):
    if name == "nw":
        return pw.metrics.nw_metric
    raise ValueError(f"Unknown metric {name!r}; currently supported: 'nw'")


def main(argv=None):
    parser = argparse.ArgumentParser(prog="tcr_tda", description="TCR topology toolkit")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # --------------------------- run ---------------------------
    run_p = subparsers.add_parser("run", help="Run the full pipeline (parallelized)")
    run_p.add_argument("--data-glob", required=True, help="Glob for input repertoires")
    run_p.add_argument("--out-dm", default="networks", help="Output directory for distance matrices/graphs")
    run_p.add_argument("--out-tda", default="Topology", help="Output directory for TDA artifacts")
    run_p.add_argument("--hom-dim", type=int, default=1, help="Homology dimension for TDA")
    run_p.add_argument("--cdr3-col", default="cdr3aa", help="Column name for AA CDR3")
    run_p.add_argument("--metric", default="nw", help="Distance metric (currently only 'nw')")
    run_p.add_argument("--build-graphs", action="store_true", help="Build graphs from distance matrices")
    run_p.add_argument("--graph-mode", default="threshold", choices=["threshold", "knn"], help="Graph construction mode")
    run_p.add_argument("--epsilon", type=float, help="Distance threshold for mode='threshold'")
    run_p.add_argument("--k", type=int, help="k for mode='knn'")
    run_p.add_argument("--parallel-graph-threshold", action="store_true", help="Parallelize threshold graph building")
    run_p.add_argument("--max-rows", type=int, help="Cap rows per dataset")
    run_p.add_argument("--no-olga", action="store_true", help="Disable OLGA Pgen")
    run_p.add_argument("--max-workers-dm", type=int, help="Workers for distance matrices")
    run_p.add_argument("--max-workers-tda", type=int, help="Workers for TDA")
    run_p.add_argument("--max-workers-node-removal", type=int, help="Workers for bulk node removal")

    # --------------------------- NRA ---------------------------
    nr_p = subparsers.add_parser(
        "nra",
        help="NRA – Node Removal Analysis (single-node removal)",
    )
    nr_p.add_argument("--dm-npz", required=True, help="NPZ file with matrices and nodes")
    nr_p.add_argument("--dataset", required=True, help="Dataset name prefix in NPZ (e.g. 'pt1_at')")
    nr_p.add_argument("--seq", required=True, help="AA CDR3 sequence to remove")
    nr_p.add_argument("--hom-dim", type=int, default=1, help="Homology dimension for TDA")
    nr_p.add_argument("--graph-mode", default="threshold", choices=["threshold", "knn"], help="Graph construction mode")
    nr_p.add_argument("--epsilon", type=float, help="Distance threshold for mode='threshold'")
    nr_p.add_argument("--k", type=int, help="k for mode='knn'")
    nr_p.add_argument("--use-olga", action="store_true", help="Use OLGA for Pgen")
    nr_p.add_argument("--parallel-graph-threshold", action="store_true", help="Parallelize threshold graph building")
    nr_p.add_argument("--out-csv", help="Optional CSV to append an NRA summary row")

    args = parser.parse_args(argv)

    # --------------------------- run ---------------------------
    if args.command == "run":
        metric_fn = _resolve_metric(args.metric)

        if args.no_olga:
            pgen_builder = None
        else:
            pgen_builder = build_olga_pgen_human_TRB

        run_pipeline(
            data_glob=args.data_glob,
            out_dir_dm=args.out_dm,
            out_dir_tda=args.out_tda,
            max_rows_per_dataset=args.max_rows,
            hom_dim=args.hom_dim,
            cdr3_col=args.cdr3_col,
            distance_metric=metric_fn,
            build_graphs=args.build_graphs,
            graph_mode=args.graph_mode,
            graph_epsilon=args.epsilon,
            graph_k=args.k,
            parallel_graph_threshold=args.parallel_graph_threshold,
            pgen_builder=pgen_builder,
            max_workers_dm=args.max_workers_dm,
            max_workers_tda=args.max_workers_tda,
            max_workers_node_removal=args.max_workers_node_removal,
        )
        return

    # --------------------------- NRA ---------------------------
    if args.command == "nra":
        npz = np.load(args.dm_npz, allow_pickle=True)
        matrix_key = f"{args.dataset}_matrix"
        nodes_key = f"{args.dataset}_nodes"

        if matrix_key not in npz or nodes_key not in npz:
            raise KeyError(
                f"Could not find keys {matrix_key!r} / {nodes_key!r} in {args.dm_npz!r}"
            )

        dm = npz[matrix_key]
        nodes = list(npz[nodes_key])

        # Build graph (possibly parallel in threshold mode)
        if args.graph_mode == "threshold" and args.parallel_graph_threshold:
            G = dm_to_graph_parallel_threshold(
                dm,
                nodes=nodes,
                epsilon=(args.epsilon if args.epsilon is not None else 0.5),
                include_weights=True,
                max_workers=None,
            )
        else:
            G = dm_to_graph(
                dm,
                nodes=nodes,
                mode=args.graph_mode,
                epsilon=args.epsilon,
                k=args.k,
            )

        pgen_fn = None
        if args.use_olga:
            try:
                pgen_fn = build_olga_pgen_human_TRB()
            except Exception as e:
                print("Warning: OLGA Pgen not available; Pgen will be NaN. Error:", e)
                pgen_fn = None

        result = node_removal_analysis(
            seq_to_remove=args.seq,
            dm=dm,
            node_labels=nodes,
            hom_dim=args.hom_dim,
            graph=G,
            pgen_fn=pgen_fn,
        )

        print(f"Dataset: {args.dataset}")
        print(f"Sequence: {result.seq}")
        print(f"Index in matrix: {result.index}")
        print(f"Delta Betti sum (H{args.hom_dim}): {result.delta_betti_sum}")
        print(f"Pgen: {result.pgen}")

        if result.network_metrics_before is not None:
            print("\nNetwork metrics BEFORE removal:")
            for k, v in result.network_metrics_before.items():
                print(f"  {k}: {v}")
        if result.network_metrics_after is not None:
            print("\nNetwork metrics AFTER removal:")
            for k, v in result.network_metrics_after.items():
                print(f"  {k}: {v}")

        if args.out_csv:
            import pandas as pd

            row = {
                "dataset": args.dataset,
                "seq": result.seq,
                "index": result.index,
                "delta_betti_sum": result.delta_betti_sum,
                "pgen": result.pgen,
            }
            if result.network_metrics_before:
                for k, v in result.network_metrics_before.items():
                    row[f"before_{k}"] = v
            if result.network_metrics_after:
                for k, v in result.network_metrics_after.items():
                    row[f"after_{k}"] = v

            if os.path.exists(args.out_csv):
                df = pd.read_csv(args.out_csv)
                df = pd.concat([df, pd.DataFrame([row])], ignore_index=True)
            else:
                df = pd.DataFrame([row])
            df.to_csv(args.out_csv, index=False)
            print("\nNRA summary row written to:", args.out_csv)
