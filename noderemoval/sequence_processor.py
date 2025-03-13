import numpy as np
import networkx as nx

def compute_distances(sequences, metric):
    """
    Compute the distance matrix using the provided metric.

    Parameters:
    sequences (numpy.ndarray): Array of DNA sequences.
    metric (str or callable): Predefined metric name ('euclidean', 'manhattan', etc.) or custom function.

    Returns:
    numpy.ndarray: A matrix of distances between sequences.
    """
    metric_function = validate_metric(metric)
    n = len(sequences)
    dist_matrix = np.zeros((n, n))

    for i in range(n):
        for j in range(i + 1, n):
            dist_matrix[i, j] = dist_matrix[j, i] = metric_function(sequences[i], sequences[j])
    return dist_matrix

def validate_metric(metric):
    """
    Validate and set the distance metric.

    Parameters:
    metric (str or callable): Metric for distance calculation.

    Returns:
    callable: The function to compute the distance between two sequences.
    """
    if isinstance(metric, str):
        if metric == 'euclidean':
            return lambda x, y: np.linalg.norm(x - y)
        elif metric == 'manhattan':
            return lambda x, y: np.sum(np.abs(x - y))
        else:
            raise ValueError("Unsupported metric name provided.")
    elif callable(metric):
        return metric
    else:
        raise TypeError("Metric must be either a string or a callable function.")

def create_network(sequences, distance_matrix, threshold=None):
    """
    Construct a graph using the computed distance matrix with an optional threshold.

    Parameters:
    sequences (numpy.ndarray): Array of DNA sequences.
    distance_matrix (numpy.ndarray): Distance matrix.
    threshold (float, optional): Distance threshold for edge creation.

    Returns:
    networkx.Graph: Graph with sequences as nodes and edges based on distance threshold.
    """
    G = nx.Graph()

    # Add nodes with sequences as attributes
    for i, seq in enumerate(sequences):
        G.add_node(i, sequence=seq)

    # Add edges
    for i in range(len(distance_matrix)):
        for j in range(i + 1, len(distance_matrix)):
            if threshold is None or distance_matrix[i, j] < threshold:
                G.add_edge(i, j, weight=distance_matrix[i, j])

    return G