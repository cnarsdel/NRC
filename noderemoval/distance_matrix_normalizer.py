import numpy as np

def normalize_distance_matrices(distance_matrices):
    """
    Normalize the distance matrices based on global min and max values.

    Parameters:
    distance_matrices (list of numpy.ndarray): List of distance matrices.

    Returns:
    list of numpy.ndarray: Normalized distance matrices.
    """
    all_values = np.concatenate([dm.flatten() for dm in distance_matrices])
    global_min = np.min(all_values)
    global_max = np.max(all_values)

    return [(dm - global_min) / (global_max - global_min) for dm in distance_matrices]