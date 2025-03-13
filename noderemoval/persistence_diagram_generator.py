from gtda.homology import VietorisRipsPersistence

class PersistenceDiagramGenerator:
    def __init__(self, distance_matrices, homology_dimensions=[0, 1]):
        """
        Initialize with a list of distance matrices and homology dimensions.

        Parameters:
        distance_matrices (list of numpy.ndarray): List of distance matrices.
        homology_dimensions (list of int, optional): Homology dimensions (default: [0, 1]).
        """
        self.distance_matrices = distance_matrices
        self.vr = VietorisRipsPersistence(metric='precomputed', homology_dimensions=homology_dimensions)

    def generate_diagrams(self):
        """
        Generate persistence diagrams for each distance matrix.

        Returns:
        list of numpy.ndarray: Persistence diagrams.
        """
        return [self.vr.fit_transform([dm])[0] for dm in self.distance_matrices]