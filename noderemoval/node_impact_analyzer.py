import numpy as np
from gtda.homology import VietorisRipsPersistence
from gtda.diagrams import PersistenceEntropy, BettiCurve
from scipy.stats import chisquare
import networkx as nx

class NodeImpactAnalyzer:
    def __init__(self, graph, distance_matrix, homology_dimensions=[0, 1]):
        """
        Initialize with a graph, distance matrix, and homology dimensions.

        Parameters:
        graph (networkx.Graph): Graph representation of sequences.
        distance_matrix (numpy.ndarray): Distance matrix of sequences.
        homology_dimensions (list of int, optional): Homology dimensions (default: [0, 1]).
        """
        self.graph = graph
        self.distance_matrix = distance_matrix
        self.vr = VietorisRipsPersistence(metric='precomputed', homology_dimensions=homology_dimensions)
        self.pe = PersistenceEntropy(normalize=True)
        self.bc = BettiCurve()

    def analyze_impact(self, node_id):
        """
        Analyze persistence diagrams with and without a specific node for Betti curves and entropy changes.

        Parameters:
        node_id (int): The ID of the node to analyze.

        Returns:
        tuple: Tuple containing Betti curves without node, Betti curves with node, entropy with node, entropy without node, chi-square value for Betti curves, and entropy change.
        """
        # Diagrams with node
        diagram_with_node = self.vr.fit_transform([self.distance_matrix])[0]

        # Diagrams without node
        subgraph = self.graph.copy()
        subgraph.remove_node(node_id)
        sub_dm = np.delete(np.delete(self.distance_matrix, node_id, axis=0), node_id, axis=1)
        diagram_without_node = self.vr.fit_transform([sub_dm])[0]

        # Compute Betti curves and entropy
        betti_with_node = self.bc.fit_transform([diagram_with_node])[0][0]
        betti_with_node1 = self.bc.fit_transform([diagram_with_node])[0][1]
        betti_without_node = self.bc.fit_transform([diagram_without_node])[0][0]
        betti_without_node1 = self.bc.fit_transform([diagram_without_node])[0][1]
        entropy_with_node = self.pe.fit_transform([diagram_with_node])[0][0]
        entropy_with_node1 = self.pe.fit_transform([diagram_with_node])[0][1]
        entropy_without_node = self.pe.fit_transform([diagram_without_node])[0][0]
        entropy_without_node1 = self.pe.fit_transform([diagram_without_node])[0][1]
        # Normalize Betti curves
        betti_with_node =betti_with_node/np.sum(betti_with_node)
        betti_without_node =betti_without_node/np.sum(betti_without_node)

        # Chi-square calculations
        betti_chi_square =( np.sum((betti_with_node- betti_without_node)**2),np.sum((betti_with_node1- betti_without_node1)**2))
        entropy_change = (entropy_with_node - entropy_without_node , entropy_with_node1 - entropy_without_node1)

        return betti_chi_square, entropy_change