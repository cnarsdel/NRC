from gtda.diagrams import PersistenceEntropy, BettiCurve

class PersistenceDiagramAnalyzer:
    def __init__(self, diagrams):
        """
        Initialize the analyzer with a list of persistence diagrams.

        Parameters:
        diagrams (list of numpy.ndarray): List of persistence diagrams.
        """
        self.diagrams = diagrams

    def compute_entropy(self):
        """
        Compute the entropy of persistence diagrams.

        Returns:
        numpy.ndarray: Entropy values.
        """
        pe = PersistenceEntropy(normalize=True)
        return pe.fit_transform(self.diagrams)

    def compute_betti_curves(self):
        """
        Compute Betti curves of persistence diagrams.

        Returns:
        numpy.ndarray: Betti curves.
        """
        bc = BettiCurve()
        return bc.fit_transform(self.diagrams)