# NodeImpact

NodeImpact is a Python package designed for the analysis and processing of DNA sequences, with a focus on utilizing network theory and topological data analysis. It enables the computation of distance matrices, creation of network graphs, generation and analysis of persistence diagrams, and evaluation of the impact of nodes on network structure.

## Features

- **Compute Distance Matrices**: Calculate pairwise distances between sequences using various predefined metrics or user-specified functions.
- **Normalize Distance Matrices**: Standardize distance matrices globally to facilitate consistent comparisons across datasets.
- **Create Network Graphs**: Construct graphs from distance matrices, allowing for visualization and analysis of relationships between sequences.
- **Generate Persistence Diagrams**: Utilize Vietoris-Rips filtration to generate persistence diagrams for topological features in data.
- **Analyze Persistence Diagrams**: Assess entropy and Betti curves to gain insights into the structural properties of persistence diagrams.
- **Node Impact Assessment**: Evaluate the influence of individual nodes on network topology, providing insights into their importance and role.

## Installation

To install the NodeImpact package, first clone the repository to your local machine. Navigate to the root directory of the package, and execute:

```bash
pip install -e .
```
This command installs the package in editable mode, allowing you to modify the source code without having to reinstall.

## Dependencies

The NodeImpact package relies on the integration of several Python libraries that provide robust functionalities for graph operations, numerical analysis, and topological data analysis. Ensure these libraries are installed in your Python environment prior to using NodeImpact:

- **networkx**: A library for the creation, manipulation, and study of the structure, dynamics, and functions of complex networks. It is used in NodeImpact for operations and analyses on network graphs.
  
- **numpy**: A powerful numerical library for Python that supports large, multi-dimensional arrays and matrices, along with a wide collection of mathematical functions. NodeImpact uses numpy for efficient handling of matrix operations and data manipulation.
  
- **gtda**: A Python library known as the `giotto-tda` package designed for topological data analysis operations. It is used primarily for generating persistence diagrams, offering insights into the geometric and topological features within data.
  
- **scipy**: A scientific computing library for Python that provides modules for statistics, optimization, integration, and more. NodeImpact uses scipy for its statistical tests, useful in evaluating node impact within graphs.

These libraries must be installed in your Python environment for NodeImpact to function correctly. You can install them using pip:

```bash
pip install networkx numpy gtda scipy
```

## Usage 

NodeImpact provides a rich toolkit to facilitate the analysis of DNA sequences through network and topological data methods. Below is an outlined example demonstrating the workflow of NodeImpact functionalities:
```python
import numpy as np
from NodeImpact import compute_distances, create_network
from NodeImpact import normalize_distance_matrices
from NodeImpact import PersistenceDiagramGenerator, PersistenceDiagramAnalyzer, NodeImpactAnalyzer
import pwseqdist as pw

# Function to generate random DNA sequences
def generate_random_dna_sequences(num_sequences, sequence_length):
    nucleotides = ['A', 'T', 'C', 'G']
    sequences = [''.join(np.random.choice(nucleotides, sequence_length)) for _ in range(num_sequences)]
    return np.array(sequences)

# Generate 10 random DNA sequences of length 16
random_dna_sequences = generate_random_dna_sequences(10, 16)
print(random_dna_sequences)

# Step 1: Compute the distance matrix for sequences using the Euclidean metric
distance_matrix = compute_distances(sequences, metric=pw.metrics.nw_metric)

# Step 2: Normalize the computed distance matrix across the dataset
normalized_matrices = normalize_distance_matrices([distance_matrix])

# Step 3: Create a network graph from the normalized distance matrix
graph = create_network(sequences, distance_matrix, threshold=1.5)
print(graph.nodes(data=True))

# Step 4: Generate persistence diagrams using topological data analysis techniques
diagram_generator = PersistenceDiagramGenerator(normalized_matrices)
diagrams = diagram_generator.generate_diagrams()

# Step 5: Analyze the persistence diagrams to derive entropy and Betti curves
diagram_analyzer = PersistenceDiagramAnalyzer(diagrams)
entropy = diagram_analyzer.compute_entropy()
betti_curves = diagram_analyzer.compute_betti_curves()

# Step 6: Perform node impact analysis to understand the influence of specific nodes
node_id = 0  # ID of the node in the graph
node_analyzer = NodeImpactAnalyzer(graph, distance_matrix)
node_impact = node_analyzer.analyze_impact(node_id)

#Betticurves, Entropies (0-D,1-D)
print("Node Impact Analysis:", node_impact)
```

The example above follows a comprehensive workflow starting from sequence input, through network graph construction, and persistence diagram analysis, concluding with node impact assessment.