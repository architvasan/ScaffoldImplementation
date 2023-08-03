# Import scaffoldgraph
import scaffoldgraph as sg
import numpy as np
# Import networkx
import networkx as nx

# Import plotting tools
import matplotlib.pyplot as plt

# Import rdkit
from rdkit.Chem import Draw
from rdkit import Chem

import random
import os

"""
Set up scaffold network
"""

def setup_network(smiles_file):
    network = sg.ScaffoldNetwork.from_smiles(smiles_file, progress=True)
    return network

"""
From a starting smiles, sample successors and predecessors + molecules in those groups. 
1. Sample all molecules with same scaffold
2. Obtain all subscaffolds from query
3. Then sample all molecules with subscaffold architecture (make sure that ) 
"""
def sample_network(network, query_smiles, number_steps):
    sampled_smiles = []
    sampled_smiles.append(query_smiles)
    
    return