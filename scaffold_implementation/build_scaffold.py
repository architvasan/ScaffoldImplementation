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
def sample_network(network, query_smiles, sampled_smiles, predec_scaffolds, succ_scaffolds):
    # sampled_smiles.append(query_smiles)
    # query_smiles should already be in sampled_smiles list
    # Just as we can find subscaffolds of a molecule we can find larger scaffolds and molecules from subscaffolds
    #query_mol = Chem.MolFromSmiles(query_smiles)

    # First lets find scaffolds in the above hierarchy (-> 2)
    # We can differentiate molecules from scaffolds using the node 'type' attribute
    # scaffold nodes have type: 'scaffold' whereas molecules have type: 'molecule'  

    if query_smiles not in predec_scaffolds:
        for pred in network.predecessors(query_smiles):
            if pred not in sampled_smiles:
                sampled_smiles.append(pred)
                if network.nodes[pred]['type'] == 'scaffold':
                    predec_scaffolds.append(pred)   

    succ_scaffolds = []
    for succ in network.successors(query_smiles):
        if succ not in sampled_smiles:
            sampled_smiles.append(succ)
            if network.nodes[succ]['type'] == 'scaffold':
                succ_scaffolds.append(succ) 

    return sampled_smiles, predec_scaffolds, succ_scaffolds

def iterate():
    return



    print('Found {} scaffolds in hierarchy 2 containing {}:'.format(len(next_scaffolds), query_smiles)) 

    mols = [Chem.MolFromSmiles(x) for x in next_scaffolds[:6]]
    Draw.MolsToGridImage(mols, highlightAtomLists=[mol.GetSubstructMatch(query_mol) for mol in mols])
    
    return