# Import scaffoldgraph
import scaffoldgraph as sg
import numpy as np
# Import networkx
import networkx as nx

# Import plotting tools
import matplotlib.pyplot as plt

import pandas as pd
# Import rdkit
from rdkit.Chem import Draw
from rdkit import Chem

import random
import os

"""
Set up scaffold network
"""

def setup_network(smiles_dataframe):
    #df = pd.read_csv(smiles_file)
#    network = sg.ScaffoldNetwork.from_dataframe(df, progress=True)
    network = sg.HierS.from_dataframe(df, progress=True)
    
    return network

"""
From a starting smiles, sample successors and predecessors + molecules in those groups. 
1. Sample all molecules with same scaffold
2. Obtain all subscaffolds from query
3. Then sample all molecules with subscaffold architecture (make sure that ) 
"""
def sample_network(network, query_smiles, sampled_smiles, predec_scaffolds, succ_scaffolds, target_size):
    # sampled_smiles.append(query_smiles)
    # query_smiles should already be in sampled_smiles list
    # Just as we can find subscaffolds of a molecule we can find larger scaffolds and molecules from subscaffolds
    #query_mol = Chem.MolFromSmiles(query_smiles)
    # We can differentiate molecules from scaffolds using the node 'type' attribute
    # scaffold nodes have type: 'scaffold' whereas molecules have type: 'molecule'  

    if query_smiles not in predec_scaffolds:
        for pred in network.predecessors(query_smiles):
            if not pred.isnumeric() and pred not in sampled_smiles:
                sampled_smiles.append(pred)
                if len(sampled_smiles)>target_size:
                    break
                if network.nodes[pred]['type'] == 'scaffold':
                    #print(pred)
                    predec_scaffolds.append(pred)   

    if query_smiles not in succ_scaffolds:
        for succ in network.successors(query_smiles):
            if not succ.isnumeric() and succ not in sampled_smiles:
                sampled_smiles.append(succ)
                if len(sampled_smiles)>target_size:
                    break                
                if network.nodes[succ]['type'] == 'scaffold':
                    succ_scaffolds.append(succ) 

    return sampled_smiles, predec_scaffolds, succ_scaffolds


def scaffold_hopping(network, target_size):
    starting_smiles = random.choice(list(network.get_scaffold_nodes()))
    sampled_smiles = []
    sampled_smiles.append(starting_smiles)

    predec_scaffolds = []
    succ_scaffolds = []

    sampled_smiles, predec_scaffolds, succ_scaffolds = sample_network(network, 
                                                                      starting_smiles, 
                                                                      sampled_smiles, 
                                                                      predec_scaffolds, 
                                                                      succ_scaffolds,
                                                                      target_size)

    while len(predec_scaffolds)==0 or len(succ_scaffolds)==0:
        starting_smiles = random.choice(list(network.get_scaffold_nodes()))
        sampled_smiles, predec_scaffolds, succ_scaffolds = sample_network(network, 
                                                                      starting_smiles, 
                                                                      sampled_smiles, 
                                                                      predec_scaffolds, 
                                                                      succ_scaffolds,
                                                                      target_size)
    it = 0
    pred_it = 0
    succ_it = 0
    while len(sampled_smiles)<target_size:
        #print(len(sampled_smiles))
        if it%2==0:
            try:
                #print(predec_scaffolds)
                query_smiles = predec_scaffolds[pred_it]
                sampled_smiles, predec_scaffolds, succ_scaffolds = sample_network(network,
                                                                              query_smiles,
                                                                              sampled_smiles,
                                                                              predec_scaffolds,
                                                                              succ_scaffolds,
                                                                              target_size)
                pred_it+=1
                #print()
                it+=1
            except:
                print("Not enough predecessors")
                starting_smiles = random.choice(list(network.get_scaffold_nodes()))
                sampled_smiles, predec_scaffolds, succ_scaffolds = sample_network(network, 
                                                                      starting_smiles, 
                                                                      sampled_smiles, 
                                                                      predec_scaffolds, 
                                                                      succ_scaffolds,
                                                                      target_size)
                it+=1
                continue
            
        elif it%2==1:
            try:
                #print(succ_scaffolds)

                query_smiles = succ_scaffolds[succ_it]
                #print(query_smiles)
                sampled_smiles, predec_scaffolds, succ_scaffolds = sample_network(network,
                                                                              query_smiles,
                                                                              sampled_smiles,
                                                                              predec_scaffolds,
                                                                              succ_scaffolds,
                                                                              target_size)
                succ_it+=1
                it+=1
            except:
                print("Not enough successors")
                starting_smiles = random.choice(list(network.get_scaffold_nodes()))
                sampled_smiles, predec_scaffolds, succ_scaffolds = sample_network(network, 
                                                                      starting_smiles, 
                                                                      sampled_smiles, 
                                                                      predec_scaffolds, 
                                                                      succ_scaffolds,
                                                                      target_size)
                
                succ_it+=1
                it+=1
                continue

    return sampled_smiles


# Test:
if False:
    smiles_file = 'data/All.sorted.Ena.CACHE.csv'
    df = pd.read_csv(smiles_file)
    df['Smiles']=[str(smi) for smi in df['SMILES']]
    df['Name']=[i for i in range(len(df))]
    df_new = df[['Smiles','Name']]
    df_new.to_csv('data/All.sorted.Ena.CACHE.smi',index=False)
if True:
    df = pd.read_csv('data/All.sorted.Ena.CACHE.smi')
    network = setup_network(df)
    #print(list(network.get_scaffold_nodes()))
    #starting_smiles = list(network.get_scaffold_nodes())[0]
    target_size = 20000

    sampled_smiles = scaffold_hopping(network, target_size)
    #print(set(sampled_smiles))
    print(len(set(sampled_smiles)))
    # list of names
    with open(r'output_smiles_scaffold_search.smi', 'w') as fp:
        fp.write('\n'.join(set(sampled_smiles)))


    # print('Found {} scaffolds in hierarchy 2 containing {}:'.format(len(next_scaffolds), query_smiles)) 

    # mols = [Chem.MolFromSmiles(x) for x in next_scaffolds[:6]]
    # Draw.MolsToGridImage(mols, highlightAtomLists=[mol.GetSubstructMatch(query_mol) for mol in mols])

    # return
