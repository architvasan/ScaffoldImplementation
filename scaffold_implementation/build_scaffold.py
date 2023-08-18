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
    network = sg.HierS.from_dataframe(smiles_dataframe, progress=True)
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

    if query_smiles not in (predec_scaffolds):
        for pred in network.predecessors(query_smiles):
            if not pred.isnumeric() and pred not in sampled_smiles:
                sampled_smiles.append(pred)
                if len(sampled_smiles)>target_size:
                    break
                if network.nodes[pred]['type'] == 'scaffold':
                    #print(pred)
                    predec_scaffolds.append(pred)   

    if query_smiles not in (succ_scaffolds ):
        for succ in network.successors(query_smiles):
            if not succ.isnumeric() and succ not in sampled_smiles:
                sampled_smiles.append(succ)
                if len(sampled_smiles)>target_size:
                    break                
                if network.nodes[succ]['type'] == 'scaffold':
                    succ_scaffolds.append(succ) 

    return sampled_smiles, predec_scaffolds, succ_scaffolds


def _scaffold_hopping_works(start_network, sampling_network, target_size):
    starting_smiles = random.choice(list(start_network.get_scaffold_nodes()))
    sampled_smiles = []
    sampled_smiles.append(starting_smiles)

    predec_scaffolds = []
    succ_scaffolds = []

    sampled_smiles, predec_scaffolds, succ_scaffolds = sample_network(sampling_network, 
                                                                      starting_smiles, 
                                                                      sampled_smiles, 
                                                                      predec_scaffolds, 
                                                                      succ_scaffolds,
                                                                      target_size)

    while len(predec_scaffolds)==0 or len(succ_scaffolds)==0:
        starting_smiles = random.choice(list(start_network.get_scaffold_nodes()))
        sampled_smiles, predec_scaffolds, succ_scaffolds = sample_network(sampling_network, 
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
                sampled_smiles, predec_scaffolds, succ_scaffolds = sample_network(sampling_network,
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
                starting_smiles = random.choice(list(start_network.get_scaffold_nodes()))
                sampled_smiles, predec_scaffolds, succ_scaffolds = sample_network(sampling_network, 
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
                sampled_smiles, predec_scaffolds, succ_scaffolds = sample_network(sampling_network,
                                                                              query_smiles,
                                                                              sampled_smiles,
                                                                              predec_scaffolds,
                                                                              succ_scaffolds,
                                                                              target_size)
                succ_it+=1
                it+=1
            except:
                print("Not enough successors")
                starting_smiles = random.choice(list(start_network.get_scaffold_nodes()))
                sampled_smiles, predec_scaffolds, succ_scaffolds = sample_network(sampling_network, 
                                                                      starting_smiles, 
                                                                      sampled_smiles, 
                                                                      predec_scaffolds, 
                                                                      succ_scaffolds,
                                                                      target_size)
                
                succ_it+=1
                it+=1
                continue

    return sampled_smiles









def scaffold_hopping(start_network, sampling_network, target_per_scaff, num_sample_extend):
    
    start_nodes = list(start_network.get_scaffold_nodes())
    sampled_smiles = []
    predec_scaffolds = []
    succ_scaffolds = []
    for it, query_smi in enumerate(start_nodes):
        sampled_smiles, predec_scaffolds, succ_scaffolds = sample_network(sampling_network,             
                                                                          query_smi, 
                                                                          sampled_smiles, 
                                                                          predec_scaffolds, 
                                                                          succ_scaffolds,       
                                                                          target_per_scaff)

        for it_pre in range(min(len(predec_scaffolds), num_sample_extend)):
            if (it_pre % 10)==0:
                print(it_pre)
            query_predec = random.choice(predec_scaffolds)
            sampled_smiles, predec_scaffolds, succ_scaffolds = sample_network(sampling_network,             
                                                                      query_predec, 
                                                                      sampled_smiles, 
                                                                      predec_scaffolds, 
                                                                      succ_scaffolds,       
                                                                      target_per_scaff)
            predec_scaffolds.remove(query_predec)


        for it_succ in range(min(len(succ_scaffolds), num_sample_extend)):
            if (it_succ % 10)==0:
                print(it_succ)
            query_succ = random.choice(succ_scaffolds)
            sampled_smiles, predec_scaffolds, succ_scaffolds = sample_network(sampling_network,
                                                                      query_succ,
                                                                      sampled_smiles,
                                                                      predec_scaffolds,
                                                                      succ_scaffolds,
                                                                      target_per_scaff)

            succ_scaffolds.remove(query_succ)

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
    df_start = pd.read_csv('data/df_top_150_ena_info.csv')
    df_start['Smiles']=df_start['SMILES']
    df_start['Name']=[i for i in range(len(df_start))]
    df_start_new = df_start[['Smiles', 'Name']]
    df_start_new.to_csv('data/df_top_150_ena.smi')
if True:
    df_sampling = pd.read_csv('data/All.sorted.Ena.CACHE.smi')
    df_start = pd.read_csv('data/df_top_150_ena.smi')
    sampling_network = setup_network(df_sampling)
    start_network = setup_network(df_start)
    #print(list(network.get_scaffold_nodes()))
    #starting_smiles = list(network.get_scaffold_nodes())[0]
    target_size = 20000

    sampled_smiles = scaffold_hopping(start_network=start_network,
                                      sampling_network=sampling_network,
                                        target_per_scaff=100,
                                        num_sample_extend=10
                                      )
    #print(set(sampled_smiles))
    print(len(set(sampled_smiles)))
    # list of names
    with open(r'output_smiles_scaffold_search.smi', 'w') as fp:
        fp.write('\n'.join(set(sampled_smiles)))


    # print('Found {} scaffolds in hierarchy 2 containing {}:'.format(len(next_scaffolds), query_smiles)) 

    # mols = [Chem.MolFromSmiles(x) for x in next_scaffolds[:6]]
    # Draw.MolsToGridImage(mols, highlightAtomLists=[mol.GetSubstructMatch(query_mol) for mol in mols])

    # return
