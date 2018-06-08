import numpy as np
import random
from ete2 import Tree
import time
import subprocess
import resource
random.seed(time.time())
seed=random.randint(0,100000000)

mutlist=[0.0000001,0.0000002,0.0000003,0.0000004,0.0000005,0.0000006,0.0000007,0.0000008,0.0000008,0.0000009,0.000001,0.000002,0.000003,0.000004,0.000005,0.000006,0.000007,0.000008,0.000009,0.00001,0.00002,0.00003,0.00004,0.00005,0.00006,0.00007,0.00008,0.00009,0.0001,0.0002,0.0003,0.0004,0.0005,0.0006,0.0007,0.0008,0.0009,0.001]
idx_mut=np.random.randint(0,len(mutlist)-1)
mut_rate=mutlist[idx_mut]

fitlist=[0.0001,0.0002,0.0003,0.0004,0.0005,0.0006,0.0007,0.0008,0.0009,0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.009,0.01]
idx_fit=np.random.randint(0,len(fitlist)-1)
fitness=fitlist[idx_fit]

maxcell=500000000
maxgene=20000
maxdriver=500

Tumor=Tree()
bastard=Tumor.get_tree_root()
bastard.name="bastard"


bastard.add_features(cell=1,driver=1,passenger=0,mutation="0\t") # number of cells in the node, number of driver in each cell, number of passengers in each cells, label of the mutations

tot_population=1
itera=0
generation=1
countnode=1

while tot_population >0 and tot_population<maxcell:

    if (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000 ) > 400000000:
        quit()
    countnode=0
    tot_population=0

    countnode=0
    G=[]
    for node in Tumor.traverse(strategy="postorder"): # traverse tree from the leaf, trick to add new child without disturbing the loop
        nb_driver=int(node.driver) # number of drivers per cell in the node
        d=0.5*np.power((1.-fitness),nb_driver) # probability to die
        max_replicate=int(node.cell)
        nb_replicate=np.random.binomial(max_replicate,1-d,1) # number of cells that replicates with a success probabilitu 1-d given max_replicate can replicate

        if nb_replicate == 0: 
            node.cell=0
            G.append(node) #store the node for further deleting with chidren connected to the next possible parent

        else:
            nb_new_clones=np.random.binomial(nb_replicate,mut_rate,1) # among nb_replicate, nb_new_clones appear, each with a probability mut_rate
            node.cell=int(2*nb_replicate-nb_new_clones) # add to the current node, the number of cells that do not aquire a new mutations
            tot_population=tot_population+int(node.cell)
            mut_list=node.mutation 
            
            for i in range(nb_new_clones):
                
                num_gene=np.random.randint(1,maxgene) # gene that is mutated
                nb_passenger=int(node.passenger) # number of passenger per cell in the node
                new_clone=node.add_child() # add a new clone
                new_clone.name=str(generation)
                generation=generation+1

                
                if num_gene <= maxdriver : # the gene is a driver, corresponding to X% of the total number of gene
                    new_clone.add_features(cell=1,driver=nb_driver+1,passenger=nb_passenger,mutation=mut_list+"%d\t"%num_gene)
                else: # the mutated gene is a passenger
                    new_clone.add_features(cell=1,driver=nb_driver,passenger=nb_passenger+1,mutation=mut_list+"%d\t"%num_gene)

        countnode=countnode+1

    for node in G: # delete nodes with no replicate
        node.delete(prevent_nondicotomic=False)
    itera=itera+1



finalnode=[]

totcell=0
for node in Tumor.traverse("postorder"):
    totcell=totcell+int(node.cell)

if totcell > 1:
    # CHECK THE PATH
    tree_fn="YOUR OWN PATH/Tree_Clonal_evolution_%d.dat"%seed

    for node in Tumor.traverse(strategy="postorder"):
        node.write(features=[], outfile=tree_fn, format=2, is_leaf_fn=None, format_root_node=True)

    # CHECK THE PATH
    temp_fn="YOUR OWN PATH/Clonal_evolution_%d.dat"%seed
    temp_f=open(temp_fn,"w")
    temp_f.write("#mut_rate\t%f\tfitness\t%f\n"%(mut_rate,fitness))
    
    for node in Tumor.traverse("postorder"):
        temp_f.write("%d\t%s\n"%(int(node.cell),node.mutation))
    temp_f.close()







