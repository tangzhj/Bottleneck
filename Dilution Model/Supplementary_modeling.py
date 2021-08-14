# Here, we simulate the Dilution Model of mitochondrial bottleneck
# Number of generations is t.
# Td: Dilution generations 
# alpha: mtDNA replication rate
# mu: mutation rate
# Date: 06/07/2021
#
# ! /usr/local/bin/python
# In[1]
import sys,os,math,random
import numpy as np
from collections import Counter
import matplotlib.pyplot as plt
import pandas as pd
import random
import time
import multiprocessing
from multiprocessing import Pool

class mtDNApop():
    def __init__(self):
        self.muts = []

def initPop(N0, initFreq = 0):
    """
    initiate the population,
    N0 represents the initial mtDNA population size, eg.500
    """
    muts_init = int(round(initFreq*N0))
    muts_asn = Counter(np.random.choice(N0,muts_init, replace = False))
    pop_init = mtDNApop()
    for i in range(N0):
        pop_init.muts.append(str(muts_asn[i]))
    return(pop_init)

def imdtPopWF(pop, n):
    """
    intermediate population with initial mutations
    pop represents the mtDNA population from the previous generation
    n refers to the new generation population size
    """
    new_choice = np.random.choice(len(pop.muts),n)
    muts_pre = [pop.muts[i] for i in new_choice]
    pop_new = mtDNApop()
    pop_new.muts = muts_pre
    return(pop_new)

def imdtPopWFmut(pop, n, mu):
    """
    intermediate population with novel mutations
    """
    # pick n offsprings from previous population
    new_choice = np.random.choice(len(pop.muts),n)
    muts_pre = [pop.muts[i] for i in new_choice]
    # mutations of offspring generation
    muts_new = np.random.poisson(mu*n)
    muts_asn = Counter(np.random.choice(n,muts_new, replace = True))
    muts_id = range(1,muts_new+1)
    pop_new = mtDNApop()
    my_id = 0
    for i in range(n):
        if muts_asn[i] >0:
            l = range(my_id+1, my_id + 1 + muts_asn[i])
            l = '_'.join(map(str, l))
            my_id = my_id + muts_asn[i]
        else:
            l = 0
        pop_new.muts.append(','.join((str(muts_pre[i]),str(l))))
    return(pop_new)

def tPopFreq(N0, initFreqs, mu, t, dilutionT, alpha, N1): # tOnly = True indicates output all generations
    """
    simulate the Nth generation
    """
    # pick frequency once a time, total k times
    df = pd.DataFrame()
    for myF in range(len(initFreqs)):
        pop_new = initPop(N0, initFreqs[myF])
        for i in range(t):
            if dilutionT is not None and i in dilutionT:
                n = int(round(len(pop_new.muts)*alpha))
                pop_new = imdtPopWF(pop_new, n)
            else:
                pop_new = imdtPopWF(pop_new, N1)
        df[str(myF)] = pop_new.muts
    # add mutations with mutation rate mu
    pop_new = initPop(N0)
    for i in range(t):
        if dilutionT is not None and i in dilutionT:
            n = int(round(len(pop_new.muts)*alpha))
            pop_new = imdtPopWFmut(pop_new, n, mu)
        else:
            pop_new = imdtPopWFmut(pop_new, N1, mu)
    mut = pd.DataFrame([x.split(',') for x in pop_new.muts])
    del mut[0]
    # concatenate existing mutations and new arising mutations
    s = pd.concat([df, mut],axis=1)
    p = []
    for myP in range(s.shape[1]):
        temp = s.iloc[:, myP]
        temp = [x.split('_') for x in temp]
        temp = [item for sublist in temp for item in sublist]
        temp = Counter(temp)
        for myMut in temp:
            if myMut != '0':
                p.append(temp[myMut]/s.shape[0])
    p = [num for num in p if num > 0.01]
    return(p)

# %% run simulation
N0 = 500
myT = 20
N1s = [500,500,300]
mu = 3e-7
numCells = [786, 2763, 1064] ## B, T, NK

def runSimul(numSimul):
    """
    process the numSimul file of simulation
    """
    print('Run task %s (%s)...' % (numSimul, os.getpid()))
    N0 = 500
    numCell = 786
    mu0 = 1e-7
    mu = mu0*16569*numCell
    N1 = 500
    np.random.seed(numSimul)
    freqTable = pd.read_csv('LMPP_mutation_freqencies_0418.txt',delimiter='\t')
    cellsID = freqTable.cellID.unique()
    i=0
    with open('Dilution_Simul_'+str(N1)+'mtDNAcopy'+str(numSimul)+'.txt','w') as f:
        while i < 10000:
            cells = np.random.choice(cellsID, numCell, replace = True)
            initFreqs = freqTable[freqTable.cellID.isin(cells)].freq.tolist()
            Td = random.randint(1,30)
            t = random.randint(10,40)
            dilutionT = range(Td)
            alpha = random.uniform(0,1)
            if t > Td and alpha**Td > 10/500:
                p = tPopFreq(N0, initFreqs, mu, t, dilutionT, alpha, N1)
                f.write(str(numCell)+'\t'+str(alpha)+'\t'+str(mu0)+'\t'+str(Td) + '\t' + str(t) + '\t' + str(round(500*alpha**Td)) + '\t' + str(p) + '\n')
                i +=1
    f.close()

if __name__=='__main__':
    print('Parent process %s.' % os.getpid())
    p = Pool(100)
    for i in range(100):
        p.apply_async(runSimul, args = (i,))
    print('Waiting for all subprocesses done...')
    p.close()
    p.join()
    print('All subprocesses done.')
