#!/usr/bin/env python
from array import*
import sys
import time
from itertools import chain
import networkx as nx
import numpy as np
import random
import itertools as it
import networkx as nx
import re
import math
import fractions
import seaborn as sns

from numpy.random import permutation
from colorama import Fore
from colorama import Style
from colorama import Back
from operator import itemgetter

from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection

from sklearn import manifold
from sklearn.metrics import euclidean_distances
from sklearn.decomposition import PCA

from sklearn.manifold import MDS
from mpl_toolkits import mplot3d


########################################################### Graph Handling

def sortByLength(set1,set2):
    if len(set1) == len(set2):
        return set1, set2
    elif len(set1) < len(set2):
        return set1, set2
    else:
        return set2,set1

def getNameSpace(names,n):
    grid0 = [range(n),names]
    num2namepoint = {j[0]:j for j in it.zip_longest(*grid0)}

    if len(names) == n:
        return num2namepoint
    else: return None

def pair2neighborhood(g,u,v):
    Nu = [e for e in nx.neighbors(g,u)]
    Nv = [f for f in nx.neighbors(g,v)]
    try:
        Nu.remove(v)
        Nv.remove(u)
    except(ValueError,TypeError):
        print("")
    lN,rN = sortByLength(Nu,Nv)
    dim = [len(lN),len(rN)]
    return lN,rN,dim

def genMunkresBasis(names,n,g,u,v):
    num2namepoint = getNameSpace(names[:n],n)
    Nu,Nv,dim = pair2neighborhood(g,u,v)
    return num2namepoint,Nu,Nv,dim

# # # # # # graph weights # # # # # # # # # # # # #

def assignweight(g,u,v):
    # negated here because of networkx requirement
    if type(u) is tuple and type(v) is tuple:
        return -abs(g.degree(u[0]) - g.degree(v[0]))    
    elif type(u) is int and type(v) is int:
        return 0
    else:
        return -1
    
def extractweight(g,u,v):
    if type(u) is tuple and type(v) is tuple:
        return abs(g.degree(u[0]) - g.degree(v[0]))
    elif type(u) is int and type(v) is int:
        return 0
    else:
        return 1
############################################################################ COMPUTATION
# capitalized

# order of ops top2down # order of names L2R # names ordered by computational dependency

def buildCostMatrix(num2namepoint,Nu,Nv,dim):
    LeftBp = [ (0,num2namepoint[idv]) if idi == 0 else (0,idv) for idi, vertices in enumerate([Nu,[i for i in range(dim[1]) ]]) for idv in vertices] # counter intuitive, but by doing the same operation twice you save more time
    RightBp =[ (1,num2namepoint[idv]) if idi == 0 else (1,idv) for idi, vertices in enumerate([Nv,[i for i in range(dim[0]) ]]) for idv in vertices]
    #print(len(LeftBp))
    #print(len(RightBp))
    return LeftBp,RightBp

def hungarianSolve(g,num2namepoint,Nu,Nv,dim):
    LeftBp,RightBp = buildCostMatrix(num2namepoint,Nu,Nv,dim)
    weighted_edges = [ (e[0],e[1],assignweight(g,e[0][1],e[1][1])) for e in it.product(LeftBp,RightBp) ] #O(n^2)
    B = nx.Graph()
    B.add_nodes_from(LeftBp,bipartite=0)
    B.add_nodes_from(RightBp,bipartite=1)
    B.add_weighted_edges_from(weighted_edges)
    M  = nx.max_weight_matching(B,maxcardinality=True)
    return B,M

def generate_P0(n,names,g):
    P0 = []
    for u,v in it.product(range(n),range(n)):
        num2namepoint,Nu,Nv,dim = genMunkresBasis(names,n,g,u,v)
        #
        B,M  = hungarianSolve(g,num2namepoint,Nu,Nv,dim)
        print(u,v," perfect matching: ",nx.is_perfect_matching(B,M))
        linsum = extractweight(g,u,v)
        for m in M:linsum+=extractweight(g,m[0][1],m[1][1])
        print(dim)
        P0.append((u,v,linsum))
    return P0

def elim2_preprocess(n,Tau,P1):
    Rho = {}
    ############################################# Refinement Algorithm 2: JimElim #######
    Rho[Tau+1] = P1
    DELTAS = {}
    killed = {}
    eliminations_0 = 0
    
    
    _worstcase = []
    images = []
    
    t = 0
    # casualties may not exceed dimensions of P
    # however since casualties increase exponentially, this should not be a huge computational burden
    
    while eliminations_0 < n*(n-1):
        epsilon = Tau- t
        # Eliminate if over budget
        Rho[epsilon] = [ (e[0],e[1],-1) if e[2]  > 2*epsilon else e for e in Rho[epsilon+1] ]

        # zero matrix
        DELTAS[epsilon] = []
        Delta = [[0]*n for _ in range(n)]

        for e in Rho[epsilon]:print(f'{Fore.BLACK}{Back.RED} {"x"} {Style.RESET_ALL}' if e[2]<0 else  f'{e[2]} {Style.RESET_ALL}',end=' ')
        print('\n')
        #
        # add delta to Rho (while pairs memoized ) O(e*n)
        for e in Rho[epsilon]:
        ######################################################## compute delta sum
            if e[2] == -1 and e[0] != e[1]:
                killed[(e[0],e[1])] = True
                killed[(e[1],e[0])] = True
                deleted = killed.keys()
                for i in range(n):
                    L = Rho[epsilon][e[0]*n+i]
                    R = Rho[epsilon][e[1]*n+i]
                    #print(L)
                    #print(R)
                    #
                    # no redundancy
                    #L
                    if not (L[0],L[1]) in deleted:Delta[L[0]][L[1]] = Delta[L[0]][L[1]]+2
                    if not (L[1],L[0]) in deleted:Delta[L[1]][L[0]] = Delta[L[1]][L[0]]+2
                    #R
                    if not (R[0],R[1]) in deleted:Delta[R[0]][R[1]] = Delta[R[0]][R[1]]+2
                    if not (R[1],R[0]) in deleted:Delta[R[1]][R[0]] = Delta[R[1]][R[0]]+2
                    # remove diagonals and duplicates
                Delta[e[0]][e[0]] = 0
                Delta[e[1]][e[1]] = 0

                # print Delta
                images.append(np.asarray(Delta))
                DELTAS[epsilon].append(np.asarray(Delta))

        print('\n')
                
       #########################################################

       # Rho+delta = lowerbound
       # 1 , 2, 3, .... .        n
       # # # # L'# ##L # # # # # # R+L=n cost matrix
       #       #
       # v2d   #  v2v
       # # # # L- # #L # # # # # # due to hungariansolve there must be exactly 2 dummies in this region 
       #       #     ^ dummies here ^ d2v after
       # d2d   #  < # # # # # # #
       #       #    #    d2v before
       #       #  < #                 cost matrix realignment, creates exactly 1 more '2'.
       #       #    #
       #       R+   R                 |v2d| = L^2, |d2V| = R^2
       # 
        eliminations_0 = len(killed.keys())
        _worstcase.append(eliminations_0)
        #print(eliminations," Eliminations at budget: ",epsilon)
        #increment counter
        t+=1
    del Rho[Tau+1]
    return Rho,DELTAS,_worstcase,images

#
def matrix_sum(base,deltas):
    for d in deltas:
        base += d;
    return base


# eliminates overbudget entries and shares info with sew
#def reap(matrix,budget,alive,dead,elimcount):
#    _new_matrix = matrix
#    _just_died = []
#    _alive = alive
#    _dead = dead
#    _elim_ctr = elimcount
    #filter non-living entries
#    for e in _alive.keys():
#        u = e[0]
#        v = e[1]
#        if matrix[u][v] > budget:
#            _new_matrix[u][v] = -1
            #record "time of death"
#            _dead[(u,v)] = _elim_ctr
            #update death count/time
#            _elim_ctr += 1
            #add to list for next phase
#            _just_died.append((u,v)]
            #update living
#            del _alive[(u,v)]
#        else:
            #nobody died
#            _new_matrix[u][v] = matrix[u][v]
#    return _new_matrix,_just_died,_dead,_alive

# 
#def sew(n,new_matrix,just_died,alive):
    #
    # delta matrices here are the partial matrices that result from 
    # adding 2 iteratively to neighborhoods. each iteration is saved as a matrix, and when the matrices are summed together they are equivalent to a single pass. 
    # this algorithm is essentially unwinding the process as much as possible to show intermediate steps
    #
#    _deltas = []
#    _just_died = just_died
#    _alive = alive.keys()
#    # empty D
#    D = np.asarray([[0]*n for _ in range(n)])
#    # those who just died
#    for e in _just_died:
#        # neighborhoods with dead filtered out (checks matrix row i and column i and row j and column j
#        for i in range(n):
#            # filters out dead
#            if (e[0],i) in _alive:
#                D[e[0]][i] = D[e[0]][i]+2
#            # filters out dead
#            if (e[1],i) in _alive:
#                [e[1]][i] = D[e[1]][i]+2
#            # filters out dead
#            if (i,e[0]) in _alive:
#                D[i][e[0]] = D[i][e[0]]+2
#            # filters out dead
#            if (i,e[1]) in _alive:
#                D[i][e[1]] = D[i][e[1]]+2
#        _deltas.append(D)
#    return _deltas

#


def main():
    #
    names = []
    with open('data/vnames.txt') as f:
        lines = f.readlines()
        for line in lines:
            for num in re.findall("\d+",line):names.append(int(num))
    print("maximum vertices",len(names))
    n = int(sys.argv[1])
    p = float(sys.argv[2])
    g = nx.erdos_renyi_graph(n,p)
    # for any counter, ctr, epsilon = Tau - t
    #
    Tau = int(sys.argv[3])
    P0 = []
    Rho = {} # Rho is indexed by budget
    #
    print("Generate P")
    P0 = generate_P0(n,names,g)
    print("P generated on ",n,"nodes")
    ############################################# Refinement Algorithm 1: DegreeDifElim #
    P1 = [ (e[0],e[1],0) if abs(g.degree(e[0]) - g.degree(e[1])) > Tau else e for e in P0 ]
    #################################################################### Elim2_preprocess 
    Rho,ELIMS,_worstcase,images = elim2_preprocess(n,Tau,P1)
    # convert to Numpy matrix
    INITS = {}
    for ctr in range(Tau):
        epsilon = Tau - ctr
        _matrix = [[0]*n for _ in range(n)] # n x n table
        for (y,x,score) in Rho[epsilon]:_matrix[y][x] = score
        INITS[epsilon] = np.asarray(_matrix)

    
    print(len(images))
    print(images[3])
    print(_worstcase)
    #fig, axes = plt.subplots(10,10, figsize=(8,8)))
    #for i,ax in enumerate(axes.flat):
    #    ax.imshow(images[i])
    plt.subplots(figsize = (n,n)) # width x height
    sns.heatmap(images[3]) # row, column, position
    #ax2 = fig.add_subplot(3, 3, 2)
    #ax3 = fig.add_subplot(3, 3, 3)
    #ax4 = fig.add_subplot(3, 3, 4)
    #ax5 = fig.add_subplot(3, 3, 5)

# We use ax parameter to tell seaborn which subplot to use for this plot
    sns.color_palette("rocket", as_cmap=True)
    #sns.heatmap(data=images[3], ax=ax1, cmap = "rocket", square=True, cbar_kws={'shrink': .3}, annot=True, annot_kws={'fontsize': 12})
    #sns.heatmap(data=subset2.corr(), ax=ax2, cmap = cmap, square=True, cbar_kws={'shrink': .3}, annot=True, annot_kws={'fontsize': 12})
    #sns.heatmap(data=subset3.corr(), ax=ax3, cmap = cmap, square=True, cbar_kws={'shrink': .3}, annot=True, annot_kws={'fontsize': 12})
    #sns.heatmap(data=subset4.corr(), ax=ax4, cmap = cmap, square=True, cbar_kws={'shrink': .3}, annot=True, annot_kws={'fontsize': 12})
    #sns.heatmap(data=subset5.corr(), ax=ax5, cmap = cmap, square=True, cbar_kws={'shrink': .3}, annot=True, annot_kws={'fontsize': 12})
    plt.show()
        
    
    
if __name__ == "__main__":
    main()
#
#
#
# print neuron-cosets based on memo2
