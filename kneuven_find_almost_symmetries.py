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
import pynauty

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

def assignweight(g,u,v):
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

def buildCostMatrix(num2namepoint,Nu,Nv,dim):
    LeftBp = [ (0,num2namepoint[idv]) if idi == 0 else (0,idv) for idi, vertices in enumerate([Nu,[i for i in range(dim[1]) ]]) for idv in vertices] 
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
    Rho[Tau+1] = P1
    DELTAS = {}
    killed = {}
    eliminations_0 = 0
    _worstcase = []
    images = []
    t = 0
    while eliminations_0 < n*(n-1):
        epsilon = Tau- t
        # Eliminate if over budget
        Rho[epsilon] = [ (e[0],e[1],-1) if e[2]  > 2*epsilon else e for e in Rho[epsilon+1] ]

        # zero matrix
        DELTAS[epsilon] = []
        Delta = [[0]*n for _ in range(n)]

        for e in Rho[epsilon]:print(f'{Fore.BLACK}{Back.RED} {"x"} {Style.RESET_ALL}' if e[2]<0 else  f'{e[2]} {Style.RESET_ALL}',end=' ')
        print('\n')
    
        for e in Rho[epsilon]:
    
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
        eliminations_0 = len(killed.keys())
        _worstcase.append(eliminations_0)
        t+=1
    return Rho,DELTAS,_worstcase,images

# Ostrowski & Kneuven
#def FindAlmostSymmetry(_n,_PR,_G,_k):
#    G = _G # set of tuples
#    PR = _PR
#    eDR = None # empty set
#    eFR = None # empty set
#    delChild = True
#    incumbentValue = n
#    incumbentSolution = None
#    NodeStack = []
#    NodeStack.append((PR,eDR,eFR,delChild))
#    while len(NodeStack) != 0:
#        node = NodeStack[-1]
#        GA = G NOT eDA
#        kA = _k - len(eDA)
#        if node[3]:
#            orbitNum =     
#    return time
# Ostrowski & Hicks
#def Algorithm1():
#    return time

#def matrix_sum(base,deltas):
#    for d in deltas:
#        base += d;
#    return base

#def save_images(images):
    #fig, axes = plt.subplots(10,10, figsize=(8,8))
    #for i,ax in enumerate(axes.flat):
    #    ax.imshow(images[i])
    #plt.subplots(figsize = (n,n)) # width x height
    #sns.heatmap(images[3]) # row, column, position
    #ax2 = fig.add_subplot(3, 3, 2)
    #ax3 = fig.add_subplot(3, 3, 3)
    #ax4 = fig.add_subplot(3, 3, 4)
    #ax5 = fig.add_subplot(3, 3, 5)
    #plt.show()

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

    #
    Tau = int(sys.argv[3])
    P0 = []
    Rho = {}
    #
    print("Generate P")
    P0 = generate_P0(n,names,g)
    print("P generated on ",n,"nodes")
    print("DegreeDifElim")
    P1 = [ (e[0],e[1],0) if abs(g.degree(e[0]) - g.degree(e[1])) > Tau else e for e in P0 ]

    Rho,ELIMS,_worstcase,images = elim2_preprocess(n,Tau,P1)
    INITS = {}
    for ctr in range(Tau):
        epsilon = Tau - ctr
        _matrix = [[0]*n for _ in range(n)] # n x n table
        for (y,x,score) in Rho[epsilon]:_matrix[y][x] = score
        INITS[epsilon] = np.asarray(_matrix)
    print(len(images))
    print(images[3])
    print(_worstcase)
            
    
    
if __name__ == "__main__":
    main()
#
#
#
# print neuron-cosets based on memo2
