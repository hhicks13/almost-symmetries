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
import igraph as ig
import re
import math
import seaborn as sns
from collections import deque
from collections import namedtuple
from numpy.random import permutation
from operator import itemgetter


#### previous

def generate_P0(n,names,g):
    P0 = []
    for u,v in it.product(range(n),range(n)):
        Nu,Nv,dim = buildCostMatrix1(n,g,u,v)
        #
        B,M  = hungarianSolve(g,Nu,Nv,dim)
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
########################################################################################################### Kneuven

def sortByLength(set1,set2):
    if len(set1) == len(set2):
        return set1, set2
    elif len(set1) < len(set2):
        return set1, set2
    else:
        return set2,set1


#    grid0 = [range(n),names]
#    num2namepoint = {j[0]:j for j in it.zip_longest(*grid0)}


def buildCostMatrix1(_GA,u,v):
    g = _GA
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
    
#left neighborhood is smaller
def buildCostMatrix2(Ni,Nj,dim,kA):
    #line 9 in BCM
    Ld = dim[1]-dim[0]+kA
    Rd = kA
    Pair = namedtuple('Pair',('isRight','isNeighbor','name'))
    LeftBp = [ Pair(0,1,idv) if idi == 0 else Pair(0,0,idv) for idi, vertices in enumerate([Ni,[i for i in range(Ld) ]]) for idv in vertices] 
    RightBp =[ Pair(1,1,idv) if idi == 0 else Pair(1,0,idv) for idi, vertices in enumerate([Nj,[i for i in range(Rd) ]]) for idv in vertices]
    return LeftBp,RightBp

#u in left set and v in right#
#u and v are named tuples
def buildCostMatrix3(i,j,u,v,_GA,_PA,_efa):
    #3 cases:
    #v2v
    _u = u.name
    _v = v.name
    if u.isNeighbor == 1 and v.isNeighbor == 1:
        if (_u,_v) in _PA.edges or (_v,_u) in _PA.edges:
            # I dont believe these two cases are possible
            if (not _u in nx.neighbors(_PA,i)) and  _GA.degree(_u) > _GA.degree(_v):
                return -2*abs(_GA.degree(_u) - _GA.degree(_v))
            elif (not _v in nx.neighbors(_PA,j)) and _GA.degree(_v) > _GA.degree(_u):
                return -2*abs(_GA.degree(_u) - _GA.degree(_v))
            else:
                return -abs(_GA.degree(_u) - _GA.degree(_v))
        else:
            return -999
    #d2d
    elif u.isNeighbor == 0 and v.isNeighbor == 0:
        return 0
    #d2v/v2d
    else:
        if (i,_u) in _efa.edges or (j,v) in _efa.edges:
            return -999
        else:
            return -2


def buildCostMatrix4(i,j,_GA,_PA,_efa,_kA):
    Ni,Nj,dim = buildCostMatrix1(_GA,i,j)
    LeftBp,RightBp = buildCostMatrix2(Ni,Nj,dim,_kA)
    # draw edge & weight, line 10
    weights = { (e[0].name,e[1].name):buildCostMatrix3(i,j,e[0],e[1],_GA,_PA,efa)) for e in it.product(LeftBP,RightBp) }
    weighted_edges = [ (e[0],e[1],weight[(e[0].name,e[1].name)]) for e in it.product(LeftBp,RightBp) ] 
    CostMatrix = nx.Graph()
    CostMatrix.add_nodes_from(LeftBp,bipartite=0)
    CostMatrix.add_nodes_from(RightBp,bipartite=1)
    CostMatrix.add_weighted_edges_from(weighted_edges)
    return CostMatrix,weights

def HungarianSolve(CostMatrix,weights):
    cost = 0
    deleteEdges = nx.max_weight_matching(CostMatrix,maxCardinality=True)
    for e in deleteEdges:cost+=weights[(e[0].name,e[1].name)]
    return cost,deleteEdges

def RefineByMatching(_GA,_PA,_EFA,kA):
    n = len(_GA.nodes)
    g = nx.Graph(_GA)
    edgeUse = {}
    changed = False
    #
    for e in _PA.edges:edgeUse[e] = 0
    #
    E_P_A = _PA.edges
    for edge in E_P_A[:]:
        cost = 0
        i = edge[0]
        j = edge[1]
        CostMatrix,weights = BuildCostMatrix4(_GA,i,j,_efa,_PA,kA)
        edgeUsed = HungarianSolve(CostMatrix)
        print(i,j," perfect matching: ",nx.is_perfect_matching(B,deleteEdges))
        for e in deleteEdges:cost+=extractweight(g,e[0][1],e[1][1])
        print(cost)
        print(deleteEdges)
        #
        if cost > 2*kA:
            EPA.remove[edge]
            changed = True
        else:
            for e in deleteEdges:edgeUse[e]+=1
        #
        pa = None # new graph based on E_P_A
            
    return changed,edgeUse,pa


# call count_automorphisms_vf2 on G.
def ComputeAutomorphisms(_GA):
    return autNum

# Need set E(PA) exposed
# _GA must be g
def DegreeDiffElim(_GA,_PA,_kA):
    return PA_refined

# dEFA(i) to be the number of fixed edges incident to i
#
def dEFA(i):
    return adjacent_fixed

def FixedDefElim(_GA,_PA,_EFA):
    return PA_refined


# G = tograph(_PA)
# nx.maximal_independent_set(G)
# return size of preceding
def GreedyIndependentSetSize(_PA):
    return lower_bound

def removeFromEPA(_PA,i,j):
    return _PA
#
# edgeUse is a dictionary
#

# in _E(GA) \ _EFA
# if neighborhood of i or of j are not empty (off diagonal)
# return the first encountered
# else prune
def FindBranchEdge(_GA,_PA,_EFA,_edgeUse):
    return branchEdge




def main():
    # used to distinguish dummy variables from vertices.
    names = []
    with open('data/vnames.txt') as f:
        lines = f.readlines()
        for line in lines:
            for num in re.findall("\d+",line):names.append(int(num))
    print("maximum vertices",len(names))
    
    n = int(sys.argv[1])
    p = float(sys.argv[2])
    G = nx.erdos_renyi_graph(n,p)

    #
    MaxBudget = int(sys.argv[3])
    Node = namedtuple('Node',('PA','EDA','EFA','delChild'))
    spacer = " "*4
    T = 0
    #

    print("Begin FindAlmostSymmetry")

    #
    print("Initialize Root Node")
    k = MaxBudget
    P_R = nx.complete_graph(n)
    EDR = nx.erdos_renyi_graph(n,0)
    EFR = nx.erdos_renyi_graph(n,0)
    delChild = True
    #
    print("Incumbents")
    incumbentValue = n #|G|
    incumbentSolution = nx.erdos_renyi_graph(n,0) #EDA
    # input to nodestack is iterable
    print("Append root to stack")
    NodeStack = deque()
    NodeStack.append(Node(P_R,EDR,EFR,delChild))
    print("Enter loop.")
    while True:
        T +=1
        try:
            #
            print(T,": Unpack Node data")
            A = NodeStack.pop()
            pa = nx.Graph(A.PA)
            eda = nx.Graph(A.EDA)
            efa = nx.Graph(A.EFA)
            delchild = A.delChild
            #
            GA = None #remove EDA from G !!
            kA = k - len(eda.edges)
            #
            print("Fix:",efa.edges,"Del",eda.edges)
            if delChild:
                print(spacer,"Node signal: Delete Edges, compute aut")
                autNum = ComputeAutomorphisms(GA)
                if autNum < incumbentValue:
                    print(spacer,spacer,"update incumbents")
                    incumbentValue = autNum
                    incumbentSolution = nx.Graph(eda)
                if len(efa.edges) == len(GA.edges) or kA == 0:
                    print(spacer,spacer,"prune: no more edges / budget depleted")
                    continue
                pa = DegreeDiffElim(GA,pa,kA) #update PA
                print(spacer,"updated edges of P_A",pa.edges)
            #
            else:
                print(spacer,"Node Signal: Fix Edges, refine P_A")
                if len(efa.edges) == len(GA.edges):
                    print(spacer,spacer,"prune: no more edges")
                    continue
                print(spacer,"update pa, FixedDegreeElim")
                pa = FixedDegreeElim(GA,pa,efa)
            #
            print(T,".0 lower bound via greedyIndependentSet")
            lowerBound = GreedyIndependentSetSize(pa)
            #
            if lowerBound >= incumbentValue:
                continue
            print(T,".1 lower bound loop")
            _prune_rbm = False
            t = 0
            print(spacer,T,".1.1: While P_A changed by RefinebyMatching")
            while True:
                changed,edgeUse,_pa = RefineByMatching(GA,pa,efa,kA)
                if not changed:
                    print(spacer,spacer,T,".1.1.0 current pa unchanged by refinement")
                    break
                else:
                    t +=1
                    pa = nx.Graph(_pa)
                    lowerBound = GreedyIndependentEdgeSet(pa)
                    print(spacer,spacer,T,".1.1.",t,"computing bound")
                    if lowerBound >= incumbentValue:
                        print(spacer,spacer,T,".1.1.",t," lowerbound >= incumbentValue: prune node")
                        _prune_rbm = True
                        break
            if _prune_rbm:
                print(spacer,T,".1.2 node pruned. lowerbound >= incumbentValue")
                continue
            else:print(T,".2 new lower bound!")
            
            print(T,".3 FindBranchEdge")
            branchEdge = FindBranchEdge(GA,pa,efa,edgeUse)
            print(spacer,T,".3.1 bounding is concluded")
            # 
            if branchEdge == -1:
                print(spacer,T,".3.2 branchEdge signal: prune")
                continue
            #
            print(T,".4 create children")
            _efa_branch = None #add new edges for fixation
            _eda_branch = None #add new edges for deletion
            print(spacer,T,".4.1: append fChild")
            NodeStack.append(Node(pa,eda,_efa_branch,False))
            print(spacer,T,".4.2: append dChild")
            NodeStack.append(Node(pa,_eda_branch,efa,True))            
        #
        except IndexError:
            print(T,".5: end of stack")
            break
        
    print("incumbentValue",incumbentValue)
    print("incumbentSolution",incumbentSolution)
        
    
    
    
    # _EFR, _EDR graphs.
    print("P generated on ",n,"nodes")
    print("DegreeDifElim")
    print("FixedDegreeElim")
    print("ComputeAutomorphisms")
    print("RefineByMatching")
    print("GreedyInependentSetSize")
    print("FindBranchEdge")
    

   
            
    
    
if __name__ == "__main__":
    main()
#
#
#
# print neuron-cosets based on memo2
