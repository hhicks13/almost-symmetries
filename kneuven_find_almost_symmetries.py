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
from operator import attrgetter
from networkx.algorithms import bipartite
import timeit


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
        pass
    lN,rN = sortByLength(Nu,Nv)
    dim = [len(lN),len(rN)]
    return lN,rN,dim
    
#left neighborhood is smaller
def buildCostMatrix2(i,j,Ni,Nj,dim,_kA):
    #line 9 in BCM
    Ld = abs(dim[1]-dim[0])+_kA
    Rd = _kA
    Pair = namedtuple('Pair',('isNeighbor','name','parent'))
    LeftBp = [ Pair(1,idv,i) if idi == 0 else Pair(0,idv,i) for idi, vertices in enumerate([Ni,[i for i in range(Ld) ]]) for idv in vertices] 
    RightBp =[ Pair(1,idv,j) if idi == 0 else Pair(0,idv,j) for idi, vertices in enumerate([Nj,[i for i in range(Rd) ]]) for idv in vertices]
    #print(LeftBp)
    #print(RightBp)
    return LeftBp,RightBp

#u in left set and v in right#
#u and v are named tuples
def buildCostMatrix3(i,j,u,v,_GA,_PA,_EFA):
    #3 cases:
    #v2v
    _u = u.name
    _v = v.name
    if u.isNeighbor == 1 and v.isNeighbor == 1:
        if (_u,_v) in _PA.edges or (_v,_u) in _PA.edges:
            # I dont believe these two cases are possible
            if  _u not in list(_PA.adj[i]) and  _GA.degree(_u) > _GA.degree(_v):
                return -2*abs(_GA.degree(_u) - _GA.degree(_v))
            elif _v not in list(_PA.adj[j]) and _GA.degree(_v) > _GA.degree(_u):
                return -2*abs(_GA.degree(_u) - _GA.degree(_v))
            else:
                return -abs(_GA.degree(_u) - _GA.degree(_v))
        else:
            #proven impossible
            return -999
    #d2d
    elif u.isNeighbor == 0 and v.isNeighbor == 0:
        return 0
    #d2v/v2d
    else:
        if (u.parent,_u) in list(_EFA.edges) or (v.parent,_v) in list(_EFA.edges):
            #fixed
            return -999
        else:
            return -2


def buildCostMatrix4(i,j,_GA,_PA,_EFA,_kA):
    #print(len(_GA.adj[i]))
    #print(len(_GA.adj[j]))
    Ni,Nj,dim = buildCostMatrix1(_GA,i,j)
    LeftBp,RightBp = buildCostMatrix2(i,j,Ni,Nj,dim,_kA)
    # draw edge & weight, line 10
    CostMatrix = nx.Graph()
    for e in it.product(LeftBp,RightBp):
        if e[0].isNeighbor == 0 and e[1].isNeighbor == 0 and e[0].name != e[1].name:
            pass
        else:
            CostMatrix.add_node(e[0],bipartite=0)
            CostMatrix.add_node(e[1],bipartite=1)
            CostMatrix.add_edge(e[0],e[1],weight=buildCostMatrix3(i,j,e[0],e[1],_GA,_PA,_EFA))
    #for u,v,weight in sorted(CostMatrix.edges.data("weight"),key=itemgetter(2,0,1)):print((u.isNeighbor,v.isNeighbor),u.name,v.name,weight)
    M = nx.max_weight_matching(CostMatrix,maxcardinality=True,weight="weight")
    assert(nx.is_perfect_matching(CostMatrix,M))
    return CostMatrix

def HungarianSolve(CostMatrix):
    cost = 0
    deleteEdges = nx.max_weight_matching(CostMatrix,maxcardinality=True)
    [cost:= cost+(-1*CostMatrix[l][r]['weight']) for l,r in deleteEdges if CostMatrix[l][r]['weight'] < 0]
    assert(cost >= 0)
    return cost, deleteEdges


def RefineByMatching(_GA,_PA,_EFA,_kA):
    print("RefineByMatching: y/n?")
    print(_kA)
    if input() == 'y':
        pass
    else:
        pass
    changed = False
    
    print(_kA)
    edgeUse = {(e[0],e[1]):0 for e in list(_GA.edges)}
    #
    PA = _PA
    for edge in list(PA.edges):
        i = edge[0]
        j = edge[1]
        CostMatrix = buildCostMatrix4(i,j,_GA,_PA,_EFA,_kA)
        cost,deleteEdges = HungarianSolve(CostMatrix)
        assert(nx.is_perfect_matching(CostMatrix,deleteEdges))
        print("COST:"," "*20,cost)
        print("2*kA"," "*20,2*_kA)
        print("DELETE EDGES:",deleteEdges)
        if input() == 'y':
            pass
        else:
            pass
        #
        if cost > 2*_kA:
            print(cost)
            print(edge)
            E_PA.remove((j,i))
            PA.remove_edge(i,j)
            print(E_PA,"!"*100)
        else:
            for u,v in deleteEdges:
                if u.name == v.name:
                    continue
                e3 = (u.name,v.name)
                if e3 in edgeUse.keys():
                    edgeUse[e3]+=1
                else:
                    print("SOMETHING IS WRONG")
                    if input() == 'y':
                        pass
                    else:
                        pass
                    
        s = list(sorted(edgeUse.items(), key=itemgetter(1), reverse=True))
        for edge in s[:25]:print(edge[1],end=", ")
        print("Done!")
        print("-"*20,"Maximum edgeUsage",max(edgeUse.items(),key=itemgetter(1))[1])
        return edgeUse,PA



# call count_automorphisms_vf2 on G.
def ComputeAutomorphisms(_GA):
    G = ig.Graph.from_networkx(_GA)
    autNum = G.count_automorphisms_vf2()
    return autNum

# Need set E(PA) exposed
# _GA must be g
def DegreeDiffElim(_GA,_PA,_kA):
    E_PA = list(_PA.edges)
    PA = _PA
    for edge in E_PA[:]:
        i = edge[0]
        j = edge[1]
        if abs(_GA.degree(i)-_GA.degree(j)) > _kA:
            E_PA.remove(edge)
            PA.remove_edge(i,j)
    
    print("-"*20,"refinement",PA)
    return PA

# dEFA(i) to be the number of fixed edges incident to i
#
def dEFA(i,_EFA):
    adjacent_fixed = _EFA.adj[i]
    return len(adjacent_fixed)

def FixedDefElim(_GA,_PA,_EFA):
    PA = _PA
    E_PA = list(_PA.edges)
    for edge in E_PA[:]:
        i = edge[0]
        j = edge[1]
        if dEFA(i,_EFA) > _GA.degree(j):
            PA.remove_edge(i,j)
    return PA


# G = tograph(_PA)
# nx.maximal_independent_set(G)
# return size of preceding
def GreedyIndependentSetSize(_PA):
    independent_set =  nx.maximal_independent_set(_PA)
    print(len(independent_set),"$$$$$$$$$")
    return len(independent_set)


#
# edgeUse is a dictionary
#

# in _E(GA) \ _EFA
# if neighborhood of i or of j are not empty (off diagonal)
# return the first encountered
# else prune
def FindBranchEdge(_GA,_PA,_EFA,_edgeUse):
    items = _edgeUse.items()
    argmax = max(items, key=itemgetter(1))[0]
    print("ARGMAX",argmax)
    if input() == 'y':
        pass
    else:
        pass
    
    if _edgeUse[argmax] > 0:
        return argmax
    fixed = set(list(_EFA.edges))
    remaining_edges = [sorted(edge,key=itemgetter(1)) for edge in list(_GA.edges) if edge not in fixed]
    print("Remaining_edges of GA",remaining_edges)
    for edge in remaining_edges:
        if i == j:
            continue
        else:
            i = edge[0]
            j = edge[1]
            _NPAi = list(_PA.adj[i])
            _NPAj = list(_PA.adj[j])
            if len(_NPAi) > 0 or len(_NPAj) > 0:
                return (i,j)
    return -1




def main():
    # used to distinguish dummy variables from vertices.
        
    n = int(sys.argv[1])
    p = float(sys.argv[2])
    G = nx.erdos_renyi_graph(n,p)

    #
    Node = namedtuple('Node',('PA','EDA','EFA','delChild'))
    spacer = " "*4
    T = 0
    t = 0
    #

    print("Begin FindAlmostSymmetry")

    #
    print("Initialize Root Node")
    k = int(sys.argv[3])
    P_R = nx.complete_graph(n)
    EDR = nx.Graph();EDR.add_nodes_from(list(G.nodes))
    EFR = nx.Graph();EFR.add_nodes_from(list(G.nodes))
    delChild = True
    #
    print("Incumbents")
    incumbentValue = n #|G|
    incumbentSolution = 0 #EDA
    # input to nodestack is iterable
    print("Append root to stack")
    NodeStack = deque()
    NodeStack.append(Node(P_R,EDR,EFR,delChild))
    print("enter loop")
    while True:
        T +=1
        try:
            #
            print(T,": Unpack Node data")
            A = NodeStack.pop();
            
            pa = nx.Graph(A.PA)
            print("EDGES OF PA",len(pa.edges))
            eda = nx.Graph(A.EDA);
            print("EDGES OF EDA",len(eda.edges))
            efa = nx.Graph(A.EFA);
            print("EDGES OF EFA",len(efa.edges))
            delChild = A.delChild;
            print("delChild",delChild)
            #
            print("constructing GA")
            GA = G
            for e in list(eda.edges):
                if e[0] != e[1]:
                    GA.remove_edge(e[0],e[1])
            
            kA = k - len(eda.edges)
            print("EDGES of GA",GA.edges)
            print("value of kA",kA)
            #
            print("Fix:",efa.edges,"Del",eda.edges)
            if delChild:
                print(spacer,"Node signal: Delete Edges, compute aut")
                autNum = ComputeAutomorphisms(GA)
                print("Autnum",autNum)
                print("incumbent",incumbentValue)
                if autNum < incumbentValue:
                    print(spacer,spacer,"update incumbents")
                    incumbentValue = autNum
                    incumbentSolution = eda
                if len(efa.edges) == len(GA.edges) or kA == 0:
                    print(spacer,spacer,"prune: no more edges / budget depleted")
                    break
                pa = DegreeDiffElim(GA,pa,kA) #update PA
                print(pa.edges)
                print(spacer,"updated edges of P_A:",len(pa.edges)," members")
            #
            else:
                print(spacer,"Node Signal: Fix Edges, refine P_A")
                if len(efa.edges) == len(GA.edges):
                    print(spacer,spacer,"prune: no more edges")
                    continue
                print(spacer,"update pa, FixedDegreeElim")
                pa = FixedDefElim(GA,pa,efa)
            #
            print(T,".0 lower bound via greedyIndependentSet")
            lowerBound = GreedyIndependentSetSize(pa)
            #
            if lowerBound >= incumbentValue:
                continue
            print(T,".1 lower bound loop")
            print(spacer,T,".1.1: While P_A changed by RefinebyMatching")
            edgeUse,_pa = RefineByMatching(GA,pa,efa,kA)
            while pa.edges != _pa.edges:
                t +=1
                lowerBound = GreedyIndependentSetSize(_pa)
                print(spacer,spacer,T,".",t,". computing bound")
                if lowerBound >= incumbentValue:
                    print(spacer,spacer,T,".",t," lowerbound >= incumbentValue: prune node")
                    break
                else:
                    print(T,".2 new lower bound!")
                    edgeUse,pa = RefineByMatching(GA,pa,efa,kA)
            print(T,".3 FindBranchEdge")
            branchEdge = FindBranchEdge(GA,pa,efa,edgeUse)
            print(spacer,T,".3.1 bounding is concluded")
            # 
            if branchEdge == -1:
                print(spacer,T,".3.2 branchEdge signal: prune")
                continue
            #
            print(T,".4 create children")
            efa_branch = efa.add_edge(branchEdge[0],branchEdge[1]);
            efa_branch = efa.add_edge(branchEdge[1],branchEdge[0])#add new edges for fixation
            eda_branch = eda.add_edge(branchEdge[0],branchEdge[1]);
            eda_branch = eda.add_edge(branchEdge[1],branchEdge[0])#add new edges for deletion
            print(spacer,T,".4.1: append fChild")
            NodeStack.appendleft(Node(_pa,eda,efa_branch,False))
            print(spacer,T,".4.2: append dChild")
            NodeStack.appendleft(Node(_pa,eda_branch,efa,True))            
        #
        except IndexError:
            print(T,".5: end of stack")
            break
        
    print("incumbentValue",incumbentValue)
    print("incumbentSolution",incumbentSolution)
    

   
            
    
    
if __name__ == "__main__":
    main()
#
#
#
# print neuron-cosets based on memo2
