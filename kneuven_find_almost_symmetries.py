#!/usr/bin/env python
from array import*
import sys
import os
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
import collections
import multiprocessing
from collections import deque
from collections import namedtuple
from numpy.random import permutation
from operator import itemgetter
from operator import attrgetter
from networkx.algorithms import bipartite
import logging
from table_logger import TableLogger
from datetime import datetime

_PA = nx.Graph()
prune = False
k = int(sys.argv[3])
counter = 0
changed = True
lowerBound = 0
incumbentValue = 0
incumbentSolution = nx.Graph()
tbl = TableLogger(columns='counter,actual_time,kA,num edges in PA, num edges in GA')

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
        # Eliminate if over budget k < max(k)
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
    Nu = [e for e in nx.neighbors(_GA,u)]
    Nv = [f for f in nx.neighbors(_GA,v)]
    try:
        Nu.remove(v)
        Nu.remove(u)
        Nv.remove(u)
        Nv.remove(v)
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
    global PA
    PA = _PA
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
    global PA
    PA = _PA
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
            CostMatrix.add_edge(e[0],e[1],weight=buildCostMatrix3(i,j,e[0],e[1],_GA,PA,_EFA))
    #for u,v,weight in sorted(CostMatrix.edges.data("weight"),key=itemgetter(2,0,1)):print((u.isNeighbor,v.isNeighbor),u.name,v.name,weight)
    #M = nx.max_weight_matching(CostMatrix,maxcardinality=True,weight="weight")
    #assert(nx.is_perfect_matching(CostMatrix,M))
    return CostMatrix

def HungarianSolve(CostMatrix):
    cost = 0
    deleteEdges = nx.max_weight_matching(CostMatrix,maxcardinality=True)
    [cost:= cost+ (-1*CostMatrix[l][r]['weight']) for l,r in deleteEdges if CostMatrix[l][r]['weight'] < 0]
    assert(cost >= 0)
    return cost, deleteEdges

# construct an array that is (n + k) by (n+ k)
def RefineByMatching(_GA,_PA,_EFA,_kA):
    global counter
    PA = _PA
    print("PA before",PA.edges)
    kA = _kA
    GA = _GA
    EFA = _EFA
    #global tbl
    # next time try 
    #global edgeUse
    #counter+=1
    #tbl(counter,datetime.now(),kA,len(PA.edges),len(GA.edges))
    _edgeUse = {(e[0],e[1]):0 for e in list(GA.edges)}
    print(list(GA.edges))
    
    for edge in list(PA.edges)[:]:
        e = tuple(sorted(edge))
        i = e[0]
        j = e[1]
        CostMatrix = buildCostMatrix4(i,j,GA,PA,EFA,kA)
        cost,deleteEdges = HungarianSolve(CostMatrix)
        assert(nx.is_perfect_matching(CostMatrix,deleteEdges))
        deletions = [(u.name,u.parent,v.name,v.parent) for u,v in deleteEdges]
        for d in deletions:
            u = d[0]
            v = d[2]
            if cost > 2*kA:
                assert(tuple(u,v) in PA.edges)
                PA.remove_edge(u,v)
            else:
                _d1 = tuple(sorted((d[0],d[1])))
                _d2 = tuple(sorted((d[2],d[3])))
                print("DELETIONS",d)
                if _d1 in GA.edges:
                    print(_d1)
                    _edgeUse[_d1] += 1
                if _d2 in GA.edges:
                    print(_d2)
                    _edgeUse[_d2] += 1
            
                    
        #s = list(sorted(edgeUse.items(), key=itemgetter(1), reverse=False))
        #print("Maximum:",s[0])
        #print("PRINTING EDGEUSE ARRAY")
        #for edge in s[:25]:print(edge[1],end=", ")
        #print("Done!")
        print("PA after:",PA.edges)
        
        return _edgeUse,PA



# call count_automorphisms_vf2 on G.
def ComputeAutomorphisms(_GA):
    g = ig.Graph.from_networkx(_GA)
    autNum = g.count_automorphisms_vf2()
    return autNum

# Need set E(PA) exposed
# _GA must be g
def DegreeDiffElim(_GA,_PA,_kA):
    PA = _PA
    GA = _GA
    edges = list(PA.edges)
    for edge in edges:
        e = tuple(sorted(edge))
        i = e[0]
        j = e[1]
        if abs(GA.degree(i)-GA.degree(j)) > _kA:
            PA.remove_edge(i,j)
    
    #print("-"*20,"refinement",PA)
    return PA

# dEFA(i) to be the number of fixed edges incident to i
#

    adjacent_fixed = nx.all_neighbors(_EFA,i)
    return len(adjacent_fixed)

def FixedDegElim(_GA,_PA,_EFA):
    PA = _PA
    EFA = _EFA
    edges = list(PA.edges)
    for edge in edges:
        e = tuple(sorted(edge))
        i = e[0]
        j = e[1]
        try:
            if len(nx.all_neighbors(EFA,i)) > GA.degree(j):
                PA.remove_edge(i,j)
        except(nx.exception.NetworkXError):
            print("FIXED DEG ELIM ERROR")
            pass
    return PA


# G = tograph(_PA)
# nx.maximal_independent_set(G)
# return size of preceding
def GreedyIndependentSetSize(_PA):
    PA = _PA
    independent_set =  nx.maximal_independent_set(PA)
    #print(len(independent_set),"$$$$$$$$$")
    return len(independent_set)


#
# edgeUse is a dictionary
#

# in _E(GA) \ _EFA
# if neighborhood of i or of j are not empty (off diagonal)
# return the first encountered
# else prune
def FindBranchEdge(_GA,_PA,_EFA,_edgeUse):
    PA = _PA
    items = _edgeUse.items()
    argmax = max(items, key=itemgetter(1))[0]
    #print("ARGMAX",argmax)
    
    if _edgeUse[argmax] > 0:
        return argmax
    fixed = set(list(_EFA.edges))
    remaining_edges = [tuple(sorted(edge)) for edge in list(_GA.edges) if edge not in fixed]
    #print("Remaining_edges of GA",remaining_edges)
    for edge in remaining_edges:
        e = tuple(sorted(edge))
        if e[0] == e[1]:
            pass
        else:
            i = e[0]
            j = e[1]
            _NPAi = list(nx.all_neighbors(_PA,i))
            _NPAj = list(nx.all_neighbors(_PA,j))
            if len(_NPAi) > 0 or len(_NPAj) > 0:
                return tuple(i,j)
    return -1




def main():
    # used to distinguish dummy variables from vertices.
    logging.basicConfig(stream=sys.stdout,level=logging.INFO)
    # _PA global prior to loop
    global _PA
    global k
    global G
    global changed
    global lowerBound
    global incumbentValue
    global incumbentSolution
    global counter
    global prune
    #
    n = int(sys.argv[1])
    p = float(sys.argv[2])
    G = nx.erdos_renyi_graph(n,p)
    #
    Node = namedtuple('Node',('PA','EDA','EFA','delChild'))
    spacer = " "*4
    T = 0
    t = 0
    #
    logging.info('Begin FindAlmostSymmetry')
    #
    logging.info('Initialize Root Node')
    k = int(sys.argv[3])
    P_R = nx.complete_graph(n)
    EDR = nx.Graph();
    EDR.add_nodes_from(list(G.nodes))
    EFR = nx.Graph();
    EFR.add_nodes_from(list(G.nodes))
    delChild = True
    #
    logging.info("Incumbents")
    incumbentValue = n #|G|
    incumbentSolution = nx.Graph(G) #EDA
    # input to nodestack is iterable
    logging.info("Append root to stack")
    NodeStack = deque()
    NodeStack.append(Node(P_R,EDR,EFR,delChild))
    logging.info("enter loop")
    while True or counter > 20:
        try:
            prune = False
            #
            logging.info('unpack Node data')
            A = NodeStack.pop();
            PA = A.PA
            eda = nx.Graph(A.EDA);
            efa = nx.Graph(A.EFA);
            delChild = A.delChild;
            #
            GA = nx.Graph(G)
            for e in list(eda.edges):
                if e[0] != e[1]:
                    GA.remove_edge(e[0],e[1])
            #
            kA = k - len(eda.edges)
            #
            
            if delChild:
                logging.info(spacer+"Node signal: Delete Edges, compute aut")
                autNum = ComputeAutomorphisms(GA)
                logging.info(spacer+"Autnum computed")
                if autNum < incumbentValue:
                    logging.info(spacer+"update incumbents")
                    incumbentValue = autNum
                    incumbentSolution = eda
                if len(efa.edges) == len(GA.edges) or kA == 0:
                    logging.info(spacer+spacer+"prune: no more edges / budget depleted")
                    break
                logging.info(spacer+"DegreeDiffElim")
                PA = DegreeDiffElim(GA,PA,kA) #update PA
                updated = str(len(PA.edges))
                logging.info(spacer+"updated edges of P_A: "+updated+" members")
            #
            else:
                logging.info(spacer+"Node Signal: Fix Edges, refine P_A")
                if len(efa.edges) == len(GA.edges):
                    logging.info(spacer+spacer+"prune: no more edges")
                    prune = True
                logging.info(spacer+"update PA, FixedDegreeElim")
                PA = FixedDegElim(GA,PA,efa)
                updated = str(len(PA.edges))
                logging.info(spacer+"updated edges of P_A:"+updated+" members")
                _PA = PA
                if prune:
                    continue
            #
            #print(T,".0 lower bound via greedyIndependentSet")
            lowerBound = GreedyIndependentSetSize(PA)
            #
            if lowerBound >= incumbentValue:
                prune = True
                continue
            logging.info(".1 lower bound loop")
            logging.info(spacer+".1.1: While P_A changed by RefinebyMatching")
            pruned_by_lb = False
            edgeUse,PA = RefineByMatching(GA,PA,efa,kA)
            # python-global _PA prior to PA
            _PA = PA
            while changed and not prune:
                lowerBound = GreedyIndependentSetSize(PA)
                logging.info(spacer+spacer+". computing bound")
                if lowerBound >= incumbentValue:
                    logging.info(spacer+spacer+" lowerbound >= incumbentValue: prune node")
                    prune = True
                else:
                    _edgeUse,PA = RefineByMatching(GA,PA,efa,kA)
                    if not len(_PA.edges) == len(PA.edges):
                        changed = True
                        edgeUse = _edgeUse
                        _PA = PA
                    else:
                        changed = False
                    #print(T,".2 new lower bound!")
                
            logging.info(".3 FindBranchEdge")
            branchEdge = FindBranchEdge(GA,_PA,efa,edgeUse)
            logging.info(spacer+".3.1 bounding is concluded")
            # 
            if branchEdge == -1:
                logging.info(space+".3.2 branchEdge signal: prune")
                prune = True
                continue
            #
            logging.info(".4 create children")
            b = tuple(sorted((branchEdge[0],branchEdge[1])))
            efa_branch = efa.add_edge(b[0],b[1])
            eda_branch = eda.add_edge(b[0],b[1])
            logging.info(spacer+".4.1: append fChild")
            NodeStack.append(Node(_PA,eda,efa_branch,False))
            logging.info(spacer+".4.2: append dChild")
            NodeStack.append(Node(_PA,eda_branch,efa,True))            
        #
        except IndexError:
            logging.info(".5: end of stack")
            break
        
    print("incumbentValue",incumbentValue)
    print("incumbentSolution",incumbentSolution)
    

   
            
    
    
if __name__ == "__main__":
    main()
