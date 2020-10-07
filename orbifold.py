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

from numpy.random import permutation
from colorama import Fore
from colorama import Style
from colorama import Back
from operator import itemgetter

def asciiTable(input_list1,rowsize,colsize,theta):
    # stackoverflow recipe #
    rows = rowsize
    cols = colsize
    content = [["."]*cols for _ in range(rows)]

    grid =input_list1
    
    for (y,x,c) in grid: content[y][x] = c

    # build frame
    width       = len(str(max(rows,cols)-1))
    contentLine = "#|values|"

    dashes      = "-".join("-"*width for _ in range(cols))
    frameLine   = contentLine.replace("values",dashes)
    frameLine   = frameLine.replace("#"," "*width)
    frameLine   = frameLine.replace("| ","+-").replace(" |","-+")

    # print grid
    #print(frameLine)
    out = []
    for i,row in enumerate(reversed(content),1):
        values = "".join(f'{Fore.BLACK}{Back.BLUE} {v} {Style.RESET_ALL}' if v%2==1 else  f'{Fore.BLACK}{Back.CYAN} {v} {Style.RESET_ALL}'for v in it.islice(it.cycle(row),theta,theta+len(row)))
        line = contentLine.replace("values",values)
        line = line.replace("#",f"{rows-i:{width}d}")
        print(30*" ",line)
    #print(frameLine)

    # x-axis numbers
    numLine = contentLine.replace("|"," ")
    numLine = numLine.replace("#"," "*width)
    colNums = " ".join(f"{i:<{width}d}" for i in range(cols))
    numLine = numLine.replace("values",colNums)
    #print(numLine)

# this should clearly display the evolution #
def permutationTable(input_list1,rowsize,colsize,theta,k):
    # stackoverflow recipe #
    margin = 0
    rows = rowsize
    cols = colsize
    content = [["."]*cols for _ in range(rows)]

    grid =input_list1
    
    for (y,x,c) in grid: content[y][x] = c

    # build frame
    width       = len(str(max(rows,cols)-1))
    contentLine = "#:values:values:values:"

    dashes      = ".".join("."*width for _ in range(cols))
    frameLine   = contentLine.replace("values",dashes)
    frameLine   = frameLine.replace("#","."*width)
    frameLine   = frameLine.replace("|",".").replace("|",".")

    # print grid
    print(margin*" ",frameLine,"budget = ",k)
    out = []
    for i,row in enumerate(reversed(content),1):
        values = "".join(f'{Fore.RED}{Back.BLACK} {"X"} {Style.RESET_ALL}' if v<0 else  f'{Fore.BLACK}{Back.WHITE} {v} {Style.RESET_ALL}'for v in it.islice(it.cycle(row),theta,theta+len(row)))
        line = contentLine.replace("values",values)
        line = line.replace("#",f"{rows-i:{width}d}")
        print(margin*" ",line)
   # print(margin*" ",frameLine)

    # x-axis numbers
    #numLine = contentLine.replace("|"," ")
    #numLine = numLine.replace("#","   "*width)
    #colNums = " ".join(f"{i:<{width}d}" for i in range(cols))
    #numLine = numLine.replace("values",colNums)
    #print((margin-3)*" ",numLine)
    
def printMatchingRows(g,dim,LeftBp,RightBp):
    for n in LeftBp:print(f"{n[0]:>3}","   ",end=" ")
    print("\n")
    for nv in reversed(LeftBp):print(f"{nv[1]}","   ",end=" ")
    print("\n")
    for m in RightBp:print(f"{m[0]:>3}","   ",end=" ")
    print("\n")
    for mv in RightBp:print(f"{mv[1]}"," ",end=" ")
    print("\n")

    for e in it.product(LeftBp,RightBp):print("          "*dim[1],"    "*dim[0],f"{Fore.WHITE}{e[0][1]}{Style.RESET_ALL}   {e[1][1]}  ",
                                              f'{Fore.MAGENTA}{abs(assignweight(g,e[0][1],e[1][1])):>5}{Style.RESET_ALL}' if type(e[0][1]) is tuple and type(e[1][1]) is tuple else f'{Fore.CYAN}{abs(assignweight(g,e[0][1],e[1][1])):>5}{Style.RESET_ALL}')

def printColumns(dim,Nu,Nv):
    print(dim[0])
    for u in Nu:print(f"{u:>3}","  ",end=" ")
    print("\n")
    print(dim[1])
    for v in Nv:print(f"{v:>3}","  ",end=" ")
    print("\n")
    print(f"{Fore.WHITE}  .  {Style.RESET_ALL}"*800,end = str(dim[0])+" + "+str(dim[1])+" by ("+str(dim[0])+"+"+str(dim[1])+")^2\n")
    
def printNeuronStructure(g,shapeLeft,shapeRight):
    for j in it.takewhile(lambda j: j[0][1]==0, it.product(shapeLeft,shapeRight)):
        w = assignweight(g,j[0][0]%2==0,j[0][1],j[1][0]%2==0,j[1][1])
        print(" "*100,f'{Fore.BLUE}{j[0][0]:>08b}{Style.RESET_ALL}' if not j[0][0]%2 else f'{Fore.RED}{j[0][0]:>08b}{Style.RESET_ALL}'," "*1,
                                                                      f'{Fore.GREEN}{j[0][1]:>08b}{Style.RESET_ALL}' if j[0][1]==0 else f'{Fore.WHITE}{j[0][1]:>08b}{Style.RESET_ALL}'," "*1,
                                                                      f'{Fore.BLUE}{j[1][0]:>08b}{Style.RESET_ALL}' if not j[1][0]%2 else f'{Fore.RED}{j[1][0]:>08b}{Style.RESET_ALL}'," "*1,
                                                                      f'{Fore.GREEN}{j[1][1]:>08b}{Style.RESET_ALL}' if j[1][1]==0 else f'{Fore.WHITE}{j[1][1]:>08b}{Style.RESET_ALL}'," "*1,f'{Fore.MAGENTA}{"+"*w}{Style.RESET_ALL}')
def printMatchingStructure(M): 
    for edge in M:print(f'{Fore.BLUE}{(lambda x: x[0][1])(e):>3}{Style.RESET_ALL}' if edge[0][0] == 0 else f'{Fore.RED}{(lambda y: y[1][1])(e):>3}{Style.RESET_ALL}','   ',end='')
    print('\n')


##################################################################################################################################################################################################################################
##################################################################################################################################################################################################################################
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

#def LoopOverPairs(names,n,g):
    # produces:
    # memo1 : {score:pair}
    # memo2 : {dim2matching}
    # Rho :
    # for uv permutations in ij :
    #    Nu,Nv,dim,num2namepoint = getMunkresBasis(names,n,g,u,v)
    #    M = hungarianflow(g,num2namepoint,Nu,Nv)
    #    score = computeWeightSum(M)
    #    
    #
    #
    #
    #

def main():
    
    ########################################################################## FILE IO # 
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
    
    
    ###################################################################### Random Graph #
    
    Tau = int(sys.argv[3])
    P0 = []
    #memo1 = {} 
    #memo2 = {}
    Rho = {} # Rho is indexed by budget
    
    ######################################################################## Generate P #
    print("Generate P")
    
    # edges sorted here O(n^2)
    for u,v in it.product(range(n),range(n)):
        num2namepoint,Nu,Nv,dim = genMunkresBasis(names,n,g,u,v)
        #
        B,M  = hungarianSolve(g,num2namepoint,Nu,Nv,dim)
        print(u,v," perfect matching: ",nx.is_perfect_matching(B,M))
        linsum = extractweight(g,u,v)
        for m in M:linsum+=extractweight(g,m[0][1],m[1][1])
        print(dim)
        P0.append((u,v,linsum))
        #memo1[linsum] = (u,v)
        #memo2[dim] = M
        # return memo1 , memo2 #
    
    print("P generated on ",n,"nodes")

    ############################################# Refinement Algorithm 1: DegreeDifElim #
    P1 = [ (e[0],e[1],0) if abs(g.degree(e[0]) - g.degree(e[1])) > Tau else e for e in P0 ]


    ############################################# Refinement Algorithm 2: JimElim #######
    print("Jim Elim")
    Rho[Tau+1] = P1
    Rho2 = {}
    killed = {}
    eliminations = 0
    
    # this encodes how many times JimElim will have to repeat itself by length/index, and the respective elimination count
    _worstcase = []
    
    t = 0
    # casualties may not exceed dimensions of P
    # however since casualties increase exponentially, this should not be a huge computational burden
    
    while eliminations < n*(n-1):
        _hc = []
        epsilon = Tau- t
        # Eliminate if over budget
        Rho[epsilon] = [ (e[0],e[1],-1) if e[2]  > 2*epsilon else e for e in Rho[epsilon+1] ]

        # zero matrix
        Rho2[epsilon] = []
        Delta = [(e[0],e[1],0) for e in Rho[epsilon]]
        for e in Rho[epsilon]:print(e[2]," ",sep=' ')
        print('\n')
        #
        # add delta to Rho (while pairs memoized ) O(e*n)
        for e in Rho[epsilon]:
            print(e[2],sep=' ')
            print('\n')
        ######################################################## compute delta sum
            if e[2] == -1:
                for i in range(n):
                    #L
                    Delta[i*n+e[1]] = (e[0],e[1],e[2]+2)
                    Delta[e[0]*n+i] = (e[0],e[1],e[2]+2)
                    #R
                    Delta[i*n + e[0]] = (e[0],e[1],e[2]+2)
                    Delta[e[1]*n+i] =(e[0],e[1],e[2]+2)
                
                killed[(e[0],e[1])] = Delta
                Rho2[epsilon].append(Delta)
                
       #########################################################

       # Rho+delta = lowerbound
       # 1 , 2, 3, .... .        n
       # # # # L'# ##L # # # # # # R+L=n cost matrix
       #       #
       # v2d   #  v2v
       # # # # L- # #L # # # # # #
       #       #  2   ^    ^ d2v after
       # d2d   #  < # # # # # # #
       #       #    #    d2v before
       #       #  < #                 cost matrix realignment, creates exactly 1 more '2'.
       #       #    #
       #       R+   R                 |v2d| = L^2, |d2V| = R^2
       
       
       # 
        eliminations = len(killed.keys())
        _worstcase.append(eliminations)
        print(eliminations," Eliminations at budget: ",epsilon)
        #increment counter
        t+=1
    del Rho[Tau+1]
    
    # Tau decreases, t increases # ######################################################## Jim Elim 
    #                                   #
    #  Rho[epsilon] contains sparse     # 
    #                                   #
    #  Rho2 contains rich               #
    #                                   #
    #  # # # # # # ## # # # # # # # # # #

    #                                               #
    #  Base image indexed by threshold              #
    #  Base image + sum of deltas  = sparse         #
    #                                               #

    #
    # BASE[epsilon] = P at budget epsilon
    #
    BASE = {}
    for ctr in range(Tau):
        epsilon = Tau - ctr
        _matrix = [[0]*n for _ in range(n)] # n x n table
        for (y,x,score) in Rho[epsilon]:_matrix[y][x] = score
        BASE[epsilon] = _matrix 
    #
    # DELTAS[epsilon] = [ Delta0, Delta1, .... ]
    #
    DELTAS = {}
    for ctr in range(Tau):
        epsilon  = Tau - ctr
        _Deltas = []
        for _DeltaList in Rho2[epsilon]:
            # make single matrix
            _matrix = [[0]*n for _ in range(n)]
            for (y,x,score) in _DeltaList:_matrix[y][x] = score
            # append to _Deltas
            _Deltas.append(_matrix)
        DELTAS[epsilon] = _Deltas
    
        


    ############################################# Algorithm 3 ###########################
    #
    # Apply multidimensional scaling
    # 
    #
    ############################################# matplot lib ############################

    from matplotlib import pyplot as plt
    from matplotlib.collections import LineCollection

    from sklearn import manifold
    from sklearn.metrics import euclidean_distances
    from sklearn.decomposition import PCA

    from sklearn.manifold import MDS
    model = MDS(n_components=2, dissimilarity='precomputed', random_state=1)
    D = LowerBoundMatrix[1]
    
    #out = model.fit_transform(D)
    #plt.scatter(out[:, 0], out[:, 1], **colorize)
    #plt.axis('equal');
    
    
        


if __name__ == "__main__":
    main()
#
#
#
# print neuron-cosets based on memo2
