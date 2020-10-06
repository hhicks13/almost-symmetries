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


# memo = 1wf = schedule

################################################################### RECONSTRUCTORS  ######################### 

#def scoreFunction():
    # reconstructs scoreTable from memo1, because scores can collide due to orbit property a*b in H implies a and b in H 

#def cosets(dim,costMatrix):
    # reconstructs cost matrices from memo2, because dims can collide, produces likelihood space
    # populates with matchings

#def orbifold(scoreTable):
    # assumed that (a,b,f(ab)) is closed:
    # plots scoreTable. 

######################################################################################## TUPLE CONSTRUCTORS

#def bipartiteCoordinates(leftSize,rightSize):
#    dims = [leftSize,rightSize]
#    bpcoords = [(idi, idv) for idi, dim in enumerate(dims) for idv in range(dim)]

#def matrix2coord(input_matrix,shape):
#    grid = [(idx[0],idx[1],input_matrix[idx[0]][idx[1]]) for idx in it.product(*[range(s) for s in shape])]
#    return grid

###########################################################################################################

def phi(n):
    amount = 0        
    for k in range(1, n + 1):
        if fractions.gcd(n, k) == 1:
            amount += 1
    return amount

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

########################################################################################### Printers  ###################################################################################################################################


#def linearsum(B,M):
#    linearsum = 0;
#    for m in M:
#        e = (m[0][1],m[1][1])
#        print(e)
        #linearsum+=abs(B.get_edge_data(*e))
        #print(linearsum)
   # return M,linearsum
        
        
    #computes the sum of weights
#
#
#
#
#
# 
#
##################################################################################################################################################################################################################################

# unfinished ideas

# memo3 -  pair2neuron(g,u,v):
# for neuron-prime
    #vertexCoords = [(idi,num2namepoint[idv]) for idi, node in enumerate([Nu,dim]) for idv in node] # one pass
    #dummyCoords = [(idi, idv) if idi == 0 else (idi, str(idv)) for idi, dim in enumerate(reversed(dim)) for idv in range(dim)] # one pass

# neuron is going to be list of tuples, reconstructed from matchings, want to get neuron cosets for CNN

##################################################################################################################################################################################################################################
##################################################################################################################################################################################################################################
##################################################################################################################################################################################################################################
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
    u = int(sys.argv[3])
    v = int(sys.argv[4])
    ###################################################################### Random Graph #
    
    Tau = int(sys.argv[5])
    P0 = []
    memo1 = {} 
    memo2 = {}
    Rho = {} # Rho is indexed by budget
    
    ######################################################################## Generate P #
    for u,v in it.product(range(n),range(n)):
        num2namepoint,Nu,Nv,dim = genMunkresBasis(names,n,g,u,v)
        #flow
        B,M  = hungarianSolve(g,num2namepoint,Nu,Nv,dim)
        #print(M)
        #print(nx.is_perfect_matching(B,M))
        linsum = extractweight(g,u,v)
        for m in M:linsum+=extractweight(g,m[0][1],m[1][1])
        #print(linsum)
        P0.append((u,v,linsum))
        memo1[linsum] = (u,v)
        #memo2[dim] = M
        # return memo1 , memo2 #

    ############################################# Refinement Algorithm 1: DegreeDifElim #
    P1 = [ (e[0],e[1],0) if abs(g.degree(e[0]) - g.degree(e[1])) > k else e for e in P0 ]


    ############################################# Refinement Algorithm 2: JimElim #######
    Rho[k+1] = P1
    Rho2 = {}
    memo = {}
    casualties = 0
    t = 0
    # change budget to survivors 
    while casualties < n*(n-1):
        
        # make permutation graph
        epsilon = Tau - t
        Rho[epsilon] = [ (e[0],e[1],-1) if e[2]  > 2*epsilon else e for e in Rho[epsilon+1] ]
        # zero matrix 
        delta = [(e[0],e[1],0) for e in Rho[epsilon]]
        #for e in Rho[budget]:print(e)
        
        #pairs memoized
        for e in Rho[epsilon]:
            if e[2] == -1:
                for i in range(n):
                    delta[i*n+e[1]] = (e[0],e[1],e[2]+2)
                    delta[e[0]*n+i] = (e[0],e[1],e[2]+2)
                memo[(e[0],e[1])] = delta
        #check to make sure copy works
        
        casualties = len(memo.keys())
        print(casualties)
        t+=1
    del Rho[k+1]


    ############################################# Algorithm 3 ###########################
    # how does the graph change if you remove an edge.
    #
    # R = BinaryCopy(Rho[epsilon])
    #
    # w_ij = 1 - |Ni intersect Nj |/| Ni union Nj | <--- distance from center
    # 
    # Distance Matrix
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

    EPSILON = np.finfo(np.float32).eps
    n_samples = 20
    seed = np.random.RandomState(seed=3)
    X_true = seed.randint(0, 20, 2 * n_samples).astype(np.float)
    X_true = X_true.reshape((n_samples, 2))
    # Center the data
    X_true -= X_true.mean()

    similarities = euclidean_distances(X_true)

    # Add noise to the similarities
    noise = np.random.rand(n_samples, n_samples)
    noise = noise + noise.T
    noise[np.arange(noise.shape[0]), np.arange(noise.shape[0])] = 0
    similarities += noise

    mds = manifold.MDS(n_components=2, max_iter=3000, eps=1e-9, random_state=seed,
                   dissimilarity="precomputed", n_jobs=1)
    pos = mds.fit(similarities).embedding_

    nmds = manifold.MDS(n_components=2, metric=False, max_iter=3000, eps=1e-12,
                    dissimilarity="precomputed", random_state=seed, n_jobs=1,
                    n_init=1)
    npos = nmds.fit_transform(similarities, init=pos)

    # Rescale the data
    pos *= np.sqrt((X_true ** 2).sum()) / np.sqrt((pos ** 2).sum())
    npos *= np.sqrt((X_true ** 2).sum()) / np.sqrt((npos ** 2).sum())

    # Rotate the data
    clf = PCA(n_components=2)
    X_true = clf.fit_transform(X_true)

    pos = clf.fit_transform(pos)

    npos = clf.fit_transform(npos)

    fig = plt.figure(1)
ax = plt.axes([0., 0., 1., 1.])

s = 100
plt.scatter(X_true[:, 0], X_true[:, 1], color='navy', s=s, lw=0,
            label='True Position')
plt.scatter(pos[:, 0], pos[:, 1], color='turquoise', s=s, lw=0, label='MDS')
plt.scatter(npos[:, 0], npos[:, 1], color='darkorange', s=s, lw=0, label='NMDS')
plt.legend(scatterpoints=1, loc='best', shadow=False)

similarities = similarities.max() / (similarities + EPSILON) * 100
np.fill_diagonal(similarities, 0)
# Plot the edges
start_idx, end_idx = np.where(pos)
# a sequence of (*line0*, *line1*, *line2*), where::
#            linen = (x0, y0), (x1, y1), ... (xm, ym)
segments = [[X_true[i, :], X_true[j, :]]
            for i in range(len(pos)) for j in range(len(pos))]
values = np.abs(similarities)
lc = LineCollection(segments,
                    zorder=0, cmap=plt.cm.Blues,
                    norm=plt.Normalize(0, values.max()))
lc.set_array(similarities.flatten())
lc.set_linewidths(np.full(len(segments), 0.5))
ax.add_collection(lc)

plt.show()
    
    
        


if __name__ == "__main__":
    main()
#
#
#
# print neuron-cosets based on memo2
