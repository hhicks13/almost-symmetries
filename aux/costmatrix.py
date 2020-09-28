#!/usr/bin/env python
from array import*
import sys
from itertools import chain
import networkx as nx
import numpy as np
import random
import itertools as it
import networkx as nx
import re
from fractions import Fraction
from numpy.random import permutation
from colorama import Fore
from colorama import Style
from operator import itemgetter




# auxiliary / plotting #
def bipartiteCoordinates(leftSize,rightSize):
    dims = [leftSize,rightSize]
    bpcoords = [(idi, idv) for idi, dim in enumerate(dims) for idv in range(dim)]

def matrix2coord(input_matrix,shape):
    grid = [(idx[0],idx[1],input_matrix[idx[0]][idx[1]]) for idx in it.product(*[range(s) for s in shape])]
    return grid

def asciiTable(input_matrix,rowsize,colsize):
    # prepare the empty content
    rows = rowsize
    cols = colsize
    content = [["."]*cols for _ in range(rows)]

    grid = matrix2coord(input_matrix,[rowsize,colsize])

    # assign values at coordinates as needed (based on your grid)

    #grid = [(4,1,"H"),(6,3,"L"),(5,2,"E"),(4,6,"R"),(7,4,"L"),(6,6,"W"),(3,6,"L"),(2,6,"D"),(5,6,"O")]
    for (y,x,c) in grid: content[y][x] = c

    # build frame
    width       = len(str(max(rows,cols)-1))
    contentLine = "# | values |"

    dashes      = "-".join("-"*width for _ in range(cols))
    frameLine   = contentLine.replace("values",dashes)
    frameLine   = frameLine.replace("#"," "*width)
    frameLine   = frameLine.replace("| ","+-").replace(" |","-+")

    # print grid
    print(frameLine)
    for i,row in enumerate(reversed(content),1):
        values = " ".join(f"{v:{width}d}" for v in row)
        line = contentLine.replace("values",values)
        line = line.replace("#",f"{rows-i:{width}d}")
        print(line)
    print(frameLine)

        # x-axis numbers
    numLine = contentLine.replace("|"," ")
    numLine = numLine.replace("#"," "*width)
    colNums = " ".join(f"{i:<{width}d}" for i in range(cols))
    numLine = numLine.replace("values",colNums)
    print(numLine)



def printlexBinaryPeriod(g,shapeLeft,shapeRight):
    for j in it.takewhile(lambda j: j[0][1]==0, it.product(shapeLeft,shapeRight)):
        w = weight(g,j[0][0]%2==0,j[0][1],j[1][0]%2==0,j[1][1])
        print(" "*100,f'{Fore.BLUE}{j[0][0]:>08b}{Style.RESET_ALL}' if not j[0][0]%2 else f'{Fore.RED}{j[0][0]:>08b}{Style.RESET_ALL}'," "*1,
                                                                      f'{Fore.GREEN}{j[0][1]:>08b}{Style.RESET_ALL}' if j[0][1]==0 else f'{Fore.WHITE}{j[0][1]:>08b}{Style.RESET_ALL}'," "*1,
                                                                      f'{Fore.BLUE}{j[1][0]:>08b}{Style.RESET_ALL}' if not j[1][0]%2 else f'{Fore.RED}{j[1][0]:>08b}{Style.RESET_ALL}'," "*1,
                                                                      f'{Fore.GREEN}{j[1][1]:>08b}{Style.RESET_ALL}' if j[1][1]==0 else f'{Fore.WHITE}{j[1][1]:>08b}{Style.RESET_ALL}'," "*1,f'{Fore.MAGENTA}{"+"*w}{Style.RESET_ALL}')

#
#
# neuron is going to be a list of tuples #
#
#

def getNameSpace(n,names):
    grid0 = [range(n),names]
    num2namepoint = {j[0]:j for j in it.zip_longest(*grid0)}

    if len(names) == n:
        return num2namepoint
    else: return None

#
# def neuron2factors(g,neuron):
#
#

def namepointlist(Nu,namespace):
    namepoints = [namespace[u] for u in Nu]
    return namepoints


def pair2neighborhood(g,u,v):
    Nu = [e for e in nx.neighbors(g,u)]
    Nv = [f for f in nx.neighbors(g,v)]
    try:
        Nu.remove(v)
        Nv.remove(u)
    except(ValueError,TypeError):
        print("- no duplicate - neighbors computed")

    lN,rN = sortByLength(Nu,Nv)
    dim = [len(lN),len(rN)]
    return lN,rN,dim

def sortByLength(Nu,Nv):
    leftN = []
    rightN = []
    if len(Nu) == len(Nv):
        return Nu, Nv
    elif len(Nu) < len(Nv):
        return Nu, Nv
    else:
        return Nv,Nu



# type inference for weighting
def weight(g,u,v):
    # both tuple
    if type(u) is tuple and type(v) is tuple:
        return abs(g.degree(u[0]) - g.degree(v[0]))    
    elif type(u) is int and type(v) is int:
        return 0
    else:
        return 1


# neither tuple

#def pair2hungarianScore(g,u,v):

#def pair2neuron(g,u,v):

#def pair2matching(g,u,v):

#def computeLR(g,u,v):

#def bcm(g,u,v):

#def get_value(v):
#    try:
#        return float(v.split(' ')[1])
#    except IndexError as e:
#        return int(v[1:])

def main():

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
    u = int(sys.argv[3])
    v = int(sys.argv[4])

    # namespace #
    num2namepoint = getNameSpace(n,names[:n])
    Nu,Nv,dim = pair2neighborhood(g,u,v)
    print(dim[0])
    for u in Nu:print(f"{u:>3}","  ",end=" ")
    print("\n")
    print(dim[1])
    for v in Nv:print(f"{v:>3}","  ",end=" ")
    print("\n")
    print(f"{Fore.WHITE}  .  {Style.RESET_ALL}"*800,end = str(dim[0])+" + "+str(dim[1])+" by ("+str(dim[0])+"+"+str(dim[1])+")^2\n")

    # build cost matrix / solve weighted assignment # Left then right
    #

    # for jim 
    LeftBp = [ (0,num2namepoint[idv]) if idi == 0 else (0,idv) for idi, vertices in enumerate([Nu,[i for i in range(dim[1]) ]]) for idv in vertices] # counter intuitive, but by doing the same operation twice you save more time

    RightBp =[ (1,num2namepoint[idv]) if idi == 0 else (1,idv) for idi, vertices in enumerate([Nv,[i for i in range(dim[0]) ]]) for idv in vertices] #
    #
    
    for n in LeftBp:print(f"{n[0]:>3}","   ",end=" ")
    print("\n")
    for nv in LeftBp:print(f"{nv[1]}"," ",end=" ")
    print("\n")
    for m in RightBp:print(f"{m[0]:>3}","   ",end=" ")
    print("\n")
    for mv in RightBp:print(f"{mv[1]}"," ",end=" ")
    print("\n")

    for e in it.product(LeftBp,RightBp):print("          "*dim[1],"    "*dim[0],f"{Fore.WHITE}{e[0][1]}{Style.RESET_ALL}   {e[1][1]}  ",
                                              f'{Fore.MAGENTA}{weight(g,e[0][1],e[1][1]):>5}{Style.RESET_ALL}' if type(e[0][1]) is tuple and type(e[1][1]) is tuple else f'{Fore.CYAN}{weight(g,e[0][1],e[1][1]):>5}{Style.RESET_ALL}')
    
    # for neuron-prime
    #vertexCoords = [(idi,num2namepoint[idv]) for idi, node in enumerate([Nu,dim]) for idv in node] # one pass
    #dummyCoords = [(idi, idv) if idi == 0 else (idi, str(idv)) for idi, dim in enumerate(reversed(dim)) for idv in range(dim)] # one pass
    #
    #

    #
    # Format Edges for networkx
    
    #weighted_edges = [ (e[0],e[1],weight(e,e[0][1],e[1][1])) for e in it.product(LeftBp,RightBp) ] #O(n^2)

    # weight graph
    #
    #LeftNode = []
    #RightNode = []
    B = nx.Graph()

    # Format Nodes for Networkx
    #for v in vertexCoords:
    #    print(v)
    #    if v[0] == 0:
    #        LeftNode.append(v[1])
    #    else:
    #        RightNode.append(v[1])

    #for d in transposeCoords:
    #    if d[0] == 0:
    #        LeftNode.append(d[1])
    #    else:
    #        RightNode.append(d[1])
    
    #B.add_nodes_from(LeftBp,bipartite=0)
    #B.add_nodes_from(RightBp,bipartite=1)
    #B.add_weighted_edges_from(
    #[(e[0][1],e[1][1],weight(g,e[0][1],e[1][1])) for e in it.product(vertexCoords,transposeCoords) ],weight='weight')
    
    #for r in RightBp:print(r)
    #for l in LeftBp:print(l)
    #for f in transposeCoords:print(f)

    
    
    
if __name__ == "__main__":
    main()
