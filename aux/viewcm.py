#!/usr/bin/env python
from array import*
import sys
from itertools import chain
import networkx as nx
import numpy as np
import random
import itertools as it
import networkx as nx
from fractions import Fraction
from numpy.random import permutation
from colorama import Fore
from colorama import Style
from operator import itemgetter
#
#
# compute
#
#
def farey_sequence(n: int, descending: bool = False) -> None:
    """Print the n'th Farey sequence. Allow for either ascending or descending."""
    (a, b, c, d) = (0, 1, 1, n)
    if descending:
        (a, c) = (1, n - 1)
    print("{0}/{1}".format(a, b))
    while (c <= n and not descending) or (a > 0 and descending):
        k = (n + b) // d
        (a, b, c, d) = (c, d, k * c - a, k * d - b)
        print("{0}/{1}".format(a, b))

def _symmetric_aspect(x):
    #euclids algorithm
    factors = []
    for i in range(1, x + 1):
        if x % i == 0:
            factors.append(i)

    return factors

def bcmblock(n):
    ar = random.randint(1,n-1)
    p = [i for i in permutation(n)]
    S = [p[0:ar],p[ar:n]]
    dim = [len(s) for s in S]
    dimz = [dim,dim]

def coords(shape):
    coords = [(idi, idv) for idi, dim in enumerate(shape) for idv in range(dim)]
    return coords

def makeLR(b):
    print("making LR")
    dim1 = b
    dim2 = dim1[::-1]
    print(dim1)
    print(dim2)
    L = coords(dim1)
    r = coords(dim2)
    R = [(node[0]+2,node[1]) for node in r]
    return L,R

def as_string(seq_of_rows):
    return '\n'.join(''.join(str(i).center(5) for i in row) for row in seq_of_rows)

def weight(g,a,b,c,d):
    if a == 1 and c ==1:
        return abs(g.degree(b)-g.degree(d))
    elif a == 0 and c == 0:
        return 0
    else:
        return 1

def pulseWidth(l,r):
    ctr = 1
    total= 0
    sequence = []
    for j in it.takewhile(lambda j: j != None, it.product(l,r)):
        total +=1
        if j[1][1] != 0:
            ctr+=1
        else:
            sequence.append(ctr)
            ctr=1
    print(sequence)

def upperBound(n,g):
    # factorize maxdegs
    degreefactors = [(idi,idv) for idi, u in enumerate(g.nodes) for idv in _symmetric_aspect(g.degree[u])]
    print("possible integer aspect ratios")
    # stack overflow recipe:
    #empty content:
    gridpoints = list(degreefactors) #list of tuples#


def neighborhoods(g,u,v):
    Nu = [e for e in nx.neighbors(g,u)]
    Nv = [f for f in nx.neighbors(g,v)]
    dim = [len(max(Nu,Nv)),len(min(Nu,Nv))]
    try:
        Nu.remove(v)
        Nv.remove(u)
    except(ValueError,TypeError):
        print("- no duplicate - neighbors computed")
    return min(Nu,Nv),max(Nu,Nv),dim

#
# plot
#

def printMC(Nx):
    for i in it.takewhile(lambda i: i != None, Nx):
        print(f'{Fore.CYAN}{i[0]}{Style.RESET_ALL}    {Fore.MAGENTA}{i[1]}{Style.RESET_ALL}' if i[0]%2 else f'{Fore.MAGENTA}{i[0]}{Style.RESET_ALL}    {Fore.CYAN}{i[1]}{Style.RESET_ALL}')

def printlexAll(g,l,r):
    for j in it.takewhile(lambda j: j != None, it.product(l,r)):
        w = weight(g,j[0][0]%2==0,j[0][1],j[1][0]%2==0,j[1][1])
        print(" "*100,f'{Fore.BLUE}{j[0][0]:>08b}{Style.RESET_ALL}' if not j[0][0]%2 else f'{Fore.RED}{j[0][0]:>08b}{Style.RESET_ALL}'," "*1,
                                                                      f'{Fore.GREEN}{j[0][1]:>08b}{Style.RESET_ALL}' if j[0][1]==0 else f'{Fore.WHITE}{j[0][1]:>08b}{Style.RESET_ALL}'," "*1,
                                                                      f'{Fore.BLUE}{j[1][0]:>08b}{Style.RESET_ALL}' if not j[1][0]%2 else f'{Fore.RED}{j[1][0]:>08b}{Style.RESET_ALL}'," "*1,
                                                                      f'{Fore.GREEN}{j[1][1]:>08b}{Style.RESET_ALL}' if j[1][1]==0 else f'{Fore.WHITE}{j[1][1]:>08b}{Style.RESET_ALL}'," "*1,f'{Fore.MAGENTA}{"+"*w}{Style.RESET_ALL}')
def printlexV2V(g,l,r):
    for j in it.takewhile(lambda j: j != None, it.product(l,r)):
        if j[0][0] == 1 or j[1][0] == 3:continue
        w = weight(g,j[0][0]%2==0,j[0][1],j[1][0]%2==0,j[1][1])
        print(" "*100,f'{Fore.BLUE}{j[0][0]}{Style.RESET_ALL}' if not j[0][0]%2 else f'{Fore.RED}{j[0][0]}{Style.RESET_ALL}'," "*1,
                                                                      f'{Fore.WHITE}{j[0][1]}{Style.RESET_ALL}' if j[0][1]%2 else f'{Fore.CYAN}{j[0][1]}{Style.RESET_ALL}'," "*1,
                                                                      f'{Fore.BLUE}{j[1][0]}{Style.RESET_ALL}' if not j[1][0]%2 else f'{Fore.RED}{j[1][0]}{Style.RESET_ALL}'," "*1,
                                                                      f'{Fore.GREEN}{j[1][1]}{Style.RESET_ALL}' if j[1][1]==0 else f'{Fore.WHITE}{j[1][1]}{Style.RESET_ALL}'," "*1,f'{Fore.MAGENTA}{"+"*w}{Style.RESET_ALL}')
def main():
    for i, arg in enumerate(sys.argv):
        print("<{}>  {}\n".format(i,arg))

    ## graph ##
    n = int(sys.argv[1])
    p = float(sys.argv[2])
    g = nx.erdos_renyi_graph(n,p)

    # compute aspect ratio for symmetric columns #
    gn = g.nodes
    ge = g.edges
    F = Fraction(len(gn),len(ge))
    numfactors = _symmetric_aspect(F.numerator)
    denomfactors = _symmetric_aspect(F.denominator)
    print(f"|G| = {len(g.nodes)}, |E| = {len(g.edges)}")
    print("|G| factors, |E| factors")
    print(min(numfactors,denomfactors))
    print(max(numfactors,denomfactors))
    print(len(gn),len(ge),F)

    ## neighborhood ##
    u = int(sys.argv[3])
    v = int(sys.argv[4])
    Nu,Nv,dim = neighborhoods(g,u,v)

    # compute aspect ratio for symmetric rows #
    print("degrees:")
    gdu = g.degree[u]
    gdv = g.degree[v]
    print(gdu)
    print(gdv)

    F2 = Fraction(gdu,gdv)
    #_symaspect returns array of factors #
    denomfactors2 = _symmetric_aspect(F2.denominator)


    ## possible integer aspect ratios ## memoization
    piar = {(idi,idv):(idi,idv) for idi, u in enumerate(g.nodes) for idv in _symmetric_aspect(g.degree[u])}
    print("possible integer aspect ratios")
    # stack overflow recipe:
    #empty content:
    rows = n
    cols = n
    content = [["."]*cols for _ in range(rows)]
    gridpoints = list(piar) #list of tuples#

    # address by farey l8r
    for (x,y) in gridpoints:
        a = min(x,y)
        b = max(x,y)
        fa = Fraction(a,b)
        content[a][b] = str(a)

    #build frame
    width = len(str(max(rows,cols)-1))
    contentLine = "# | values |"
    dashes = "-".join("-"*width for _ in range(cols))
    frameLine   = contentLine.replace("values",dashes)
    frameLine   = frameLine.replace("#"," "*width)
    frameLine   = frameLine.replace("| ","+-").replace(" |","-+")

    #print frame
    print(frameLine)
    for i,row in enumerate(reversed(content),1):
        values = " ".join(f"{v:{width}s}" for v in row)
        line   = contentLine.replace("values",values)
        line   = line.replace("#",f"{rows-i:{width}d}")
        print(line)

    #dots = " ".join(f'{".":{width}s}' for _ in range(cols))
    #symLine = contentLine.replace("values",dots)
    #symLine = symLine.replace("#"," "*width)
    # do not delete or clean any of this #
    #print(symLine)

    print(frameLine)

    numLine = contentLine.replace("|"," ")
    numLine = numLine.replace("#"," "*width)
    colNums = " ".join(f"{i:<{width}d}" for i in range(cols))
    numLine = numLine.replace("values",colNums)
    print(numLine)

    ## plot piar in build frame F: nodes x degree-factors ##

    # reshape #

    # print nodes of g #
    print(f'nodes of g: {" ".join(str(n) for n in gn)}')

    #print('nodes of g:',f"{*gn, }",'\n')
    #f'unpack a list: {" ".join(str(x) for x in a)}'
    #print("edges of g:"," ,".join(map("{:1}".format,*ge)),'\n')

    # uncomment me this prints edge tuples in lexographical #
    #print('\n'.join(" ".join(map("{:>3}".format, e)) for e in ge))

    ## loop ##
    print("Nu\n")
    print(Nu)
    print("Nv\n")
    print(Nv)
    print("dim")
    print(dim)

    #g matrix print#
    #for j in it.product

    ####cost matrix print#
    #block = bcmblock(int(n))
    #blockstring = as_string(block)
    #
    #print('nodes of g:',f"{*gn, }",'\n')
    #f'unpack a list: {" ".join(str(x) for x in a)}'
    #print("edges of g:"," ,".join(map("{:1}".format,*ge)),'\n')

    # uncomment me this prints edge tuples in lexographical #
    #print('\n'.join(" ".join(map("{:>3}".format, e)) for e in ge))

    ## loop ##
    print("Nu\n")
    print(Nu)
    print("Nv\n")
    print(Nv)
    print("dim")
    print(dim)

    #g matrix print#
    #for j in it.product

    ####cost matrix print#
    block = bcmblock(int(n))
    print(block)
    
    #
    #print("makeLR input")
    #print(block[1])
    Nu, Nv = makeLR(dim)
    #print(*Nu)
    #print(*Nv)
    #pulseWidth(Nu,Nv)
    #printlexV2V(g,Nu,Nv)
    print(" . "*250)
    printlexAll(g,Nu,Nv)

    #cost matrix build, put in networkx, solve using in-built#

    #solve using hungarian#

    #plot results using graph 2 ascii#

    #plot results on point-group lattice " "*n#


if __name__ == "__main__":
    main()

