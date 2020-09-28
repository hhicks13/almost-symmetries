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
import itertools as it
from scipy.spatial.distance import squareform
import numpy as np
import re
import math
import fractions
from numpy.random import permutation
from colorama import Fore
from colorama import Style
from colorama import Back
from operator import itemgetter

def matrix2coord(input_matrix,shape):
    grid = [(idx[0],idx[1],input_matrix[idx[0]][idx[1]]) for idx in it.product(*[range(s) for s in shape])]
    return grid

def asciiTable(input_list1,rowsize,colsize,theta):
    # stackoverflow recipe #
    rows = rowsize
    cols = colsize
    content = [["."]*cols for _ in range(rows)]

    grid =input_list1

    for (y,x,c) in grid: content[y][x] = c

    # build frame
    width       = len(str(max(rows,cols)-1))
    contentLine = "valuesvaluesvaluesvaluesvalues"

    dashes      = "-".join("-"*width for _ in range(cols))
    frameLine   = contentLine.replace("values",dashes)
    frameLine   = frameLine.replace("#"," "*width)
    frameLine   = frameLine.replace("| ","+-").replace(" |","-+")

    # print grid
    #print(frameLine)
    out = []
    for i,row in enumerate(reversed(content),1):
        values = "".join(f'{Fore.BLACK}{Back.BLUE} {v} {Style.RESET_ALL}' if v%3==0 else  f'{Fore.BLACK}{Back.CYAN} {v} {Style.RESET_ALL}'for v in it.islice(it.cycle(row),theta,theta+len(row)))
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

def main():
    n = math.comb(10,5)+1
    print(n)
    a = squareform(np.arange(1,11))
    dim = a.shape
    symran = a.tolist()
    print(dim)
    symrancoords = matrix2coord(symran,dim)

    i = 0
    while(i != -1):
        i = (i+1)%n
        asciiTable(symrancoords,dim[0],dim[1],i)
        time.sleep(0.1)
        asciiTable(symrancoords,dim[0],dim[1],i+1)
        time.sleep(0.1)
        asciiTable(symrancoords,dim[0],dim[1],i+2)
        time.sleep(0.1)
        asciiTable(symrancoords,dim[0],dim[1],i+3)
        
    
    
    
    


if __name__ == "__main__":
    main()
