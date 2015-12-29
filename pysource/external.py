# -*-Python-*-
################################################################################
#
# File:         external.py
# RCS:          $Header: $
# Description:  Routines for communicating with external yagm program.
# Author:       Staal Vinterbo
# Created:      Tue Oct 13 11:14:31 2015
# Modified:     Tue Dec 29 10:56:16 2015 (Staal Vinterbo) staal@klump.gateway.pace.com
# Language:     Python
# Package:      N/A
# Status:       Experimental
#
# (c) Copyright 2015, Staal Vinterbo, all rights reserved.
#
################################################################################

__all__ = ['external', 'write_adjmat']

from .graphmatch import GraphMatch, TreeMatch, hungarian, strength, randtree
import networkx as nx
from itertools import combinations, permutations
import sys
from logging import debug, info, warning
import time
from time import clock, time
import subprocess
import tempfile
import os

def write_adjmat(G, f, sep=',', header=True, sort=True):
    #print(list(w['weight'] for (a,b,w) in G.edges(data=True)))
    nl = sorted(G.nodes()) if sort else G.nodes()
    m = nx.to_numpy_matrix(G, nodelist=nl).getA()
    if header:
        f.write(sep.join(map(str, nl)) + '\n')
    if G.node[nl[0]].has_key('weight'):
        f.write(sep.join(str(G.node[n]['weight']) for n in nl) + '\n')
    else:
        f.write(sep.join(["1.0"] * len(nl)) + '\n')
    for row in m:
        f.write(sep.join(map(str, row)) + '\n')

def printawmatrix(G1, G2, similarity = True, to1=sys.stdout, to2=sys.stderr):

    if not similarity:
        write_adjmat(G1, to1, header=False)
        write_adjmat(G2, to2, header=False)
        return
        
    gm = GraphMatch(G1, G2)
    ediff = lambda h1, h2 : 0.0 if h1 == {} or h2 == {} else gm.w(h1, h2)
    ndiff = gm.v
    
    v1 = sorted(G1.nodes(data=True))
    v2 = sorted(G2.nodes(data=True))

    debug('v1: ' + ','.join(map(str, (a for (a,_) in v1))) + '\n')
    debug('v2: ' + ','.join(map(str, (a for (a,_) in v2))) + '\n')

    n1 = dict((n,i) for i,(n,_) in enumerate(v1))
    n2 = dict((n,i) for i,(n,_) in enumerate(v2))    

    for a,ha in v1:
        to1.write(','.join(map(str, (ndiff(ha,hb) for b,hb in v2))) + '\n')
                                   

    m,n = len(G1), len(G2)

    ifun = lambda a,b : (lambda (x,y) : (x*(x - 1)/2 + y))((a,b) if a > b else (b,a))
    ifun1 = lambda a,b : (lambda c : c - (c/(m + 1) + 1))(a * m + b)
    ifun2 = lambda a,b : (lambda c : c - (c/(n + 1) + 1))(a * n + b)

    f1 = ifun1 if nx.is_directed(G1) else ifun
    f2 = ifun2 if nx.is_directed(G2) else ifun

    geth = lambda G, a, b : G[a][b] if a in G and b in G[a] else {}

    e1 =[(a,b,geth(G1, a, b)) for (a,b) in (permutations(G1, 2) if nx.is_directed(G1) else combinations(G1, 2))]
    e1 = sorted(e1, key = lambda (a,b,c) : f1(n1[a],n1[b]))
    e2 =[(a,b,geth(G2, a, b)) for (a,b) in (permutations(G2, 2) if nx.is_directed(G2) else combinations(G2, 2))]
    e2 = sorted(e2, key = lambda (a,b,c) : f2(n2[a],n2[b]))

    #e1 = sorted(G1.edges(data=True), key = lambda (a,b,c) : f1(n1[a],n1[b]))
    #e2 = sorted(G2.edges(data=True), key = lambda (a,b,c) : f2(n2[a],n2[b]))


    debug('e1: ' + ','.join(map(str, ((i,(a,b,c)) for i, (a,b,c) in enumerate(e1)))) + '\n')
    debug('e2: ' + ','.join(map(str, ((i,(a,b,c))  for i, (a,b,c) in enumerate(e2)))) + '\n')    

    for (fa,ta,ha) in e1:
        to2.write(','.join(map(str, (ediff(ha,hb) for (fb,tb,hb) in e2))) + '\n')


###### listing of trees using a fold
# nx DiGraph interface
npreorder = lambda t : list(nx.dfs_preorder_nodes(t, min(t.in_degree().items(), key = lambda (x,y) : y)[0]))
getroot = lambda t: npreorder(t)[0]
maketree = lambda T: (T, getroot(T))
label = lambda (T, root) : root
subtrees = lambda (T, root) : map(lambda x : (T, x), T[root])
def treefold(f, g, a, tree): return f(label(tree), treefolds(f, g, a, subtrees(tree)))
def treefolds(f, g, a, strees):
    return a if strees == [] else g(treefold(f, g, a, strees[0]), treefolds(f,g,a, strees[1:]))
def listtree(T,trans):
    f = lambda label, subtrees : (trans[label], subtrees)
    g = lambda subtree, subtrees : [subtree] + subtrees
    z = []
    return treefold(f, g, z, maketree(T))
def mktrans(G):
    flip = lambda (x,y) : (y,x)
    return dict(map(flip, enumerate(sorted(G.nodes()))))

def external(G1, G2, trees, eps = 0.2, gabound = False, program = 'yagm', seed = None, verbose=False, similarity=False, strict=False):
    #f0 = '/tmp/colorinput985876.0.txt'
    #f1 = '/tmp/colorinput985876.1.txt'    
    #f2 = '/tmp/colorinput985876.2.txt'

    fd0, f0 = tempfile.mkstemp(text=True)
    fd1, f1 = tempfile.mkstemp(text=True)    
    fd2, f2 = tempfile.mkstemp(text=True)
    fd0 = os.fdopen(fd0, 'w')
    fd1 = os.fdopen(fd1, 'w')
    fd2 = os.fdopen(fd2, 'w')    
    
    #fd0 = open(f0, 'w')
    
    for tree in trees:
        thetree = tree()
        if len(thetree) > 0:
            fd0.write(str(listtree(thetree, mktrans(G1))) + '\n')
    fd0.close()
    #fd1 = open(f1, 'w')    
    #fd2 = open(f2, 'w')    
    printawmatrix(G1, G2, similarity = similarity, to1=fd1, to2=fd2)
    fd1.close()
    fd2.close()
    start = time()
    arglist = ([program, f0, f1, f2, "-e", str(eps)] +
               (["-v"] if verbose else []) +
               (["-g"] if gabound else []) +
               (["-strict"] if strict else []) +               
               (["-m"] if not similarity else []) +
               (["-seed", str(seed)] if seed != None else []))
    debug('calling external program: ' + str(arglist))
    output = subprocess.Popen(arglist,
                              stdout=subprocess.PIPE).communicate()[0]
    mat = map(lambda s : map(lambda s : float(s[:-1]), s.strip().split(',')),
              output.strip().split('\n'))
    hist = dict(zip(sorted(G1), (dict(zip(sorted(G2), row)) for row in mat)))
    end = time()

    # remove tempfiles
    os.unlink(f0)
    os.unlink(f1)
    os.unlink(f2)    
    
    return (hungarian(hist), strength(hist), end - start)


