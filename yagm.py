#!/usr/env python
################################################################################
#
# File:         yagm.py
# RCS:          $Header: $
# Description:  Yet Another Graph Matcher: a random subtree matching
#               based approach for the assignment problem.
# Author:       Staal Vinterbo
# Created:      Fri Jun 10 14:31:48 2011
# Modified:     Fri Apr 22 11:21:25 2016 (Staal Vinterbo) staal@klump.gateway.pace.com
# Language:     Python
# Package:      N/A
# Status:       Experimental
#
# (c) Copyright 2011, Staal Vinterbo, all rights reserved.
#
#
# yagm.py is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# yagm.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with yagm.py; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
################################################################################

from itertools import product, combinations, takewhile, islice, permutations
from yagma import GraphMatch, TreeMatch, gaalign, hungarian, strength, randtree
from yagma import external, write_adjmat
from yagma.graphmatch import dsim
import sys
import time
from time import clock, time
import re
from random import random, sample, seed, randint
from collections import defaultdict
import networkx as nx
import networkx.algorithms.isomorphism as iso
from logging import info, debug, error, warning
import os
import subprocess

import gc


def printproblem(G1, G2, to=sys.stdout):
    gm = GraphMatch(G1, G2)
    ediff = gm.w
    ndiff = gm.v
    
    v1 = sorted(G1)
    v2 = sorted(G2)

    vmap1 = dict((k, i) for i,k in enumerate(v1))
    vmap2 = dict((k, i) for i,k in enumerate(v2))

    #to.write('c ' + str(dict(enumerate(v1))) + '\n')
    #to.write('c ' + str(dict(enumerate(v2))) + '\n')    

    e1 = sorted(G1.edges())
    e2 = sorted(G2.edges())

    asses = list(enumerate(sorted(product(v1,v2))))
    # # edges between assignments in tensor product
    nedges = 0
    for (i,(v,w)),(j, (x,y)) in combinations(asses, 2):
        if v == x or w == y:
            continue
        nedges += 1

    # problem parameters
    al = len(asses)
    to.write(' '.join(map(str, ['p', len(v1), len(v2), al, nedges])) + '\n')

    # assignments (vertices in tensor product)
    for i, (v,w) in asses:
        to.write(' '.join(map(str, ['a', i, str(vmap1[v]), str(vmap2[w]),
                                    ndiff(G1.node[v], G2.node[w])])) + '\n')

    # edges between assignments in tensor product
    for (i,(v,w)),(j, (x,y)) in combinations(asses, 2):
        if v == x or w == y:
            continue
        to.write(' '.join(map(str, ['e', i, j,  ediff(G1[v][x], G2[w][y])])) + '\n')

    #for i,j in combinations(range(len(v1)), 2):
    #    to.write(' '.join(['n0', str(i), str(j)]) + '\n')
    #for i,j in combinations(range(len(v2)), 2):
    #    to.write(' '.join(['n1', str(i), str(j)]) + '\n')

    return dict(enumerate(v1)), dict(enumerate(v2))
    



class ubi:
    def __init__(self, T, G, url='http://127.0.0.1:20738/RPC2'):
        self.G = G
        self.T = T
        self.U = ubigraph.Ubigraph(URL=url)
        self.U.clear()
        self.Tedges = {}
        self.Tnodes = dict(
            (n, self.U.newVertex(shape="sphere", color="#008B00",
                        fontcolor="#ffff00", size=0.4,
                        label=str(n), fontsize=20)) for n in T)
        for a,b,d in T.edges(data=True):
            self.Tedges[(a,b)] = self.U.newEdge(self.Tnodes[a], self.Tnodes[b], label=str(round(d['weight'],3)),
                      strength=d['weight'], fontsize=10,
                      fontcolor ='#8B6508')
        self.Gedges = {}
        self.Gnodes = dict(
            (n, self.U.newVertex(shape="sphere", color="#ff0000",
                            size= 0.1,
                            label=str(n), fontsize=20)) for n in G)
        for a,b,d in G.edges(data=True):
            self.Gedges[(a,b)] = self.U.newEdge(self.Gnodes[a], self.Gnodes[b],
                  strength=d['weight'], fontsize=10,
                  fontcolor ='#8B6508')
        self.preva = []
        self.Gweighted = []
        self.prevb = []

    def __call__(self, ass, strength):
        # delete previous assignment
        for ua in self.preva:
            ua.destroy()
        self.preva = []            
        for e in self.Gweighted:
            e.set(label=None)
        self.Gweighted = []
        for b in self.prevb:
            b.set(size = 0.1)

        for a,b in self.T.edges():
            ga, gb = ass[a], ass[b]
            if (ga, gb) in self.Gedges and 'weight' in self.G[ga][gb]:
                self.Gedges[(ga, gb)].set(label=str(self.G[ga][gb]['weight']))
                self.Gweighted.append(self.Gedges[(ga, gb)])

        for f,t in ass.items():
            self.preva.append(self.U.newEdge(self.Tnodes[f], self.Gnodes[t], stroke="solid",
                      arrow_radius = 0.2 + strength[f], width=1.3,
                      color='#CAE1FF',
                      arrow=True))
            self.Gnodes[t].set(size=0.4)
            self.prevb.append(self.Gnodes[t])
        

        
        
            

        
        
        

        
        

def ubiviz(T, G, ass, ranks, url='http://127.0.0.1:20738/RPC2'):
    '''vizualize graphs and assignments using ubigraph'''

    U = ubigraph.Ubigraph(URL=url)
    U.clear()

    assrev = defaultdict(lambda : None, [(b,a) for a,b in ass.items()])

    Tnodes = dict(
        (n, U.newVertex(shape="sphere", color="#008B00",
                        fontcolor="#ffff00", size=0.4,
                        label=str(n), fontsize=20)) for n in T)

    for a,b,d in T.edges(data=True):
        U.newEdge(Tnodes[a], Tnodes[b], label=str(round(d['weight'],3)),
                  strength=d['weight'], fontsize=10,
                  fontcolor ='#8B6508')

    Gnodes = dict(
        (n, U.newVertex(shape="sphere", color="#ff0000",
                        size= 0.4 if n in assrev else 0.1,
                        label=str(n), fontsize=20)) for n in G)

    for a,b,d in G.edges(data=True):
        lab = None
        if T.has_edge(assrev[a], assrev[b]) or T.has_edge(assrev[b], assrev[a]):
            lab = str(round(d['weight'], 3))
        U.newEdge(Gnodes[a], Gnodes[b],
                  label = lab, strength=d['weight'], fontsize=10,
                  fontcolor ='#8B6508')

    for f,t in ass.items():
        U.newEdge(Tnodes[f], Gnodes[t], stroke="solid",
                  arrow_radius = 0.2 + ranks[f], width=1.3,
                  color='#CAE1FF',
                  arrow=True)

def randoms(G1, G2, check, times = 100):
    '''random search algorithm'''
    bv, ba = 0, None
    for i in range(times):
        a = dict(zip(G1, sample(G2, len(G1))))
        v = check(G1, G2, a)
        if v > bv:
            bv, ba = v, a
    return a


def graphmlread(f):
    '''read a graphml graph representation from f

    changes graphml attribute 'name' to become node id,
    and applies 'float()' to all edge attribute values.'''
    g = nx.read_graphml(f)
    return g                          # XXXX FIXME
    trans = dict(g.nodes(data=True))
    edges = []
    for a,b,d in g.edges(data=True):
        edges.append((trans[a]['name'],
                      trans[b]['name'],
                      dict(map(lambda (k,v): (k,float(v)), d.items()))))
    gg = type(g)()
    gg.add_edges_from(edges)
    # copy vertex data
    gg.add_nodes_from((w['name'],dict((k,v) for k,v in w.items() if k!='name'))
                      for n,w in trans.items())
    return gg

def adjmatread(fin, sep=' ,;:\t', gtype = nx.Graph):
    f = open(fin) if type(fin) == str else fin
    mat = [re.split('[ \t,;]+', line.strip()) for line in f]
    #debug('1')
    names = mat[0]
    #debug(str((len(mat), len(names), mat, names)))
    vweights = len(mat) > len(names) + 1
    start = 1
    vw = defaultdict(lambda : 1)
    if vweights:
        #debug('2')
        start += 1
        vw.update(zip(names, map(float, mat[1])))
        #debug('3')        
    try:
        #debug('4')
        amat = [map(float, row) for row in mat[start:]]
        #debug('5')        
    except:
        error('Non float value in adjacency matrix!')
        sys.exit(1)
    G = gtype()
    #debug('6')            
    G.add_nodes_from((n, {'weight': vw[n] }) for n in names)
    #debug('7')            
    G.add_edges_from((n1, n2, {'weight': amat[i][j] }) for
                     (i, n1), (j, n2) in product(enumerate(names),
                                                 enumerate(names))
                     if i != j)
    #debug('Read graph...')
    return G
    
    

def taufilter(g, tau):
    '''remove edges with weight < tau'''
    gg = g.copy()
    edges = filter(lambda (a,b,w) : w['weight'] < tau, gg.edges(data=True))
    gg.remove_edges_from((a,b) for (a,b,_) in edges)
    rnodes = filter(lambda n : gg.degree(n) == 0, gg.nodes())
    gg.remove_nodes_from(rnodes)
    return gg

 
# test code below

def gencol(G, T):
    '''generate coloring that yields optimal solution for
       known subgraph T of G'''
    a = [0]*len(G)
    for i,t in enumerate(T):
        a[t] = i
    return a

def genhist(v):
    '''generate a histogram of values in a dict'''
    h = defaultdict(lambda : 0)
    for w in v:
        h[w] += 1
    w = h.values()
    ws = sum(w)
    wh = reduce(max, w)
    wl = reduce(min, w)
    info('edge weight histogram #cat/min/max/ave: ' + str((len(w), wl, wh, ws/len(w))))
    return h

def gentest(n = 10, m = 5, p=0.5, randw = False, prec=100, directed=True):
    '''generate test graphs'''
    #G = nx.complete_graph(n) if p == 1 else nx.gnp_random_graph(n, p)
    G = nx.complete_graph(n) if p == 1 else nx.powerlaw_cluster_graph(n, 3, p)
    debug('generated graph...\n')
    for a,b in G.edges():
        weight = round(prec*random())/prec + 1.0/prec
        debug('setting edge weight ' + str(weight))
        G[a][b]['weight'] = weight


    genhist([w['weight'] for _,_,w in G.edges(data=True)])


    root = sample(G.nodes(), 1)[0]
    T = randtree(G, m, root)

    if randw:
        for (a,b) in T.edges():
            T[a][b]['weight'] += random()/5
            
    debug('T: ' + str(T.edges(data=True)))
    genhist([w['weight'] for _,_,w in T.edges(data=True)])
    info('gentest: G: ' + str((len(G), len(G.edges()))) + ',  T: ' + str((len(T), len(T.edges()))))
    genhist([w['weight'] for _,_,w in G.edges(data=True)])        

    for n in G.nodes():
        G.node[n]['weight'] = 1

    return(G, T, gencol(G,T))


def gentest2(n=20, p=0.1, s = 1):
    G = nx.gnm_random_graph(n, round((n**2)*p))
    rem = [n for n in G if G.degree(n) == 0]
    info('removing ' + str(rem))    
    G.remove_nodes_from(rem)
    for a,b in G.edges():
        G[a][b]['weight'] = random()
    G2 = G.copy()
    for a,b in G2.edges():
        G[a][b]['weight'] += uniform(0,s)
    return G2, G

def nsim(h1, h2):
    '''node value histogram compatibility

       assumes monotonicity in node information and graph size.
    '''
    lh1 = len(h1)
    lh2 = len(h2)
    if lh1 + lh2 == 0: return 1 # no information
    if lh1 > lh2: return 0 # smaller graph has more info on node
    return dsim(h1, h2)


############### Experiments

srank = lambda a,s : len(list(takewhile(lambda (_,i,j) : i == j,
                                                sorted(((s[i],i,j) for i,j in a.items()),
                                                       reverse=True))))
correct = lambda h : len([(a,b) for a,b in h.items() if a == b])

class Experiments:
    def __init__(self, seedvalue = 1, tests=0):
        self.correct = lambda h : len([(a,b) for a,b in h if a == b])
        self.gares = []
        self.tgres = []
        self.bigres = []
        self.seed = seedvalue
        self.tests = tests
        self.srank = lambda a,s : len(list(takewhile(lambda (_,i,j) : i == j,
                                                sorted(((s[i],i,j) for i,j in a),
                                                       reverse=True))))
        seed(seedvalue)

    def nextseed(self):
        seed(seed)
        self.seed = randint(-sys.maxint, sys.maxint)
        
    def genGG(self, m, n, prec=20):
        G2, T, cols = gentest(n, m, prec=prec)
        G1 = nx.subgraph(G2, sample(G2.nodes(), m)).copy()
        over005 = 0
        over0025 = 0
        if self.tests > 0: # add noise to G1 edges
            for _,_,w in G1.edges(data=True):
                noise = npr.normal(scale=self.tests)
                if abs(noise) > 0.05:
                    over005 += 1
                if abs(noise) > 0.025:
                    over0025 += 1
                w['weight'] += noise
            info('Edge Noise e: |e| > 0.025 : ' + str(over0025) + ', |e| > 0.05: ' + str(over005))

        return (G1,G2)
    
    def genTG(self, m, n, prec=20): 
        G, T, cols = gentest(n, m, prec=prec)
        return (T,G)
    
    def gaonly(self, reps = 10, ms = [5,10,20], ns = [20,30], uiw=True):
        for m, n in product(ms, ns):
            for i in range(reps):
                #self.nextseed()
                G1, G2 = self.genGG(m, n)
                start = clock()
                s, a = gaalign(G1, G2, GraphMatch(G1, G2).check)
                end = clock()
                self.printrow((m,n, end - start, s , self.correct(a.items())))


    def tg(self, reps = 10, ms = [5,9], ns = [15,40], uiw=True):
        for m, n in product(ms, ns):
            gc.collect()
            debug('TG: garbage length: ' + str(len(gc.garbage)))
            for i in range(reps):
                info('TG: ' + str((m,n,i)))
                #self.nextseed()
                T, G = self.genTG(m, n)
                tm1 = TreeMatch(T, G)
                tm2 = TreeMatch(T, G)                
                gm = GraphMatch(T,G)
                start = clock()
                sbound,_ = gaalign(T, G, gm.check)
                end = clock()
                tsb = end - start
                start = clock()
                anb,snb = tm1(uisim=uiw)
                end = clock()
                #sys.stderr.write(str((m,n, anb, snb)) + '\n')
                tnb = end - start
                start = clock()
                ab,sb = tm2(sbound=sbound - 1e-10, uisim=uiw)
                end = clock()
                #sys.stderr.write(str((m,n,ab,sb)) + '\n')
                tb = end - start
                self.printrow((m,n, tnb, snb, self.correct(anb[0]) if anb != None else 0,
                                   tb, sb, self.correct(ab[0]) if ab != None else 0, tsb, sbound))
    
    def big(self, reps = 5, ms = [10,20], ns=[20,40], ks = [4,6], cs=[2,5], uiw=True, save=None):
        for m, n, k, coverage in product(ms, ns, ks, cs):
            info('Experiment(big): ' + str((m,n,k,coverage)))
            gc.collect()
            debug('BIG: garbage length: ' + str(len(gc.garbage)))
            for i in range(reps):
                #self.nextseed()
                G1, G2 = self.genGG(m, n)
                gm = GraphMatch(G1, G2, k = k)
                start = clock()
                a,s = gm(coverage=coverage, uisim=uiw)
                end = clock()
                meta = (m,n, k, coverage, end - start, self.srank(a.items(),s), self.correct(a.items()), gm.check(G1, G2, a))
                self.printrow(meta)
                if save != None:
                    fnprefix = save + '/' + '-'.join(map(str, (m,n,k,coverage,i)))
                    f = open(fnprefix + '-G1.txt', 'w')
                    write_adjmat(G1, f)
                    f.close()
                    f = open(fnprefix + '-G2.txt', 'w')
                    write_adjmat(G2, f)
                    f.close()
                    f = open(fnprefix + '-META.txt', 'w')
                    f.write(','.join(map(str, meta)) + '\n' + str(a)+ '\n' + str(s) + '\n')
                    f.close()
                    

    def printrow(self, row):
        print(' '.join(map(str, row)))

    

class reporter:
    def __init__(self, T, G):
        self.T = T
        self.G = G
        self.check = GraphMatch(T,G).check
        self.score = lambda ass : self.check(self.T, self.G, ass)
        self.start = clock()
    def __call__(self, T, flist, s, cover, B):
        ass = hungarian(B)
        ranks = strength(B)
        print(' '.join(map(str,
                           [len(T), s, self.score(ass), correct(ass), srank(ass, ranks),
                            clock() - self.start, len(flist)])))



def readfrom(fin):
    f = fin
    if type(fin) == str:
        f = sys.stdin if fin == '-' else open(fin)
    ass = dict()
    for line in f:
        pairs = re.split('[,]+', line.strip())
        for pair in pairs:
            #debug(pair)
            s = re.split('[:\->]+', pair.strip())
            if s == None or len(s) != 2:
                continue
            #debug('s: ' + str(s))
            #try:
            #    a,b = int(s[0].strip()), int(s[1].strip())
            #    debug('inttype')
            #except:
            #    debug('stringtype')
            a,b = s[0].strip(), s[1].strip()
            #debug('a,b = ' + str((a,b)))
            ass[a] = b
    return ass

    


def degrees(G):
    degs = G.degree().values()
    return (min(degs), sum(degs)/float(len(degs)), max(degs))

def update(G1, G2, pair):
    '''remove matched pair from pair of graphs'''
    a, b = pair
    for n in G1.neighbors(a):
        G1.node[n]['weight'] += G1[a][n]['weight']
    for n in G2.neighbors(b):
        G2.node[n]['weight'] += G2[a][n]['weight']

    G1.remove_node(a)
    G2.remove_node(b)

    return (G1, G2)
        

###############


if __name__ == "__main__":

    import logging
    from logging import debug, warning, error, info

    import networkx as nx
    from random import sample, seed
    from collections import defaultdict
    import yagma.ubigraph as ubigraph

    from yagma import hungarian, strength, treeiterator, randtree, randstar
    from yagma import vf2match, gridvf2

    import numpy.random as npr

    formats = {'graphml': graphmlread,
               'pajek': nx.read_pajek,
               #'dot': nx.read_dot,
               'edgelist' : nx.read_weighted_edgelist,
               'amat' : adjmatread}

    wformats = {'graphml': nx.write_graphml,
               'pajek': nx.write_pajek,
               #'dot': nx.write_dot,
               'edgelist' : nx.write_weighted_edgelist,
                'amat' : write_adjmat}

    from optparse import OptionParser
    Version = 2.71
    parser = OptionParser(usage = "".join(
        ["%prog [options] URL1 URL2\n",
         "\n Version: ", str(Version), " SAV (C) 2011\n",
         " Compute a vertex correspondence between two graphs.\n",
         " The graphs can have numeric vertex and edge attributes in [0,1].\n",
         " At least one edge attribute must be 'weight'.\n",
         " URLS can also be a simple filenames or '-' for stdin.\n"]),
                          version = ''.join(["%prog ",str(Version)," (C) SAV 2007"]))
    parser.add_option("-I", "--inputformat", dest="iformat", metavar="NAME",
		      default = 'graphml',
		      choices = formats.keys(),
                      help="graph input/output format. " + "Choice: " + str(formats.keys()) + " Default: %default")
    parser.add_option("-O", "--outputformat", dest="oformat", metavar="NAME",
		      default = 'graphml',
		      choices = formats.keys(),
                      help="graph input/output format. " + "Choice: " + str(formats.keys()) + " Default: %default")
    parser.add_option("-v", "--verbose", dest="verbose", action="count",
                       default = 0,
                       help="print more informational "
                      "stuff than strictly needed. Given"
                       " twice results in debug info.")
    parser.add_option('-k', "--subsize", dest="k", metavar="INT",
		      default = 5,
		      type = 'int',
                      help="Random subtree size (default: %default)")
    parser.add_option("--seed", dest="seed", metavar="INT",
		      default = None,
		      type = 'int',
                      help="Set random seed (default: %default)")
    parser.add_option("--eps", dest="eps", metavar="FLOAT",
		      default = 0.2,
		      type = 'float',
                      help="color coding target probability (default: %default)")
    parser.add_option("--tau", dest="tau", metavar="FLOAT",
		      default = None,
		      type = 'float',
                      help="filter edges with weight smaller than (default: %default)")
    parser.add_option("--maxit", dest="maxit", metavar="INT",
		      default = sys.maxint,
		      type = 'int',
                      help="Maximal number of iterations to wait for improvement (default: %default)")
    parser.add_option("--coverage", dest="coverage", metavar="INT",
		      default = 2,
		      type = 'int',
                      help="Coverage of each node required (default: %default)")
    parser.add_option("--test", dest="test",
		      default = False,
                      action='store_true',
                      help="generate random test graphs.")
    parser.add_option("--testtree", dest="testtree",
		      default = False,
                      action='store_true',
                      help="generate random test graphs, but let G1 be a tree (to me matched with color coding).")
    parser.add_option("--testN", dest="testn", metavar="INT",
		      default = 20,
		      type = 'int',
                      help="size of complete test graph (default: %default)")
    parser.add_option("--testM", dest="testm", metavar="INT",
		      default = 10,
		      type = 'int',
                      help="size of random test subgraph (default: %default)")
    parser.add_option("--tests", dest="tests", metavar="INT",
		      default = 0,
		      type = 'float',
                      help="the std. dev. of  0-mean normal noise added to edges in G1  (default: %default)")
    parser.add_option("--testr", dest="testr", metavar="B",
		      default = 20,
		      type = 'int',
                      help="weights rounded into B + 1 bins (default: %default)")
    parser.add_option("--testp", dest="testp", metavar="Probability",
		      default = 1,
		      type = 'float',
                      help="probability of adding an edge in test graph (default: %default)")
    parser.add_option("--random", dest="random",
                          default = False,
                          action='store_true',
                          help="do random search instead with -k NUM trials.")
    parser.add_option("--external", dest="external",
                          default = False,
                          action='store_true',
                          help="do random search instead with -k NUM trials.")
    parser.add_option("--external_adjacency", dest="adjacency",
                          default = False,
                          action='store_true',
                          help="Output adjacency matrices (weight attribute only) to external program.")
    parser.add_option("--external_program", dest="program",
                          default = "yagm",
                          help="name of external executable (default: %default).")
    parser.add_option("--ga", dest="ga",
                      default = False,
                      action='store_true',
                      help="run ga on full instance.")
    parser.add_option("--nostar", dest="nostar",
                      default = False,
                      action='store_true',
                      help="don't restrict random trees to stars.")
    parser.add_option("--nostrict", dest="single",
                      default = True,
                      action='store_false',
                      help="have treematch filter less aggressively.")
    parser.add_option("--unit", dest="unit",
                      default = False,
                      action='store_true',
                      help="Assume unit interval similarities to determine early stopping.")
    parser.add_option("--gabound", dest="gabound",
                      default = False,
                      action='store_true',
                      help="run ga to establish lower bounds.")
    parser.add_option("--nocc", dest="nocc",
                      default = False,
                      action='store_true',
                      help="match trees using ga only.")
    parser.add_option("--xout", dest="xout",
                      default = False,
                      action='store_true',
                      help="output one line of experimental stats.")
    parser.add_option("--xoutonly", dest="xoutonly",
                      default = False,
                      action='store_true',
                      help="only print --xout line.")
    parser.add_option("--vf2", dest="vf2",
                      default = False,
                      action='store_true',
                      help="Run VF2 algorithm (works on (subgraph) isomorphisms only).")
    parser.add_option("--grid", dest="grid",
                      default = False,
                      action='store_true',
                      help="Run VF2 algorithm with grid search on parameters.")
    parser.add_option("--gridn", dest="gridn", metavar="INT",
		      default = 10,
		      type = 'int',
                      help="use INT * INT grid for vf2 grid (default: %default)")
 
    parser.add_option("--vf2a", dest="vf2a",
                      metavar = 'FLOAT',
                      default = 1.0e-9,
                      type = 'float',
                      help="VF2 algorithm absolute tolerance.")
    parser.add_option("--vf2r", dest="vf2r",
                      metavar = 'FLOAT',
                      default = 10.0e-7,
                      type = 'float',
                      help="VF2 algorithm relative tolerance.")
    parser.add_option("-e", "--morehelp", dest="morehelp",
                      default = False,
                      action='store_true',
                      help="print more help.")
    parser.add_option("--printproblem", dest="printproblem",
                          default = False,
                          action='store_true',
                          help="Print problem spec in terms of"
                      " tensor product graph (to output1).")
    parser.add_option("--printawmatrix", dest="printawmatrix",
                          default = False,
                          action='store_true',
                      help="Print A and W graphs (to output1 and output2).")
    parser.add_option("--output1", dest="output1",
                          default = None,
                          help="Filename to print to (default: stdout)")
    parser.add_option("--output2", dest="output2",
                          default = None,
                          help="Filename to print to (default: stderr)")
    parser.add_option("--printgraphs", dest="printgraphs",
                          default = False,
                          action='store_true',
                          help="Print graphs (to output1 and output2).")
    parser.add_option("--ubigraph", dest="ubigraph",
                      default = False,
                      action='store_true',
                      help="Visualize assignment using ubigraph.")
    parser.add_option("--animate", dest="animate",
                      default = False,
                      action='store_true',
                      help="'animate' algorithm. With --ubigraph,"
                      " sends updates to the server, otherwise prints --xout type lines.")
    parser.add_option("--ubiurl", dest="ubiurl", metavar="URL",
		      default = 'http://127.0.0.1:20738/RPC2',
                      help="RPC URL of ubigraph server (default: %default).")
    parser.add_option("--experiment", dest="experiment", metavar="INT",
		      default = None,
		      type = 'int',
                      help="run experiment number INT (default: %default)")
    parser.add_option("-R", "--readfrom", dest="readfrom", metavar="URL",
		      default = None,
                      help="Read assignment from input instead of computing."
                      " URL can be a url, simple filename, or '-' for stdin."
                      " assignment format can be pairs a(:|[ ]*->[ ]*)b. When"
                      " more than one on each line, the must be separated by ','.")
    

    (opt, args) = parser.parse_args()
    logging.basicConfig(level=max(logging.DEBUG,
                                  logging.WARNING - 10*opt.verbose),
                        format='%(levelname)-5s %(message)s')

    if opt.morehelp:
        import yagma.docstring as docstring
        print(docstring.manpage)
        sys.exit(0)

    if opt.experiment != None:
        E = Experiments(opt.seed if opt.seed != None else 1, tests=opt.tests)
        E2 = Experiments(opt.seed if opt.seed != None else 1, tests=opt.tests)        
        el = [E.gaonly, E.tg, E.big]
        try:
            el[opt.experiment](uiw=opt.unit)
        except Exception as e:
            error('Illigal experiment number.' + str(e))
            sys.exit(1)
        sys.exit(0)

    if len(args) < 2 and not opt.test:
        error("Missing url/files. Try -h for help.")
        sys.exit(1)

    if opt.seed != None:
        seed(opt.seed)
        npr.seed(opt.seed)

    gabound = lambda a,b,c : (0,None)
    if (opt.gabound) or opt.nocc:
	gabound = gaalign

    if opt.testtree:
        info('--testtree is set, also setting --strict and --test')
        opt.test = True
        opt.single = True

    if not opt.test:
        # read input graphs
        info('fetching graphs from ' + args[0] + ', ' + args[1])
        try:
            G1 = formats[opt.iformat](args[0])
            G2 = formats[opt.iformat](args[1])
            debug('G1.edges:' + str(G1.edges(data=True)))
        except Exception, e:
            error('could not read graphs!')
            error(str(e))
            sys.exit(1)
    else:
        G2, T, cols = gentest(opt.testn, opt.testm, opt.testp, prec=opt.testr)
        if opt.testtree:
            G1 = T
        else:
            G1 = nx.subgraph(G2, sample(G2.nodes(), opt.testm)).copy()
        over005 = 0
        over0025 = 0
        if opt.tests > 0: # add noise to G2 edges
            for _,_,w in G1.edges(data=True):
                noise = npr.normal(scale=opt.tests)
                if abs(noise) > 0.05:
                    over005 += 1
                if abs(noise) > 0.025:
                    over0025 += 1
                w['weight'] += noise
            info('Edge Noise e: |e| > 0.025 : ' + str(over0025) + ', |e| > 0.05: ' + str(over005))
            

    info('G1 # edges: ' + str(len(G1.edges())))
    info('G1 average clustering: ' + str(nx.average_clustering(G1.to_undirected())))
    info('G1 degree min/ave/max: ' + str(degrees(G1)))
    info('G1 # nodes: ' + str(len(G1.nodes())))    
    info('G2 # edges: ' + str(len(G2.edges())))
    info('G2 average clustering: ' + str(nx.average_clustering(G2.to_undirected())))
    info('G2 degree min/ave/max: ' + str(degrees(G2)))    
    info('G2 # nodes: ' + str(len(G2.nodes())))    

    G1f = G1
    G2f = G2
    if(opt.tau != None):
        G1f = taufilter(G1, opt.tau)
        G2f = taufilter(G2, opt.tau)        
        info('filtered G1 # edges: ' + str(len(G1f.edges())))
        info('filtered G1 # nodes: ' + str(len(G1f)))
        info('filtered G2 # edges: ' + str(len(G2f.edges())))    
        info('filtered G2 # nodes: ' + str(len(G2f)))

    if opt.printgraphs:
        wformats[opt.oformat](G1f, sys.stdout if opt.output1 == None else open(opt.output1, 'w'))
        wformats[opt.oformat](G2f, sys.stderr if opt.output2 == None else open(opt.output2, 'w'))        
        sys.exit(0)
        
    if opt.printproblem:
        printproblem(G1f, G2f, sys.stdout if opt.output1 == None else open(opt.output1, 'w'))
        sys.exit(0)

    if opt.printawmatrix:
        printawmatrix(G1f, G2f,
                         sys.stdout if opt.output1 == None else open(opt.output1, 'w'),
                         sys.stderr if opt.output2 == None else open(opt.output2, 'w'))
        sys.exit(0)

    match = GraphMatch(G1f, G2f,
                       quickbound=gabound, v=nsim,
                       k=opt.k, eps=opt.eps)

    addsec = 0
    stampstart = clock()
    if opt.readfrom != None:
        ass = readfrom(opt.readfrom)
        ranks = defaultdict(lambda : 0.5)
    elif opt.vf2:
        if not opt.grid:
            val, ass = vf2match(G1f, G2f, atol=opt.vf2a, rtol=opt.vf2r)
        else:
            val, ass = gridvf2(G1f, G2f, opt.gridn)
        ranks = defaultdict(lambda : 0.5)        
    elif opt.random:
        ass = randoms(G1f, G2f, match.check, opt.k)
        ranks = defaultdict(lambda : 0.5)
    elif opt.ga:
        val, ass = gaalign(G1f, G2f, match.check)
        ranks = defaultdict(lambda : 0.5)
    elif opt.external:
        if not opt.testtree:
            ntrees = opt.coverage * len(G1f)
            info('Generating ' + str(ntrees) + ' trees')
            trees = islice(treeiterator(G1f, opt.k, randtree if opt.nostar else randstar), ntrees)
        else:
            trees = [lambda : G1f]
        ass, ranks, add = external(G1f.to_undirected(), G2f, trees,
                                eps = opt.eps,
                                program=opt.program,
                                gabound=opt.gabound,
                                verbose = opt.verbose > 0,
                                seed = opt.seed, similarity=not opt.adjacency,
                                strict=opt.single)
        addsec = add
    elif opt.testtree:
        ass, val = TreeMatch(G1f, G2f, eps=opt.eps, sbound=0, v=nsim, w=dsim)(uisim=False)
        ranks = defaultdict(lambda : 0.5)
        ass = defaultdict(lambda : 0, ass[0])
    else:
        callback = None
        if opt.animate:
            if opt.ubigraph:
                uv = ubi(G1, G2, opt.ubiurl)
                callback = lambda T, flist, s, cover, B: uv(hungarian(B), strength(B))
            else:
                callback = reporter(G1, G2)
        ass, ranks = match(coverage=opt.coverage, nocc=opt.nocc, maxit=opt.maxit,
                           star=not opt.nostar, single=opt.single, uisim=opt.unit,
                           callback=callback)

        if False:
            #####
            bass = dict()
            branks = dict()
            mval = -1

            for ii in xrange(len(G1f) - 2):
                ass, ranks = match(coverage=opt.coverage, nocc=opt.nocc, maxit=opt.maxit,
                                   star=not opt.nostar, single=opt.single, uisim=opt.unit,
                                   callback=callback)


                info('# name correspondences: ' + str(sum( x == y for (x,y) in ass.items())))
                ncorrscore = match.check(G1, G2, dict((a,a) for a in G1))
                info('name correspondence score: ' +
                     str(ncorrscore))
                debug('G1nodes: ' + str(G1.nodes()))
                debug('G2nodes: ' + str(G2.nodes()))    
                tmp = match.check(G1, G2, ass)
                info('Sum (1 - |overlap edge weight diff|): ' +
                     str(match.check(G1, G2, ass)) + ' (# edges: ' +
                     str(len(G1.edges())) + ')')
                tkeys = sorted(ass.keys(), key=lambda x : ranks[x], reverse=True)

                bass[tkeys[0]] = ass[tkeys[0]]
                branks[tkeys[0]] = ranks[tkeys[0]]

                for k in tkeys:
                    print(str(k) + ' ' + str(ass[k]) + ' ' + str(ranks[k]))

                print(' ASSIGNED:' + str(tkeys[0]) + ' --> ' + str(ass[tkeys[0]]))
                print('---')

                G1f, G2f = update(G1f, G2f, (tkeys[0], ass[tkeys[0]]))

            bass[tkeys[1]] = ass[tkeys[1]]
            branks[tkeys[1]] = ranks[tkeys[1]]

            ass = bass
            ranks = branks
            ####
    

        
        
    stampstop = clock() + addsec

    info(' ranks: ' + str(sorted(ranks.items(), key = lambda (x,y): y)))
    info('length of assignment: ' + str(len(ass)))
    debug('assignment: ' + str(ass))
    info('# name correspondences: ' + str(sum( x == y for (x,y) in ass.items())))
    ncorrscore = match.check(G1, G2, dict((a,a) for a in G1))
    info('name correspondence score: ' +
         str(ncorrscore))
    debug('G1nodes: ' + str(G1.nodes()))
    debug('G2nodes: ' + str(G2.nodes()))    
    tmp = match.check(G1, G2, ass)
    info('Sum (1 - |edge diff|) + Sum(1 - |node diff|): ' +
         str(match.check(G1, G2, ass)) + ' (# edges: ' +
         str(len(G1.edges())) + ', '+ str(len(G2.edges())) + ')')

    if opt.xout or opt.xoutonly:
        print(','.join(map(str, [len(G1), len(G2), opt.k, opt.coverage, stampstop - stampstart,
                                 srank(ass, ranks), correct(ass), match.check(G1, G2, ass),
                                 ncorrscore])))
        if opt.xoutonly:
            sys.exit(0)

    tkeys = sorted(ass.keys(), key=lambda x : ranks[x], reverse=True)
    for k in tkeys:
        print(str(k) + ' ' + str(ass[k]) + ' ' + str(ranks[k]))

    if opt.ubigraph:
        try:
            maxr = reduce(max, ranks.values(), 0)
            minr = reduce(min, ranks.values(), 1)
            if minr != maxr:
                ranks = defaultdict(lambda : 1,
                                [(k, (v - minr)/(maxr - minr)) for
                                 k,v in ranks.items()])
            if not opt.animate:
                uv = ubi(G1,G2,opt.ubiurl)
            uv(ass, ranks)
        except Exception as e:
            error('Vizualization failed! ' + str(type(e)))
            error(str(e))
