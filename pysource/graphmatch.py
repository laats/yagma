# -*-Python-*-
################################################################################
#
# File:         graphmatch.py
# RCS:          $Header: $
# Description:  Color Coding based graph matching in simple attributed graphs
# Author:       Staal Vinterbo
# Created:      Sun Jul 17 15:46:11 2011
# Modified:     Mon Dec 19 17:36:16 2011 (Staal Vinterbo) staal@dink
# Language:     Python
# Package:      yagma
# Status:       Experimental
#
# graphmatch.py is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# graphmatch.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with graphmatch.py; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
# (c) Copyright 2011, Staal Vinterbo, all rights reserved.
#
################################################################################
#
# Revisions:
#
# Mon Dec 19 16:07:14 2011 (Staal Vinterbo) staal@dink
#  Fixed random generation of trees and stars
# Mon Sep 26 12:23:13 2011 (Staal Vinterbo) staal@dink
#  use the next root as the prefix end (for storage table)
################################################################################
'''graphmatch: color coding based matching in simple attributed graphs.'''

__all__ = ['GraphMatch', 'TreeMatch', 'hungarian', 'strength', 'randtree', 'randstar', 'treeiterator']

from random import sample, shuffle, randint, random, choice
import math
import sys
from itertools import product, repeat, chain, islice, ifilter, dropwhile, cycle, imap, takewhile, combinations
import networkx as nx
from collections import defaultdict
from operator import or_, add
from logging import debug, info, warning, error
from munkres import Munkres

_tolerance = 2*sys.float_info.epsilon

def mad(x,y):
    '''mean abs difference'''
    xy = zip(x,y)
    return float(sum(abs(a - b) for a,b in xy))/len(xy)

def dsim(h1,h2):
    '''1-mad for two numeric dicts'''
    return (1-mad([v for _,v in sorted(h1.items())], [v for _,v in sorted(h2.items())])
            if min(len(h1), len(h2)) > 0 else 1)


def findroot(T):
    debug('findroot ' + str(T.edges()))
    assert(len(T) > 0)
    l = sorted(T.in_degree_iter(), key=lambda (n,d): d)
    (root, deg) = l[0]
    assert(deg == 0)
    assert(len(l) > 1 and l[1][1] > 0)
    debug('findroot -> ' + str((root, deg)))
    return root

# ---------------------------------

def treeheight(tree, root):
    acc = []
    for child in tree[root]:
        acc += treeheight(tree, child)
    return [(root, len(acc) + 1)] + acc

def preordernodes(tree, root):
    acc = [root]
    for st in tree[root]:
        acc += preordernodes(tree, st)
    return acc

def randtree(G, k, root):
    '''return a random order k subtree of G rooted at root'''
    if k < 1: return nx.DiGraph()
    B,S,V = set((root, x) for x in G[root]), [], set()
    for i in range(k-1):
        if len(B) == 0: break
        (p,c) = choice(list(B))
        V |= set([p,c])
        B = set((a,b) for (a,b) in B if b not in V) 
        B |= set([(c,x) for x in G[c] if x not in V])
        S += [(p,c)]
    return nx.DiGraph((a,b,G[a][b]) for (a,b) in S)

def randstar(G, k, root):
    '''return a random order k minimum depth tree from G rooted at root.'''
    if k < 1: return nx.DiGraph()
    B,S,V = [(root, x) for x in G[root]], [], set()
    for i in range(k-1):
        if len(B) == 0:
            B = reduce(add, [[(c, x) for x in G[c] if x not in V]
                       for (p,c) in S], [])
            if len(B) == 0:
                break
        (p,c) = choice(list(B))
        V |= set([p,c])
        S += [(p,c)]
        B = set((a,b) for (a,b) in B if b not in V)
    return nx.DiGraph((a,b,G[a][b]) for (a,b) in S)

def treeiterator(G, k, method=randstar):
    return cycle(imap(lambda x : lambda : method(G,k,x), G))


class TreeMatch:
    def __init__(self, T, G, eps=0.1, sbound=0, v=dsim, w=dsim):
        self.T, self.G, self.eps, self.v, self.w = T, G, eps, v, w
        self.n, self.m, self.sbound = len(G), len(T), sbound
        if not nx.is_directed(T):
            error('TreeMatch: T is not directed!')
            raise ValueError('Input Tree T is not a directed graph!')
        if self.m >= math.log(sys.maxint, 2):
            error('TreeMatch: T is too large, colorset storage exceeded!')
            raise ValueError('Input Tree T too large!')

        if len(T) < 1:
            error('TreeMatch: T has no vertices!')
            raise ValueError('Input Tree T empty!')
        

        self.root = findroot(T)
        self.norder = preordernodes(self.T, self.root)
        debug('TreeMatchLF: norder = ' + str(self.norder))        
        nind = dict(zip(self.norder, range(self.m)))
        debug('TreeMatchLF: nind = ' + str(nind))                
        self.pre = [0] * (self.m + 1)
        for (a,b) in T.edges():
            self.pre[nind[b]] = nind[a]
        debug('TreeMatchLF: pre = ' + str(self.pre))
        self.times = int(math.ceil(math.exp(self.m) * math.log(1/self.eps)))
        self.color = self.randcolor()
        # filtering info
        vvals = sorted((v(T.node[a], G.node[b]) for a,b in product(T, G)), reverse=True)[:self.m]
        wvals = sorted((w(te, ge) for te, ge in product((d for _,_,d in T.edges(data=True)),
                                                        (d for _,_,d in G.edges(data=True)))),
                       reverse=True)[:self.m]
        svals = [a + b for a,b in zip(vvals, wvals)]
        maxes,_ = reduce(lambda (l, s), x : (l + [s + x], s + x), svals, ([0], 0)) # cumsum
        self.ibound = list(islice(chain(maxes, repeat(2)), 0, self.m))
        debug('TreeMatchLF.ibound = ' + str(self.ibound))
        self.table = None # defaultdict(lambda : [(0,[])] * (2**self.m - 1))
        self.keys = set()
        self.pruned = 0
        self.keep = self.keep_strict
        
    def keep_multi(self, val, i):
        return val + self.ibound[self.m - (i + 1)] >= self.sbound

    def keep_strict(self, val, i):
        return val + self.ibound[self.m - (i + 1)] > self.sbound
    
    def randcolor(self):
        tmp = range(self.m) + [randint(0, self.m - 1) for i in range(self.n - self.m)]
        shuffle(tmp)
        return dict(zip(self.G, tmp))

    def __call__(self, sbound=None,  maxit=sys.maxint, uisim = True, single = True):
        self.bestkey, self.bestval, self.sbound = None, 0, self.sbound if sbound == None else sbound
        self.keep = self.keep_strict if single else self.keep_multi
        
        def init():
            self.color = self.randcolor()
            self.bestkey = None
            self.bestval = 0
            self.table = defaultdict(lambda : [(0,[])] * (2**self.m - 1))
            self.keys = set()
            self.pruned = 0
            for u in self.G:
                val = self.v(self.T.node[self.norder[0]], self.G.node[u])
                if self.keep(val, 0):
                    cs = 1 << self.color[u]
                    self.table[u][cs - 1] = (val,[u])
                    self.keys.add((u, cs))
                    if self.bestval < val - _tolerance:
                        self.bestval = val
                        self.bestkey = (u, cs)
                else:
                    #debug('Pruned ' + str(val) + ', ' + str((self.sbound, self.ibound[self.m - 1])))
                    self.pruned += 1
            #debug('TreeMatchLF.init: len(keys): ' + str(len(self.keys)))

        def step(i):
            keys = set()
            tmp = False # debug
            for (ku,s) in self.keys:
                q, l = self.table[ku][s-1]
                u = l[self.pre[i]] # contact point
                for v in self.G[u]: # all extensions
                    cv = 1 << self.color[v]
                    if cv & s:
                        continue
                    
                    cc = cv | s
                    val = (q
                           + self.v(self.T.node[self.norder[i]], self.G.node[v])
                           + self.w(self.T[self.norder[self.pre[i]]][self.norder[i]],
                                    self.G[u][v]))
                    
                    if not self.keep(val, i):
                        self.pruned += 1
                        continue
                    
                    # next contact point
                    newl = l + [v]
                    contact = newl[self.pre[i + 1]]
                    
                    if tmp:
                        debug('step ' + str(i) + ' contact '
                              + str(self.norder[self.pre[i]])
                              + ' connecting ' + str(self.norder[i])
                              + ' next contact '
                              + str(self.norder[self.pre[i + 1]]))
                        tmp = False

                    
                    qq,ll = self.table[contact][cc - 1]
                    if qq < val - _tolerance:
                        self.table[contact][cc - 1] = (val, newl)
                        keys.add((contact, cc))
                        if self.bestval < val - _tolerance:
                            self.bestval = val
                            self.bestkey = (contact, cc)
            self.keys = keys


        bests, bestf = 0, []
        csi = 2**self.m - 2
        for i in range(self.times):
            init()
            map(step, range(1,self.m)) # go through all the steps
            
            if self.bestval > bests + _tolerance:
                bests = self.bestval
                tmp = [self.table[u][csi] for u in self.G]
                bestf = [l for (s,l) in [self.table[u][csi] for u in self.G]
                         if s >= bests - _tolerance]
                continue
            if self.bestval >= bests - _tolerance:
                bestf += [l for (s,l) in [self.table[u][csi] for u in self.G]
                          if s >= bests - _tolerance]
            if uisim and bests == self.m * 2 - 1: # weights from [0,1] and 1-mad => max(bests) = 2m-1
                break
            if bests > self.sbound:
                self.sbound = bests
            debug('TreeMatchLF: ' + str(bests) + ', '
                  + str(self.sbound) + ', ' + str(self.pruned)
                  + ' uisim: ' + str((uisim, self.m * 2 - 1, bests == self.m*2 - 1)))

        info('TreeMatchLF: '+str((bests,len(bestf) if bestf != None else 0))
             + ' after ' + str(i+1) + ' trials.')
        bestf = [zip(self.norder, f) for f in bestf]
        debug('TreeMatchLF: bestf: ' + str(bestf))
        return (bestf, bests)

        
class GraphMatch:
    def __init__(self, G1, G2, k=5, eps=0.1,
                 quickbound = lambda T, G, f : (_tolerance, None), v=dsim, w=dsim):
        self.G1, self.G2, self.eps, self.v, self.w = G1, G2, eps, v, w
        self.m, self.n = len(G1), len(G2)
        self.quickbound, self.k = quickbound, k
        self.check = lambda T, G, a : (sum(v(T.node[i], G.node[a[i]]) for i in T if a[i] in G)
                                          +
                                          sum(w(d, G[a[i]][a[j]]) for i,j,d in T.edges(data=True) if
                                              a[i] in G and a[j] in G[a[i]]))
    def __call__(self, coverage=1, maxit = sys.maxint, callback=None, nocc=False,
                 star=True, uisim=False, single = False):
        B = defaultdict(lambda : defaultdict(lambda : 0))
        cover = dict(zip(self.G1, repeat(0)))
        minit = -1
        iteration = 0
        method = self.randomstar if star else self.randomsubtree
        treeit = cycle(imap(lambda x : lambda : method(x), self.G1))
        minit = self.m * coverage 
        while (iteration < minit):
            T = treeit.next()()
            if len(T) < 1:
                continue
            sbound, qf = self.quickbound(T, self.G2, self.check) 
            flist, s = (TreeMatch(T, self.G2, eps=self.eps, v=self.v,
                              w=self.w,sbound=sbound - _tolerance)(maxit = maxit, uisim=uisim, single=single)
                    if not nocc else ([qf.items()], sbound))
            if flist == None:
                continue
            flistlen = len(flist)
            debug('GraphMatch: len(flist) = ' + str(flistlen))
            for f in flist:
                for a,b in f:
                    B[a][b] += s/float(flistlen)
                    cover[a] += 1/float(flistlen)
            if callback != None:
                callback(T, flist, s, cover, B)
            iteration += 1
            info('GraphMatch: done with iteration ' + str(iteration) + ((' of ' + str(minit)) if minit > 0 else '')
                 + ' computed ' + str(flistlen) + ' assignments (' 
                 + str(reduce(min, cover.values())) + ').')
        return hungarian(B), strength(B)

    def randomsubtree(self, root):
        return randtree(self.G1, self.k, root)

    def randomstar(self, root):
        return randstar(self.G1, self.k, root)

def hungarian(hist):
    '''find optimal bipartite matching in bigraph represented by hist.'''
    rown = list(sorted(hist.keys()))
    coln = list(sorted(reduce(or_,
                              [set(hh.keys()) for hh in hist.values()],
                              set())))
    n = max(len(rown), len(coln))
    # create weight matrix
    mat = [[0]*n for i in range(n)]
    maxv = 0
    for i, row in enumerate(rown):
        for j, col in enumerate(coln):
            if hist[row].has_key(col):
                v = hist[row][col]
                mat[i][j] = hist[row][col]
                if maxv < v:
                    maxv = v
    #make costs out of weight
    mat = [[maxv - v for v in row] for row in mat]
    m = Munkres()
    hass = m.compute(mat)
    rdict = dict((rown[i],
                  coln[j] if j < len(coln) else None)
                 for (i,j) in hass if i < len(rown))
    return defaultdict(lambda : 0, rdict)

        

def npairs(s):
    '''computes sum_{i=0}^{len(s) - 2} sum_{j=i+1}^{len(s) -1} s[i]*s[j].

       if s is the list of sizes of equivalence classes induced by a function
       on a set, the above is the number of pairs in that set that the
       function discerns between.
    '''
    acc, prev = 0, 0
    for x in s:
        acc += x * prev
        prev = x + prev
    return acc

def flatness(h):
        '''a "sharpness" type measure for histograms.

        The fewer pairs it discerns between, the "sharper"'''
        vals = h.values()
        s = sum(vals)
        if s == 0:
            return 1 # uniform is really flat
        totalpairs = (s*(s - 1))/2.0
        return npairs(vals)/max(1.0, totalpairs)
        
def strength(B):
    return defaultdict(lambda : 0, ((key, 1-flatness(value)) for (key, value) in B.items()))
        




