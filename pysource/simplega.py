# -*-Python-*-
################################################################################
#
# File:         simplega.py
# RCS:          $Header: $
# Description:  A simple genetic algorithm that for finding injections
#               between two sets.
# Author:       Staal Vinterbo
# Created:      Mon May 18 11:17:09 2009
# Modified:     Wed Jul 20 15:33:07 2011 (Staal Vinterbo) staal@dink
# Language:     Python
# Package:      N/A
# Status:       Experimental
#
# (c) Copyright 2009, Staal Vinterbo, all rights reserved.
#
# simplega.py is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# simplega.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with simplega.py; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
################################################################################
__all__ = ['sga', 'gaalign']

from random import sample, choice, random, shuffle
from collections import defaultdict
import logging
import sys
from logging import debug, info, error, warning
from itertools import dropwhile, islice
from operator import add

def mix(a,b):
    '''randomly mix two injections, creating an injection

    such that the resulting injection c is such that
       c(x) in {a(x), b(x)}.

    Inputs are two lists of length m representing two injections
    a,b:{0,1,...,m} -> {0,1,...,n} for some n >= m (a sends i to a[i]).

    Algorithm sketch:
    The algorithm is an extension of cycle crossover with random selection
    of cycles (not altenating). This is extended to the injection case
    as explained below.

    If n = m we can view each of a and b as a perfect bipartite matching
    between Domain and Codomain.
    Let B be the bipartite graph constructed by the union of a and b.
    Each vertex in B has degree at most two, and B can be decomposed into
    disconnected components in which either all vertices have degree 1 or 2.
    All components where vertices have degree 1 are single edges, and
    represent the set of shared assigments. These are kept as they are.
    The remaining components are of vertice degree two, and are cycles.
    Of these cycles, we are only allowed to include every other edge for the
    result to be a perfect matching. Starting at a Domain vertex we can choose
    between the two directions we want to go, one direction includes
    edges from a, the other includes edges from b. We do this until we have
    traversed all cycles.

    If n > m, we can "ignore" cycle parts that have vertices outside
    of Domain and the union of Ranges. We simply stop processing an alternating
    path when we venture into the "ignored" part. The idea is to sample among
    all possible extensions of a and b into bijections and then processing
    as in the case where n = m.
    '''
    used = [False]*len(a)
    new = [None]*len(a)
    start = 0

    # compute back edges
    backa = defaultdict(lambda : None, [(j,i) for i,j in enumerate(a)])
    backb = defaultdict(lambda : None, [(j,i) for i,j in enumerate(b)])

    # step through each alternating path in random direction
    while start != None and not used[start]:
        # randomly decide path direction
        there, back = (a, backb) if random() > 0.5 else (b, backa)
        # follow alternating path until it ends or comes back to start
        x = start
        while x != None and not used[x]:
            new[x] = there[x]
            used[x] = True
            x = back[there[x]] 
        # find next unused start
        while start < len(a) and used[start]:
            start += 1
        if start == len(a): 
            start = None # no more possible starts, we are done
    return new
        
def wheel(weights, rnum = random()):
    '''roulette wheel sampling from a histogram'''
    choice = rnum * sum(weights)
    for i, w in enumerate(weights):
        choice -= w
        if choice < 0:
            return i

def sus(n, weights):
    '''stochastic universal sampling: roulette wheel with n markers'''
    if n > len(weights):
        warning('n > len(weights)')
    sw = float(sum(weights))
    choice = random()/n * sw
    ws = sorted([(w, i) for i,w in enumerate(weights)], reverse=True)
    j = 0
    for k in range(n):
        if j >= len(ws) or len(ws[j]) < 1:
            error(str(j) + ' >= ' + str(len(ws)) + ', ' + str((n,weights)))
        while choice - ws[j][0] >= 0:
            choice -= ws[j][0]
            j += 1
        yield ws[j][1]
        choice += sw/n

def split(what, lefts):
    '''partition list "what" into two lists (left, right)

    where left contains the elements (in order) that are indexed by
    "lefts", and right contains the remainder (in stable order).'''
    tmp = [[1, x] for x in what]
    for i in lefts:
        tmp[i][0] = 0
    l = [[],[]]
    for (i, x) in map(tuple,tmp):
        l[i].append(x)
    return tuple(l)


class sga:
    '''simple assignment search maximizing fitness genetic algorithm'''
    def __init__(self, fitness, m, n, psize = 128):
        self.psize = psize
        self.fitness = fitness
        self.m = m
        self.n = n
        self.pop = sorted(self.randpop(), reverse=True)

    def randpop(self):
        p = [sample(range(self.n), self.m) for i in xrange(self.psize)]
        return [(self.fitness(x), x) for x in p]

    def __xover(self, a, b):
        return [(self.fitness(x),x) for x in [mix(a,b), mix(a,b)]] 

    def __mutate(self, a):
        x = mix(a, sample(range(self.n), self.m))
        return (self.fitness(x), x)

    def crossover(self):
        '''selects 1/4 of population for parenthood and does mating'''
        idx = list(sus(self.psize/4, [w for (w,x) in self.pop]))
        shuffle(idx) # don't want to mate good with good and bad with bad
        offspring = []
        for k, l in [(idx[i], idx[i+1]) for i in range(0, len(idx), 2)]:
            offspring += self.__xover(self.pop[k][1], self.pop[l][1])
        return offspring

    def mutate(self, divisor=10):
        '''selects 10% for mutation and performs it'''
        idx = sus(self.psize/divisor, [w for (w,x) in self.pop])
        return map(self.__mutate, [self.pop[i][1] for i in idx])

    def generation(self, method = 1, mdiv=10):
        '''apply operations and recombine: generational, elitist'''
        ox = self.crossover()
        om = self.mutate(mdiv)
        cands = sorted(self.pop + ox + om, reverse=True)
        if method == 0:
            return cands[:self.psize]
        if method == 1:
            idx = sorted(sus(self.psize - 1, [w for (w,x) in cands[1:]]))
            newpop = [cands[0]]
            for i in idx:
                newpop.append(cands[i+1])
            return newpop # cands[:self.psize]
        if method == 3:
            return sorted(self.randpop(), reverse=True)
        newpop = [cands[0]]
        ws = [w for (w,x) in cands[1:]]
        for i in range(self.psize):
            newpop.append(cands[wheel(ws) + 1])
        return newpop
            

    def run(self, delay = 100, maxgen = 100000, method=1, mdiv=10):
        '''run the search.'''
        best = self.pop[0]
        gen = 0
        wait = delay
        verbose = logging.getLogger().getEffectiveLevel() <=  logging.DEBUG
        div = mdiv
        while wait and gen < maxgen:
            if wait == delay/5: # getting desperate, upping mutation rate
                div /= 2
            if wait == delay/10:# even more desperate, upping mutation rate
                div /= 2
            self.pop = self.generation(method, mdiv=div)
            if best[0] < self.pop[0][0]:
                best = self.pop[0]
                wait = delay + 1
                div = mdiv
            wait -= 1
            gen += 1
            if verbose:
                print >> sys.stderr, "\033[80D\033[K", gen, wait, best[0],
        if verbose:
            print >> sys.stderr
        debug('Best score: ' + str(best[0]))
        return best
            
        

if __name__ == "__main__":
    import sys
    import numpy as np

    wait = int(sys.argv[1])
    psize = int(sys.argv[2])

    data = map(lambda line : map(float, line.split(' ')), sys.stdin)

    m = len(data[0])
    mmat = np.array(data[0:m])
    nmat = np.array(data[m:])

    
    logging.basicConfig(level=max(logging.DEBUG,
                                  logging.WARNING - 10*2),
                        format='%(levelname)-5s %(message)s')
    
    info('m: ' + str(m))
    info('n: ' + str(len(nmat)))

    f = lambda a : m*m - abs(mmat - nmat[np.ix_(a,a)]).sum()
    ga = sga(f, m, len(nmat), psize=64)
    (score, sol) = ga.run(delay=200, method=3)
    print (score,sol)


#########
    
def gaalign(G1, G2, check):
    '''use a genetic algorithm to find a vertex assignment.'''
    #if type(G1) == type(nx.DiGraph()):
    #    warning('The GA algorithm is not designed for directed graphs...')
    transinv1 = dict(enumerate(sorted(G1)))
    trans1 = dict((b,a) for a,b in transinv1.items())
    GG = type(G1)()
    GG.add_edges_from((trans1[a], trans1[b], w) for a,b,w in G1.edges(data=True))
    GG.add_nodes_from((n, G1.node[transinv1[n]]) for n in GG)
    transinv2 = dict(enumerate(sorted(G2)))
    trans2 = dict((b,a) for a,b in transinv2.items())
    GG2 = type(G2)()
    GG2.add_edges_from((trans2[a], trans2[b], w) for a,b,w in G2.edges(data=True))
    GG2.add_nodes_from((n, G2.node[transinv2[n]]) for n in GG2)    
    debug('creating ga...')
    fitness = lambda a : check(GG, GG2, dict(enumerate(a))) 
    ga = sga(fitness, len(GG), len(GG2), 64)
    info('Running GA...')
    assi = ga.run(delay=300, method=1, mdiv=10)[1]
    debug('translating result....' + str(assi))
    ass = dict((transinv1[i], transinv2[j]) for i,j in enumerate(assi))
    return (check(G1, G2, ass), ass)
