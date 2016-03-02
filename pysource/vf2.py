# -*-Python-*-
################################################################################
#
# File:         vf2.py
# RCS:          $Header: $
# Description:  networkx vf2 helpers.
# Author:       Staal Vinterbo
# Created:      Tue Oct 13 12:09:20 2015
# Modified:     Tue Oct 13 12:13:02 2015 (Staal Vinterbo) staal@klump.gateway.pace.com
# Language:     Python
# Package:      N/A
# Status:       Experimental
#
# (c) Copyright 2015, Staal Vinterbo, all rights reserved.
#
################################################################################
'''vf2: networkx vf2 helpers.'''

__all__ = ['vf2match', 'gridvf2']

import networkx.algorithms.isomorphism as iso
from logging import debug, info
from itertools import islice, combinations
from collections import defaultdict

from .graphmatch import GraphMatch


def vf2match(G1, G2, nit=1, atol = 1e-9, rtol=10e-7):
    ef = iso.numerical_edge_match('weight', 1.0, rtol, atol)
    nf = iso.numerical_node_match('weight', 1.0, rtol, atol)
    gm = iso.GraphMatcher(G2, G1, nf, ef)
    it = list(islice(gm.subgraph_isomorphisms_iter(), 1))
    debug('it ' + str(it))
    defret = G2.nodes()[0]
    ass = defaultdict(lambda : defret)
    if len(it) > 0:
        for x in G1:
            if it[0].has_key(x):
                ass[x] = it[0][x]
        debug('vf2: ' + str(ass))
    else:
        info('VF2: could not find isomorphism!')
    val = 0 if len(ass) == 0 else GraphMatch(G1, G2).check(G1, G2, ass)
    return val, ass

def gridvf2(G1, G2, nstep = 10):
    info('gridvf2: nstep = ' + str(nstep))
    stepw = 1.0/(nstep + 1)
    parmv = [0.000001 + stepw*i for i in xrange(nstep)]
    grid = list(combinations(parmv, 2))
    val1, ass1 = 0, {}
    for (a,b) in grid:
        debug('gridvf2: trying ' + str((a,b)))
        val, ass = vf2match(G1, G2, a, b)
        if len(ass) == 0:
            next
        if val > val1:
            val1 = val
            ass1 = ass
            debug('gridvf2: new best:' + str(val))
    return val1, ass1
