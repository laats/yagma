# yagma
Yet another graph matching approach

Author
  ~ Staal A. Vinterbo

Copyright
  ~ 2011--2015 Staal A. Vinterbo

Version
  ~ generated for 1.24

Availability
  ~ [GPL](http://www.gnu.org/copyleft/gpl.html)

Homepage
  ~ [<http://laats.github.io/sw/yagma/>](http://laats.github.io/sw/yagma/)

Download
  ~ [<http://laats.github.io/sw/yagma/dist>](http://laats.github.io/sw/yagma/dist)

Synopsis
--------

Python:

    from yagma import GraphMatch, TreeMatch
    m = GraphMatch(G1, G2, k=5, eps=0.1, 
                   quickbound = lambda T, G, f : (_tolerance, None),  
                   v=dsim, w=dsim)
    assignment, strength = m(coverage=1, maxit = sys.maxint,
                             callback=None, nocc=False, star=True, 
                             uisim=False, single = False)
    tm = TreeMatch(T, G, eps=0.1, sbound=0, v=dsim, w=dsim)
    assignments, score = tm(sbound=None,  maxit=sys.maxint,
                      uisim = True, single = True)

Standalone programs:

    $ python yagm.py -h
    $ yagm -help

**Yagma** is Yet Another Graph Matching Approach. It is both an
algorithm/approach, a [Python](http://www.python.org/) module with an
accompanying command line script, and an [OCaml](http://caml.inria.fr/)
implementation of the algorithm.

