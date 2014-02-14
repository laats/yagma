################################################################################
#
# File:         setup.py
# RCS:          $Header: $
# Description:  
# Author:       Staal Vinterbo
# Created:      Thu Jun  7 18:27:56 2007
# Modified:     Wed Oct  5 12:11:23 2011 (Staal Vinterbo) staal@dink
# Language:     Python
# Package:      sssm
# Status:       Experimental
#
# (c) Copyright 2007, Staal Vinterbo, all rights reserved.
#
# setup.py is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# setup.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with setup.py; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
################################################################################
"""yagma: Yet Another Graph Matching Approach."""

classifiers = """\
Development Status :: 3 - Alpha
Intended Audience :: Developers
Intended Audience :: Science/Research
License :: OSI Approved :: GNU General Public License (GPL)
Programming Language :: Python
Topic :: Scientific/Engineering :: Artificial Intelligence
Operating System :: Microsoft :: Windows
Operating System :: POSIX
Environment :: Console
"""

#from distutils.core import setup
from setuptools import setup


import pysource.docstring as docstring
version = docstring.version
scriptname = docstring.script
modulename = docstring.module
url = docstring.url

doclines = docstring.docstring.split("\n")


versionshort = version
pname = modulename

setup(name=pname,
      version=versionshort,
      author="Staal A. Vinterbo",
      author_email="sav@ucsd.edu",
      url = url,
      license = "http://www.gnu.org/copyleft/gpl.html",
      platforms = ["any"],
      description = doclines[0],
      classifiers = filter(None, classifiers.split("\n")),
      long_description = "\n".join(doclines[2:]),
      package_dir = {pname:'./pysource'},
      packages = [pname],
      install_requires = ['setuptools',
                          'munkres>=1.0.5', 'networkx>=1.5'],
      scripts=['./' + scriptname]
      )
