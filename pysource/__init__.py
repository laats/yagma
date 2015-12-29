# -*-Python-*-
################################################################################
#
# File:         __init__.py
# RCS:          $Header: $
# Description:  Graph matching by color coding
# Author:       Staal Vinterbo
# Created:      Tue Jul 12 14:24:45 2011
# Modified:     Tue Oct 13 12:14:05 2015 (Staal Vinterbo) staal@klump.gateway.pace.com
# Language:     Python
# Package:      N/A
# Status:       Experimental
#
# __init__.py is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# __init__.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with __init__.py; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
# (c) Copyright 2011, Staal Vinterbo, all rights reserved.
#
################################################################################
'''Graph matching by color coding'''
import graphmatch
import docstring
import external
import vf2

from simplega import gaalign
__all__ = graphmatch.__all__ + external.__all__ + vf2.__all__ + ['gaalign']
__doc__ = docstring.docstring
from external import *
from vf2 import *
from graphmatch import *
