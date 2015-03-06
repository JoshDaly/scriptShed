#!/usr/bin/env python
###############################################################################
#
# __script_name__.py - description!
#
###############################################################################
# #
# This program is free software: you can redistribute it and/or modify #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation, either version 3 of the License, or #
# (at your option) any later version. #
# #
# This program is distributed in the hope that it will be useful, #
# but WITHOUT ANY WARRANTY; without even the implied warranty of #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the #
# GNU General Public License for more details. #
# #
# You should have received a copy of the GNU General Public License #
# along with this program. If not, see <http://www.gnu.org/licenses/>. #
# #
###############################################################################

__author__ = "Josh Daly"
__copyright__ = "Copyright 2014"
__credits__ = ["Josh Daly"]
__license__ = "GPL3"
__version__ = "0.0.1"
__maintainer__ = "Josh Daly"
__email__ = ""
__status__ = "Development"

###############################################################################

import argparse
import sys

from multiprocessing import Pool
from subprocess import Popen, PIPE

#import os
#import errno

#import numpy as np
#np.seterr(all='raise')

#import matplotlib as mpl
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import axes3d, Axes3D
#from pylab import plot,subplot,axis,stem,show,figure


###############################################################################
###############################################################################
###############################################################################
###############################################################################

  # classes here

###############################################################################
###############################################################################
###############################################################################
###############################################################################

def runCommand(cmd):
    """Run a command and take care of stdout

expects 'cmd' to be a string like "foo -b ar"

returns (stdout, stderr)
"""
    p = Popen(cmd.split(' '), stdout=PIPE)
    return p.communicate()

def doWork( args ):
    """ Main wrapper"""
    #Colour array
    currentphylum="No"
    colours=("black","darkgreen","red","orange","yellow","mediumaquamarine","darkslategrey","grey")
    colourindex=-1
    colourmax=len(colours)
    phyla_parsing = False
    # read in OTU table
    with open(args.otu_table,"r") as fh:
        header = fh.readline().rstrip()
        if args.phyla_file == "T" or args.phyla_file == "True" or args.phyla_file == "TRUE":
            phyla_parsing = True
        else:
            print header
        lowest_rank = "NA"
        passes_cut_off = False
        for l in fh:
            tabs = l.split("\t")
            taxon_string = tabs[0]
            rel_abund = tabs[1:]
            ranks = taxon_string.split(";")
            phylum = ranks[1]
            for i in tabs[1:]:
                i = float(i)
                #print i
                if i >= args.cut_off:
                    passes_cut_off = True
                    #print "yep"
                    continue
                else: 
                    pass
            #print "\t".join([tabs[i] for i in range(1,len(tabs[1:]))])
            
            if phyla_parsing:
                if passes_cut_off:
                    if currentphylum != phylum:
                        colourindex+=1
                        if colourindex >= colourmax:
                                colourindex=0
                        
                        currentphylum=phylum
        
                    print colours[colourindex]
                    #print taxon_string+"\t"+("\t".join([rel_abund[i].rstrip() for i in range(len(rel_abund))]))
        
            else:    
                if passes_cut_off:
                    if ranks[5] != "Other":
                        print ranks[5]+"\t"+("\t".join([rel_abund[i].rstrip() for i in range(len(rel_abund))]))
                    else:
                        if ranks[4] != "Other":
                            print ranks[4]+"\t"+("\t".join([rel_abund[i].rstrip() for i in range(len(rel_abund))]))
                        else:
                            if ranks[3] != "Other":
                                print ranks[3]+"\t"+("\t".join([rel_abund[i].rstrip() for i in range(len(rel_abund))]))
                            else:
                                if ranks[2] != "Other":
                                    print ranks[2]+"\t"+("\t".join([rel_abund[i].rstrip() for i in range(len(rel_abund))]))
                                else: 
                                    if ranks[1] != "Other":
                                        print ranks[1]+"\t"+("\t".join([rel_abund[i].rstrip() for i in range(len(rel_abund))]))
                                    else:
                                        if ranks[0] != "Other":
                                            print ranks[0]+"\t"+("\t".join([rel_abund[i].rstrip() for i in range(len(rel_abund))]))
            passes_cut_off = False
            
            
            
                
            
            
            
            
            
    
    

    """
# run somethign external in threads
pool = Pool(6)
cmds = ['ls -l', 'ls -alh', 'ps -ef']
print pool.map(runCommand, cmds)
"""

    """
# parse a file
try:
with open(filename, "r") as fh:
for line in fh:
print line
except:
print "Error opening file:", filename, exc_info()[0]
raise
"""

    """
fig = plt.figure()

#-----
# make a 3d plot
ax = fig.add_subplot(111, projection='3d')
ax.scatter(points[:,0],
points[:,1],
points[:,2],
#edgecolors='none',
#c=colors,
#s=2,
#marker='.'
)

#-----
# make a 2d plot
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(points[:,0],
points[:,1],
'*g')

#-----
# show figure
plt.show()
# or save figure
plt.savefig(filename,dpi=300,format='png')

#-----
# clean up!
plt.close(fig)
del fig
"""

    return 0

###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--otu_table', help="QIIME style OTU table")
    parser.add_argument('-c','--cut_off',type=float,help="Relative abundance cut off %")
    parser.add_argument('-p','--phyla_file',default=False,help="Phyla file output i.e. T, True or TRUE: default False")
    #parser.add_argument('-o','--otu_file',default=False,help="OTU table output i.e. T or F: default True")
    #parser.add_argument('output_file', help="output file")
    #parser.add_argument('positional_arg3', nargs='+', help="Multiple values")
    #parser.add_argument('-X', '--optional_X', action="store_true", default=False, help="flag")

    # parse the arguments
    args = parser.parse_args()

    # do what we came here to do
    doWork(args)

###############################################################################
###############################################################################
###############################################################################
###############################################################################
