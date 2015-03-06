#!/usr/bin/env python
###############################################################################
#
# __template__.py - Description!
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
import glob

from multiprocessing import Pool
from subprocess import Popen, PIPE

#from Bio import SeqIO
#from Bio.Seq import Seq

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

# put classes here 

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
    # objects
    pairsLookup     = {}
    workingGIDs    = []
    workingPIDs     = []
    
    # read in pairs file
    with open(args.pairs, 'r') as fh:
        header = fh.readline()
        for l in fh:
            commas = l.rstrip().split(",")
            gid1 = commas[0]
            gid2 = commas[1]
            pid = commas[2]
            try:
                pairsLookup[gid1][gid2] = pid
            except KeyError:
                pairsLookup[gid1] = {gid2:pid}
            try:
                pairsLookup[gid2][gid1] = pid
            except KeyError:
                pairsLookup[gid2] = {gid1:pid}
    # read in list of gids
    with open(args.gid_list, 'r') as fh:
        for l in fh:
            gid = l.rstrip()
            workingGIDs.append(gid)
            
    for i in range(len(workingGIDs)):
        for j in range(i+1, len(workingGIDs)):
            pid = pairsLookup[workingGIDs[i]][workingGIDs[j]]
            print ">"+pid+"_"
            workingPIDs.append(pid)
    
    
     
            
            
    """
# parse fasta file using biopython
for accession,sequence in SeqIO.to_dict(SeqIO.parse(c_file,"fasta")).items():
if accession in genomes_dict:
pass
else:
#print accession
genomes_dict[accession] = [len(sequence),img_id, sequence.seq
"""  
    


###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-p','--pairs', help="...")
    parser.add_argument('-g','--gid_list', help="...")
    #parser.add_argument('input_file2', help="gut_img_ids")
    #parser.add_argument('input_file3', help="oral_img_ids")
    #parser.add_argument('input_file4', help="ids_present_gut_and_oral.csv")
    #parser.add_argument('output_file', help="output file")
    #parser.add_argument('positional_arg3', nargs='+', help="Multiple values")
    #parser.add_argument('-X', '--optional_X', action="store_true", default=False, help="flag")
    #parser.add_argument('-X', '--optional_X', action="store_true", type=int,default=False, help="flag")

    # parse the arguments
    args = parser.parse_args()

    # do what we came here to do
    doWork(args)

###############################################################################
###############################################################################
###############################################################################
###############################################################################
