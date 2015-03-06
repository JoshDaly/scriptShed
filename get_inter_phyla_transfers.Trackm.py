#!/usr/bin/env python
###############################################################################
#
# __get_inter_phyla_transfers.Trackm__.py - Description!
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

class HitDataParser(object):
    def __init__(self,l):
        self.readHitData(l)
        
    def readHitData(self,l):
        tabs=l.rstrip().split(",")
        pass
    
class DirtyDataParser(object):
    def __init__(self,l):
        self.readDirtyData(l)
        
    def readDirtyData(self,l):
        underscore = l.rstrip().split("_")
        self.pid = underscore[0] 
        self.sqid = underscore[1]

class DBDataParser(object):
    def __init__(self,l):
        self.readDBData(l)
        
    def readDBData(self,l):
        commas = l.rstrip().split(',')
        self.pid        = commas[0]
        self.hid        = commas[1]
        self.len_1      = commas[2]
        self.len_2      = commas[3]
        self.ident      = commas[4]
        self.strand_1   = commas[5]
        self.cid_1      = commas[6]
        self.cid_2      = commas[7]
        self.sqid_1     = commas[8]
        self.sqid_2     = commas[9]
        self.start_1    = commas[10]
        self.start_2    = commas[11]
        self.strand_2   = commas[12]
        self.ani_comp   = commas[13]
        self.gid_1      = commas[17]
        self.gid_2      = commas[18]

class ContigLengthParser(object):
    def __init__(self,l):
        self.readContigLengthFile(l)
        
    def readContigLengthFile(self,l):
        tabs = l.rstrip().split("\t")
        self.cid            = tabs[0]
        self.contig         = tabs[1]
        self.contigLength   = tabs[2]
        
class MetaDataParser(object):
    def __init__(self,l):
        self.readMetaData(l)
        
    def readMetaData(self,l):
        tabs=l.rstrip().split("\t")
        self.gid = tabs[0]
        self.status = tabs[3]
        self.sequencingCentre = tabs[6]
        self.phylum = tabs[7]
        self.genus  = tabs[11]
        self.bodysite = tabs[57]
        self.genomeSize = tabs[71]
        self.scaffoldCount = tabs[73]
        self.horizontalTransferred = tabs[119]
        self.sequencingMethod = tabs[132]
        

class HitData(object):
    def __init__(self):
        self.hitData        = {}
        self.cidToPID       = {}
        self.contigLength   = {}
        self.metaData       = {}
        self.dirtyData      = {}
        
    def addDBHitData(self,file):
        #DBDataParser = DBDataParser()
        with open(file,'r') as fh:
            for l in fh: 
                if l[0:3] != "pid": # header 
                    DBP = DBDataParser(l)
                    DBP.readDBData(l)
                    self.hitData[DBP.hid] = [DBP.pid,       #0
                                             DBP.ident,     #1
                                             DBP.ani_comp,  #2
                                             DBP.gid_1,     #3
                                             DBP.len_1,     #4
                                             DBP.strand_1,  #5
                                             DBP.cid_1,     #6
                                             DBP.sqid_1,    #7
                                             DBP.start_1,   #8
                                             DBP.gid_2,     #9
                                             DBP.len_2,     #10
                                             DBP.strand_2,  #11
                                             DBP.cid_2,     #12
                                             DBP.sqid_2,    #13
                                             DBP.start_2    #14
                                             ]
        
    def addMetaData(self,file):
        #MDP = MetaDataParser()
        with open(file,'r') as fh:
            header = fh.readline() # capture header 
            for l in fh:
                MDP = MetaDataParser(l)
                MDP.readMetaData(l)
                self.metaData[MDP.gid] = [MDP.bodysite,                 #0
                                          MDP.phylum,                   #1
                                          MDP.genus,                    #2
                                          MDP.status,                   #3
                                          MDP.sequencingMethod,         #4
                                          MDP.sequencingCentre,         #5
                                          MDP.horizontalTransferred,    #6
                                          MDP.genomeSize,               #7
                                          MDP.scaffoldCount             #8
                                          ]
        
    def addDirtyData(self,file):
        #DDP = DirtyDataParser()
        with open(file,'r') as fh:
            for l in fh:
                DDP = DirtyDataParser(l)
                DDP.readDirtyData(l)
                # lookup table for dirty sqids
                self.dirtyData[DDP.sqid] = 1
    
    def addContigLength(self,file):
        #CLP = ContigLengthParser()
        with open(file,'r') as fh:
            for l in fh:
                CLP = ContigLengthParser(l)
                CLP.readContigLengthFile(l)
                self.contigLength[CLP.cid] = [CLP.contig,       #0
                                              CLP.contigLength  #1
                                              ]
    
    def isDirty(self):
        for hid in self.hitData.keys():
            if self.hitData[hid][7] in self.dirtyData or self.hitData[hid][13] in self.dirtyData: 
                self.hitData[hid] += '1' 
            else: 
                self.hitData[hid] += '0'
               
    def printOutAllTheThings(self):
        for hid in self.hitData.keys():
            print "\t".join([hid,
                             self.hitData[hid][0],
                             self.hitData[hid][2],
                             self.hitData[hid][1],
                             self.hitData[hid][3],
                             self.metaData[self.hitData[hid][3]][0],
                             self.metaData[self.hitData[hid][3]][1],
                             self.metaData[self.hitData[hid][3]][2],
                             self.metaData[self.hitData[hid][3]][3],
                             self.metaData[self.hitData[hid][3]][4],
                             self.metaData[self.hitData[hid][3]][5],
                             self.metaData[self.hitData[hid][3]][6],
                             self.metaData[self.hitData[hid][3]][7],
                             self.metaData[self.hitData[hid][3]][8],
                             self.hitData[hid][4],
                             self.hitData[hid][5],
                             self.hitData[hid][6],
                             self.contigLength[self.hitData[hid][6]][0],
                             self.contigLength[self.hitData[hid][6]][1],
                             self.hitData[hid][7],
                             self.hitData[hid][8],
                             self.hitData[hid][9],
                             self.metaData[self.hitData[hid][9]][0],
                             self.metaData[self.hitData[hid][9]][1],
                             self.metaData[self.hitData[hid][9]][2],
                             self.metaData[self.hitData[hid][9]][3],
                             self.metaData[self.hitData[hid][9]][4],
                             self.metaData[self.hitData[hid][9]][5],
                             self.metaData[self.hitData[hid][9]][6],
                             self.metaData[self.hitData[hid][9]][7],
                             self.metaData[self.hitData[hid][9]][8],
                             self.hitData[hid][10],
                             self.hitData[hid][11],
                             self.hitData[hid][12],
                             self.contigLength[self.hitData[hid][12]][0],
                             self.contigLength[self.hitData[hid][12]][1],
                             self.hitData[hid][13],
                             self.hitData[hid][14],
                             self.hitData[hid][15]
                             ])
    
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

def printHeader():
    print "\t".join(["hid",
                     "pid",
                     "ani_comp",
                     "ident",
                     "gid_1",
                     "bodysite_1",
                     "phylum_1",
                     "genus_1",
                     "status_1",
                     "sequencingMethod_1",
                     "sequencingCentre_1",
                     "horizontalTransferred_1",
                     "genomeSize_1",
                     "scaffoldCount_1",
                     "len_1",
                     "strand_1",
                     "cid_1",
                     "contig_1",
                     "contigLength_1",
                     "sqid_1",
                     "start_1",
                     "gid_2",
                     "bodysite_2",
                     "phylum_2",
                     "genus_2",
                     "status_2",
                     "sequencingMethod_2",
                     "sequencingCentre_2",
                     "horizontalTransferred_2",
                     "genomeSize_2",
                     "scaffoldCount_2",
                     "len_2",
                     "strand_2",
                     "cid_2",
                     "contig_2",
                     "contigLength_2",
                     "sqid_2",
                     "start_2",
                     "dirty"
                     ])

def doWork( args ):
    """ Main wrapper"""  
    # Objects
    HD = HitData() 
        
    # Add data to HitData object
    printHeader()
    HD.addDBHitData(args.hitData)
    HD.addDirtyData(args.dirtyData)
    HD.isDirty()
    HD.addContigLength(args.contigLength)
    HD.addMetaData(args.metadata)
    HD.printOutAllTheThings()
    
        
    """    
    # read in hitData file
    with open(args.hitData, 'r') as fh:
        header = fh.readline() # capture header
        for l in fh:
            tabs        =  l.rstrip().split("\t")
            pid         = tabs[0]  
            gid1        = tabs[1]
            gid2        = tabs[10]
            phylum1     = tabs[4]
            phylum2     = tabs[13]
            hitsPost    = tabs[22]
            
            if phylum1 != phylum2:
                print pid
    """
            
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
    parser.add_argument('hitData', help="...")
    parser.add_argument('dirtyData', help="...")
    parser.add_argument('contigLength', help="...")
    parser.add_argument('metadata', help="...")
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
