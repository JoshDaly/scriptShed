#!/usr/bin/env python
###############################################################################
#
# __parse_hit_data.trackm__.py - description!
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
import networkx as nx

#from multiprocessing import Pool
#from subprocess import Popen, PIPE

#import os
#import errno

#import numpy as np
#np.seterr(all='raise')

#import matplotlib as mpl
import matplotlib.pyplot as plt
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
        tabs = l.rstrip().split("\t")
        self.hid=                        tabs[0]
        self.pid=                        tabs[1]
        self.ani_comp=                   tabs[2]
        self.ident=                      tabs[3]
        self.gid_1=                      tabs[4]
        self.bodysite_1=                 tabs[5]
        self.phylum_1=                   tabs[6]
        self.genus_1=                    tabs[7]
        self.status_1=                   tabs[8]                       
        self.sequencingMethod_1=         tabs[9]
        self.sequencingCentre_1=         tabs[10]
        self.horizontalTransferred_1=    tabs[11]
        self.genomeSize_1=               tabs[12]
        self.scaffoldCount_1=            tabs[13]
        self.len_1=                      tabs[14]
        self.strand_1=                   tabs[15]
        self.cid_1=                      tabs[16]
        self.contig_1=                   tabs[17]
        self.contigLength_1=             tabs[18]
        self.sqid_1=                     tabs[19]
        self.start_1=                    tabs[20]
        self.gid_2=                      tabs[21]
        self.bodysite_2=                 tabs[22]
        self.phylum_2=                   tabs[23]
        self.genus_2=                    tabs[24]
        self.status_2=                   tabs[25]
        self.sequencingMethod_2=         tabs[26]
        self.sequencingCentre_2=         tabs[27]
        self.horizontalTransferred_2=    tabs[28]
        self.genomeSize_2=               tabs[29]
        self.scaffoldCount_2=            tabs[30]
        self.len_2=                      tabs[31]
        self.strand_2=                   tabs[32]
        self.cid_2=                      tabs[33]
        self.contig_2=                   tabs[34]
        self.contigLength_2=             tabs[35]
        self.sqid_2=                     tabs[36]
        self.start_2=                    tabs[37]
        self.dirty=                      tabs[38]

class PhiXDataParser(object):
    def __init__(self,l):
        self.readPhixData(l)
    
    def readPhixData(self,l):
        underscore = l.rstrip().split("_")
        self.pid  = underscore[0] 
        self.sqid = underscore[1]
        
class HitData(object):
    def __init__(self):
        self.interPhylaData =       {}
        self.scatterX =             [] 
        self.scatterY =             []
        self.scatterMarkerColour  =  []
        self.scatterContigvContigX= []
        self.scatterContigvContigY= []
        self.scatterDirtyTransfers= []
        self.PhiXHits             = {}
        self.colourByphix         = []
        self.phixSize             = []
        self.dirtySize            = []
        
        # Contig cutoff analysis
        self.phylum_interaction_map_cutoff        = {}
        self.phylum_interaction_map_cutoff_nodes  = {}
        self.keepTrackOfPhyla                     = {}
        self.gids                                 = {}
        
    def phixContaminatedHits(self,file):
        with open(file,'r') as fh:
            for l in fh:
                PDP = PhiXDataParser(l)
                self.PhiXHits[PDP.sqid] = 1
                
    def getPhylaCounts(self,file):
        with open(file,'r') as fh:
            for l in fh:
                HDP = HitDataParser(l)
                if l[0:3] != 'hid':
                    HDP.readHitData(l)
                    if HDP.gid_1 in self.gids: 
                        pass
                    else: 
                        self.gids[HDP.gid_1] = 1
                        try:
                            self.keepTrackOfPhyla[HDP.phylum_1] += 1 
                        except KeyError:
                            self.keepTrackOfPhyla[HDP.phylum_1] = 1
                    
                    if HDP.gid_2 in self.gids: 
                        pass
                    else: 
                        self.gids[HDP.gid_2] = 1
                        try:
                            self.keepTrackOfPhyla[HDP.phylum_2] += 1 
                        except KeyError:
                            self.keepTrackOfPhyla[HDP.phylum_2] = 1
        
        for phylum in self.keepTrackOfPhyla.keys():
            print "%s: %d" % (phylum,self.keepTrackOfPhyla[phylum])
    
    
    def getDataWithSizeCutoff(self,file,contig_size_cutoff):
        with open(file,'r') as fh:
            for l in fh:
                HDP = HitDataParser(l)
                if l[0:3] != 'hid':
                    HDP.readHitData(l)
                    if HDP.phylum_1 != HDP.phylum_2: #inter phyla 
                        if int(HDP.contigLength_1) >= contig_size_cutoff and int(HDP.contigLength_2) >= contig_size_cutoff:
                            if HDP.gid_1 in self.gids: 
                                pass
                            else: 
                                self.gids[HDP.gid_1] = 1
                                try:
                                    self.keepTrackOfPhyla[HDP.phylum_1] += 1 
                                except KeyError:
                                    self.keepTrackOfPhyla[HDP.phylum_1] = 1
                            
                            if HDP.gid_2 in self.gids: 
                                pass
                            else: 
                                self.gids[HDP.gid_2] = 1
                                try:
                                    self.keepTrackOfPhyla[HDP.phylum_2] += 1 
                                except KeyError:
                                    self.keepTrackOfPhyla[HDP.phylum_2] = 1
                            
                            self.phylum_interaction_map_cutoff_nodes[HDP.phylum_1] = 1
                            self.phylum_interaction_map_cutoff_nodes[HDP.phylum_2] = 1
                            if HDP.phylum_1 > HDP.phylum_2:
                                try:
                                    self.phylum_interaction_map_cutoff[HDP.phylum_1][HDP.phylum_2] += 1
                                except KeyError:
                                    try:
                                        self.phylum_interaction_map_cutoff[HDP.phylum_1][HDP.phylum_2] = 1 
                                    except KeyError:
                                        self.phylum_interaction_map_cutoff[HDP.phylum_1] = {HDP.phylum_2: 1}
                                        
                            else:
                                try:
                                    self.phylum_interaction_map_cutoff[HDP.phylum_2][HDP.phylum_1] += 1
                                except KeyError:
                                    try:
                                        self.phylum_interaction_map_cutoff[HDP.phylum_2][HDP.phylum_1] = 1 
                                    except KeyError:
                                        self.phylum_interaction_map_cutoff[HDP.phylum_2] = {HDP.phylum_1: 1}
        for phylum in self.keepTrackOfPhyla.keys():
            print "%s: %d" % (phylum,self.keepTrackOfPhyla[phylum])
    
    def PlotNetwork(self,showPlot,outfile):
        """Create network graph"""
        print "Building network plot:"
        print "Creating %d nodes" % (len(self.phylum_interaction_map_cutoff_nodes))
        G=nx.Graph()
        
        # add genes (nodes)
        G.add_nodes_from(self.phylum_interaction_map_cutoff_nodes.keys()) 
        
        # add edges
        for phylum1 in self.phylum_interaction_map_cutoff:
            for phylum2 in self.phylum_interaction_map_cutoff[phylum1]:
                    G.add_edge(phylum1,
                               phylum2,
                               label = str(self.phylum_interaction_map_cutoff[phylum1][phylum2]),
                               line_weight = self.phylum_interaction_map_cutoff[phylum1][phylum2])
                               #weight= self.phylum_interaction_map_cutoff[phylum1][phylum2])
                    
        
        # Build network plot
        #pos= nx.spring_layout(G,iterations=500)
        pos= nx.circular_layout(G)
        #fig = plt.figure(figsize=(21,10),dpi=300)
        
        #plt.subplot(1,1,1,axisbg='black',autoscale_on=False, aspect='equal', xlim=[-0.2,1.2], ylim=[-0.2,1.2])
        nx.draw_networkx_nodes(G,pos,linewidths=0, alpha=1, node_size=1000)
        #edge_width=dict([((u,v,),d['line_weight']) for u,v,d in G.edges(data=True)])
        nx.draw_networkx_edges(G, pos, edge_color = "#373737")
        edge_labels=dict([((u,v,),d['label']) for u,v,d in G.edges(data=True)])
        nx.draw_networkx_edge_labels(G,pos,edge_labels=edge_labels, font_size=7)
        nx.draw_networkx_labels(G,pos,font_size=4,font_family='sans-serif')
        plt.tick_params(axis='both',
                        which='both',
                        bottom='off',
                        top='off',
                        labelbottom='off',
                        labelleft='off')
        
        if showPlot:
            plt.show() 
        else: 
            plt.savefig(outfile,dpi=400,format='png')
        
    def getData(self,file,dataType, dataStatus):
        with open(file,'r') as fh:
            for l in fh:
                HDP = HitDataParser(l)
                if l[0:3] != 'hid':
                    HDP.readHitData(l)
                    if dataType.lower() == "inter":
                        if dataStatus.lower() == "all":
                            if HDP.phylum_1 != HDP.phylum_2:
                                print l.rstrip()
                        elif dataStatus.lower() == "clean":
                            if HDP.phylum_1 != HDP.phylum_2 and HDP.dirty == '0':
                                print l.rstrip()
                        elif dataStatus.lower() == "dirty":
                            if HDP.phylum_1 != HDP.phylum_2 and HDP.dirty == '1':
                                print l.rstrip()
                        else: 
                            print "Please indicate what status of data is to be used: all, clean or dirty!"
                            break
                    elif dataType.lower() == "intra":
                        if dataStatus.lower() == "all":
                            if HDP.phylum_1 == HDP.phylum_2:
                                print l.rstrip()
                        elif dataStatus.lower() == "clean":
                            if HDP.phylum_1 == HDP.phylum_2 and HDP.dirty == '0':
                                print l.rstrip()
                        elif dataStatus.lower() == "dirty":
                            if HDP.phylum_1 == HDP.phylum_2 and HDP.dirty == '1':
                                print l.rstrip()
                        else: 
                            print "Please indicate what status of data is to be used: all, clean or dirty!"
                            break
                    elif dataType.lower() == "both":
                        if dataStatus.lower() == "all":
                            print l.rstrip()
                        elif dataStatus.lower() == "clean":
                            if HDP.dirty == '0':
                                print l.rstrip()
                        elif dataStatus.lower() == "dirty":
                            if HDP.dirty == '1':
                                print l.rstrip()
                    else: 
                        print "Please indicate which data type you would prefer: intra or inter phyla!"
                        break

    def buildScatterPlot(self,file,dataType, dataStatus):
        # build two arrays, x = transfer length, y = contig length.
        with open(file,'r') as fh:
            for l in fh:
                HDP = HitDataParser(l)
                if l[0:3] != 'hid':
                    HDP.readHitData(l)
                    if dataType.lower() == "inter":
                        if dataStatus.lower() == "all":
                            if HDP.phylum_1 != HDP.phylum_2:
                                #print l.rstrip()
                                self.scatterX.append(int(HDP.len_1))
                                self.scatterY.append(int(HDP.contigLength_1))
                                self.scatterX.append(int(HDP.len_2))
                                self.scatterY.append(int(HDP.contigLength_2))
                                self.scatterContigvContigX.append(int(HDP.contigLength_1))
                                self.scatterContigvContigY.append(int(HDP.contigLength_2))
                                #########################
                                # Colour by genome status
                                #########################
                                if HDP.status_1.lower() == 'finished' and HDP.status_2.lower() == 'finished':
                                    self.scatterMarkerColour.append('red')
                                elif HDP.status_1.lower() == 'finished' and HDP.status_2.lower() != 'finished':
                                    self.scatterMarkerColour.append('yellow')
                                elif HDP.status_1.lower() != 'finished' and HDP.status_2.lower() == 'finished':    
                                    self.scatterMarkerColour.append('yellow')
                                else: 
                                    self.scatterMarkerColour.append('blue')
                                    
                                #########################
                                # Colour by dirty status
                                #########################    
                                if int(HDP.dirty) == 1:
                                    self.scatterDirtyTransfers.append('blue')
                                    self.dirtySize.append(20)
                                else: 
                                    self.scatterDirtyTransfers.append('yellow')
                                    self.dirtySize.append(0)  
                                
                                #############################
                                # Colour by Phix contamination
                                #############################
                                if HDP.sqid_1 in self.PhiXHits or HDP.sqid_2 in self.PhiXHits:
                                    self.colourByphix.append('blue')
                                    self.phixSize.append(20)
                                else:
                                    self.colourByphix.append('yellow')
                                    self.phixSize.append(0)   
                                    
                                    
                                
                        elif dataStatus.lower() == "clean":
                            if HDP.phylum_1 != HDP.phylum_2 and HDP.dirty == '0':
                                #print l.rstrip()
                                self.scatterX.append(int(HDP.len_1))
                                self.scatterY.append(int(HDP.contigLength_1))
                                self.scatterX.append(int(HDP.len_2))
                                self.scatterY.append(int(HDP.contigLength_2))
                                self.scatterContigvContigX.append(int(HDP.contigLength_1))
                                self.scatterContigvContigY.append(int(HDP.contigLength_2))
                                #########################
                                # Colour by genome status
                                #########################
                                if HDP.status_1.lower() == 'finished' and HDP.status_2.lower() == 'finished':
                                    self.scatterMarkerColour.append('red')
                                    self.scatterMarkerColour.append('red')
                                elif HDP.status_1.lower() == 'finished' and HDP.status_2.lower() != 'finished':
                                    self.scatterMarkerColour.append('red')
                                    self.scatterMarkerColour.append('blue')
                                elif HDP.status_1.lower() != 'finished' and HDP.status_2.lower() == 'finished':    
                                    self.scatterMarkerColour.append('blue')
                                    self.scatterMarkerColour.append('red')
                                else: 
                                    self.scatterMarkerColour.append('blue')
                                    self.scatterMarkerColour.append('blue')
                        elif dataStatus.lower() == "dirty":
                            if HDP.phylum_1 != HDP.phylum_2 and HDP.dirty == '1':
                                #print l.rstrip()
                                self.scatterX.append(int(HDP.len_1))
                                self.scatterY.append(int(HDP.contigLength_1))
                                self.scatterX.append(int(HDP.len_2))
                                self.scatterY.append(int(HDP.contigLength_2))
                                self.scatterContigvContigX.append(int(HDP.contigLength_1))
                                self.scatterContigvContigY.append(int(HDP.contigLength_2))
                                
                                #########################
                                # Colour by genome status
                                #########################
                                if HDP.status_1.lower() == 'finished' and HDP.status_2.lower() == 'finished':
                                    self.scatterMarkerColour.append('red')
                                    self.scatterMarkerColour.append('red')
                                elif HDP.status_1.lower() == 'finished' and HDP.status_2.lower() != 'finished':
                                    self.scatterMarkerColour.append('red')
                                    self.scatterMarkerColour.append('blue')
                                elif HDP.status_1.lower() != 'finished' and HDP.status_2.lower() == 'finished':    
                                    self.scatterMarkerColour.append('blue')
                                    self.scatterMarkerColour.append('red')
                                else: 
                                    self.scatterMarkerColour.append('blue')
                                    self.scatterMarkerColour.append('blue')
                                    
                                #########################
                                # Colour by dirty status
                                #########################    
                                if int(HDP.dirty) == 1:
                                    self.scatterDirtyTransfers.append('blue')
                                else: 
                                    self.scatterDirtyTransfers.append('yellow')
                                
                                #############################
                                # Colour by Phix contamination
                                #############################
                                if HDP.sqid_1 in self.PhiXHits or HDP.sqid_2 in self.PhiXHits:
                                    self.colourByphix.append('blue')
                                    ###############################################################################################################333
                                else:
                                    self.colourByphix.append('yellow')
                                    
                                    
                        else: 
                            print "Please indicate what status of data is to be used: all, clean or dirty!"
                            break
                    elif dataType.lower() == "intra":
                        if dataStatus.lower() == "all":
                            if HDP.phylum_1 == HDP.phylum_2:
                                #print l.rstrip()
                                self.scatterX.append(int(HDP.len_1))
                                self.scatterY.append(int(HDP.contigLength_1))
                                self.scatterX.append(int(HDP.len_2))
                                self.scatterY.append(int(HDP.contigLength_2))
                                self.scatterContigvContigX.append(int(HDP.contigLength_1))
                                self.scatterContigvContigY.append(int(HDP.contigLength_2))
                                #########################
                                # Colour by genome status
                                #########################
                                if HDP.status_1.lower() == 'finished' and HDP.status_2.lower() == 'finished':
                                    self.scatterMarkerColour.append('red')
                                    self.scatterMarkerColour.append('red')
                                elif HDP.status_1.lower() == 'finished' and HDP.status_2.lower() != 'finished':
                                    self.scatterMarkerColour.append('red')
                                    self.scatterMarkerColour.append('blue')
                                elif HDP.status_1.lower() != 'finished' and HDP.status_2.lower() == 'finished':    
                                    self.scatterMarkerColour.append('blue')
                                    self.scatterMarkerColour.append('red')
                                else: 
                                    self.scatterMarkerColour.append('blue')
                                    self.scatterMarkerColour.append('blue')
                        elif dataStatus.lower() == "clean":
                            if HDP.phylum_1 == HDP.phylum_2 and HDP.dirty == '0':
                                #print l.rstrip()
                                self.scatterX.append(int(HDP.len_1))
                                self.scatterY.append(int(HDP.contigLength_1))
                                self.scatterX.append(int(HDP.len_2))
                                self.scatterY.append(int(HDP.contigLength_2))
                                self.scatterContigvContigX.append(int(HDP.contigLength_1))
                                self.scatterContigvContigY.append(int(HDP.contigLength_2))
                                #########################
                                # Colour by genome status
                                #########################
                                if HDP.status_1.lower() == 'finished' and HDP.status_2.lower() == 'finished':
                                    self.scatterMarkerColour.append('red')
                                    self.scatterMarkerColour.append('red')
                                elif HDP.status_1.lower() == 'finished' and HDP.status_2.lower() != 'finished':
                                    self.scatterMarkerColour.append('red')
                                    self.scatterMarkerColour.append('blue')
                                elif HDP.status_1.lower() != 'finished' and HDP.status_2.lower() == 'finished':    
                                    self.scatterMarkerColour.append('blue')
                                    self.scatterMarkerColour.append('red')
                                else: 
                                    self.scatterMarkerColour.append('blue')
                                    self.scatterMarkerColour.append('blue')
                        elif dataStatus.lower() == "dirty":
                            if HDP.phylum_1 == HDP.phylum_2 and HDP.dirty == '1':
                                #print l.rstrip()
                                self.scatterX.append(int(HDP.len_1))
                                self.scatterY.append(int(HDP.contigLength_1))
                                self.scatterX.append(int(HDP.len_2))
                                self.scatterY.append(int(HDP.contigLength_2))
                                self.scatterContigvContigX.append(int(HDP.contigLength_1))
                                self.scatterContigvContigY.append(int(HDP.contigLength_2))
                                #########################
                                # Colour by genome status
                                #########################
                                if HDP.status_1.lower() == 'finished' and HDP.status_2.lower() == 'finished':
                                    self.scatterMarkerColour.append('red')
                                    self.scatterMarkerColour.append('red')
                                elif HDP.status_1.lower() == 'finished' and HDP.status_2.lower() != 'finished':
                                    self.scatterMarkerColour.append('red')
                                    self.scatterMarkerColour.append('blue')
                                elif HDP.status_1.lower() != 'finished' and HDP.status_2.lower() == 'finished':    
                                    self.scatterMarkerColour.append('blue')
                                    self.scatterMarkerColour.append('red')
                                else: 
                                    self.scatterMarkerColour.append('blue')
                                    self.scatterMarkerColour.append('blue')
                        else: 
                            print "Please indicate what status of data is to be used: all, clean or dirty!"
                            break
                    elif dataType.lower() == "both":
                        if dataStatus.lower() == "all":
                            #print l.rstrip()
                            self.scatterX.append(int(HDP.len_1))
                            self.scatterY.append(int(HDP.contigLength_1))
                            self.scatterX.append(int(HDP.len_2))
                            self.scatterY.append(int(HDP.contigLength_2))
                            self.scatterContigvContigX.append(int(HDP.contigLength_1))
                            self.scatterContigvContigY.append(int(HDP.contigLength_2))
                            
                            #########################
                            # Colour by genome status
                            #########################  
                            if HDP.status_1.lower() == 'finished' and HDP.status_2.lower() == 'finished':
                                self.scatterMarkerColour.append('red')
                                self.scatterMarkerColour.append('red')
                            elif HDP.status_1.lower() == 'finished' and HDP.status_2.lower() != 'finished':
                                self.scatterMarkerColour.append('red')
                                self.scatterMarkerColour.append('blue')
                            elif HDP.status_1.lower() != 'finished' and HDP.status_2.lower() == 'finished':    
                                self.scatterMarkerColour.append('blue')
                                self.scatterMarkerColour.append('red')
                            else: 
                                self.scatterMarkerColour.append('blue')
                                self.scatterMarkerColour.append('blue')
                            
                            """                       
                            if HDP.status_1.lower() == 'finished' and HDP.status_2.lower() == 'finished':
                                self.scatterMarkerColour.append('red')
                            elif HDP.status_1.lower() == 'finished' and HDP.status_2.lower() != 'finished':
                                self.scatterMarkerColour.append('yellow')
                            elif HDP.status_1.lower() != 'finished' and HDP.status_2.lower() == 'finished':    
                                self.scatterMarkerColour.append('yellow')
                            else: 
                                self.scatterMarkerColour.append('blue')
                            """   
                            #########################
                            # Colour by dirty status
                            #########################    
                            if int(HDP.dirty) == 1:
                                self.scatterDirtyTransfers.append('blue')
                            else: 
                                self.scatterDirtyTransfers.append('yellow')
                            
                            #############################
                            # Colour by Phix contamination
                            #############################
                            if HDP.sqid_1 in self.PhiXHits or HDP.sqid_2 in self.PhiXHits:
                                self.colourByphix.append('blue')
                                self.phixSize.append(20)
                            else:
                                self.colourByphix.append('yellow')
                                self.phixSize.append(0)
                                
                            """
                            if HDP.status_1.lower() == 'finished' and HDP.status_2.lower() == 'finished':
                                self.scatterMarkerColour.append('red')
                                self.scatterMarkerColour.append('red')
                            elif HDP.status_1.lower() == 'finished' and HDP.status_2.lower() != 'finished':
                                self.scatterMarkerColour.append('red')
                                self.scatterMarkerColour.append('blue')
                            elif HDP.status_1.lower() != 'finished' and HDP.status_2.lower() == 'finished':    
                                self.scatterMarkerColour.append('blue')
                                self.scatterMarkerColour.append('red')
                            else: 
                                self.scatterMarkerColour.append('blue')
                                self.scatterMarkerColour.append('blue')
                            """
                        elif dataStatus.lower() == "clean":
                            if HDP.dirty == '0':
                                #print l.rstrip()
                                self.scatterX.append(int(HDP.len_1))
                                self.scatterY.append(int(HDP.contigLength_1))
                                self.scatterX.append(int(HDP.len_2))
                                self.scatterY.append(int(HDP.contigLength_2))
                                self.scatterContigvContigX.append(int(HDP.contigLength_1))
                                self.scatterContigvContigY.append(int(HDP.contigLength_2))
                                #########################
                                # Colour by genome status
                                #########################
                                if HDP.status_1.lower() == 'finished' and HDP.status_2.lower() == 'finished':
                                    self.scatterMarkerColour.append('red')
                                    self.scatterMarkerColour.append('red')
                                elif HDP.status_1.lower() == 'finished' and HDP.status_2.lower() != 'finished':
                                    self.scatterMarkerColour.append('red')
                                    self.scatterMarkerColour.append('blue')
                                elif HDP.status_1.lower() != 'finished' and HDP.status_2.lower() == 'finished':    
                                    self.scatterMarkerColour.append('blue')
                                    self.scatterMarkerColour.append('red')
                                else: 
                                    self.scatterMarkerColour.append('blue')
                                    self.scatterMarkerColour.append('blue')
                        elif dataStatus.lower() == "dirty":
                            if HDP.dirty == '1':
                                #print l.rstrip()
                                self.scatterX.append(int(HDP.len_1))
                                self.scatterY.append(int(HDP.contigLength_1))
                                self.scatterX.append(int(HDP.len_2))
                                self.scatterY.append(int(HDP.contigLength_2))
                                self.scatterContigvContigX.append(int(HDP.contigLength_1))
                                self.scatterContigvContigY.append(int(HDP.contigLength_2))
                                #########################
                                # Colour by genome status
                                #########################
                                if HDP.status_1.lower() == 'finished' and HDP.status_2.lower() == 'finished':
                                    self.scatterMarkerColour.append('red')
                                    self.scatterMarkerColour.append('red')
                                elif HDP.status_1.lower() == 'finished' and HDP.status_2.lower() != 'finished':
                                    self.scatterMarkerColour.append('red')
                                    self.scatterMarkerColour.append('blue')
                                elif HDP.status_1.lower() != 'finished' and HDP.status_2.lower() == 'finished':    
                                    self.scatterMarkerColour.append('blue')
                                    self.scatterMarkerColour.append('red')
                                else: 
                                    self.scatterMarkerColour.append('blue')
                                    self.scatterMarkerColour.append('blue')
                    else: 
                        print "Please indicate which data type you would prefer: intra or inter phyla!"
                        break
        
    def plotScatter(self,showPlot, outfile):
        plt.subplot(1,1,1,autoscale_on=True)#, xlim=[-0.2,10000], ylim=[-0.2,10000])#, xlim=[-0.2,1.2], ylim=[-0.2,1.2])
        plt.axis([-5000, max(self.scatterX)+5000,-100000,max(self.scatterY)+500000])
        plt.xlabel('Transferred Sequence Length (bp)')
        plt.ylabel('Contig Size (bp)')
        #plt.subplot(1,1,1,axisbg='black')
        plt.scatter(self.scatterX,
                    self.scatterY,
                    marker='o',
                    alpha = 0.1,
                    s = 10,
                    linewidths = 0,
                    c= self.scatterMarkerColour)
        if showPlot:
            plt.show() 
        else: 
            plt.savefig(outfile,dpi=400,format='png')     
            
    def plotScatterContigLengthvContigLength(self,showPlot, outfile):   
        plt.subplot(1,1,1,autoscale_on=True)#, xlim=[-0.2,10000], ylim=[-0.2,10000])#, xlim=[-0.2,1.2], ylim=[-0.2,1.2])
        plt.axis([-250000, max(self.scatterContigvContigX)+5000,-250000,max(self.scatterContigvContigY)+500000])
        plt.plot([749965,749965],[-250000,max(self.scatterContigvContigY)+500000],color='grey',linestyle='--')
        plt.plot([-250000,max(self.scatterContigvContigY)+500000],[749965,749965],color='grey',linestyle='--')
        plt.xlabel('Contig Size (bp)')
        plt.ylabel('Contig Size (bp)')
        #plt.subplot(1,1,1,axisbg='black')
         
        
        plt.scatter(self.scatterContigvContigX,
                    self.scatterContigvContigY,
                    marker='o',
                    alpha = 0.3,
                    s = 40,
                    linewidths = 0,
                    c= self.scatterMarkerColour)
        if showPlot:
            plt.show() 
        else: 
            plt.savefig(outfile,dpi=400,format='png')
        
    def plotScatterContigLengthvContigLengthSideBySide(self,showPlot, outfile):
        plt.subplot(1,2,1,autoscale_on=True)#, xlim=[-0.2,10000], ylim=[-0.2,10000])#, xlim=[-0.2,1.2], ylim=[-0.2,1.2])
        plt.axis([-250000, max(self.scatterContigvContigX)+5000,-250000,max(self.scatterContigvContigY)+500000])
        plt.plot([749965,749965],[-250000,max(self.scatterContigvContigY)+500000],color='grey',linestyle='--')
        plt.plot([-250000,max(self.scatterContigvContigY)+500000],[749965,749965],color='grey',linestyle='--')
        plt.xlabel('Contig Size (bp)')
        plt.ylabel('Contig Size (bp)')
        #plt.subplot(1,1,1,axisbg='black')
         
        
        plt.scatter(self.scatterContigvContigX,
                    self.scatterContigvContigY,
                    marker='o',
                    alpha = 0.3,
                    s = 40,
                    linewidths = 0,
                    c= self.scatterMarkerColour)
        
        plt.subplot(1,2,2,autoscale_on=True)#, xlim=[-0.2,10000], ylim=[-0.2,10000])#, xlim=[-0.2,1.2], ylim=[-0.2,1.2])
        plt.axis([-250000, max(self.scatterContigvContigX)+5000,-250000,max(self.scatterContigvContigY)+500000])
        plt.plot([749965,749965],[-250000,max(self.scatterContigvContigY)+500000],color='grey',linestyle='--')
        plt.plot([-250000,max(self.scatterContigvContigY)+500000],[749965,749965],color='grey',linestyle='--')
        plt.xlabel('Contig Size (bp)')
        plt.ylabel('Contig Size (bp)')
        #plt.subplot(1,1,1,axisbg='black')
         
        
        plt.scatter(self.scatterContigvContigX,
                    self.scatterContigvContigY,
                    marker='o',
                    alpha = 0.3,
                    s = 40,
                    linewidths = 0,
                    c= self.scatterDirtyTransfers)

        if showPlot:
            plt.show() 
        else: 
            plt.savefig(outfile,dpi=400,format='png')
    
    def plotScatterContigLengthvContigLengthSideBySideAndPhiXcontamination(self,showPlot, outfile):
        plt.subplot(1,3,1,autoscale_on=True)#, xlim=[-0.2,10000], ylim=[-0.2,10000])#, xlim=[-0.2,1.2], ylim=[-0.2,1.2])
        plt.axis([-250000, max(self.scatterContigvContigX)+5000,-250000,max(self.scatterContigvContigY)+500000])
        plt.plot([749965,749965],[-250000,max(self.scatterContigvContigY)+500000],color='grey',linestyle='--')
        plt.plot([-250000,max(self.scatterContigvContigY)+500000],[749965,749965],color='grey',linestyle='--')
        plt.xlabel('Contig Size (bp)')
        plt.ylabel('Contig Size (bp)')
        #plt.subplot(1,1,1,axisbg='black')
         
        
        plt.scatter(self.scatterContigvContigX,
                    self.scatterContigvContigY,
                    marker='o',
                    alpha = 0.3,
                    s = 40,
                    linewidths = 0,
                    c= self.scatterMarkerColour)
        
        plt.subplot(1,3,2,autoscale_on=True)#, xlim=[-0.2,10000], ylim=[-0.2,10000])#, xlim=[-0.2,1.2], ylim=[-0.2,1.2])
        plt.axis([-250000, max(self.scatterContigvContigX)+5000,-250000,max(self.scatterContigvContigY)+500000])
        plt.plot([749965,749965],[-250000,max(self.scatterContigvContigY)+500000],color='grey',linestyle='--')
        plt.plot([-250000,max(self.scatterContigvContigY)+500000],[749965,749965],color='grey',linestyle='--')
        plt.xlabel('Contig Size (bp)')
        plt.ylabel('Contig Size (bp)')
        #plt.subplot(1,1,1,axisbg='black')
         
        
        plt.scatter(self.scatterContigvContigX,
                    self.scatterContigvContigY,
                    marker='o',
                    alpha = 0.3,
                    s = self.dirtySize,
                    linewidths = 0,
                    c= self.scatterDirtyTransfers)
        
        plt.subplot(1,3,3,autoscale_on=True)#, xlim=[-0.2,10000], ylim=[-0.2,10000])#, xlim=[-0.2,1.2], ylim=[-0.2,1.2])
        plt.axis([-250000, max(self.scatterContigvContigX)+5000,-250000,max(self.scatterContigvContigY)+500000])
        plt.plot([749965,749965],[-250000,max(self.scatterContigvContigY)+500000],color='grey',linestyle='--')
        plt.plot([-250000,max(self.scatterContigvContigY)+500000],[749965,749965],color='grey',linestyle='--')
        plt.xlabel('Contig Size (bp)')
        plt.ylabel('Contig Size (bp)')
        #plt.subplot(1,1,1,axisbg='black')
        
        
        plt.scatter(self.scatterContigvContigX,
                    self.scatterContigvContigY,
                    marker='o',
                    alpha = 0.3,
                    s = self.phixSize,
                    linewidths = 0,
                    c= self.colourByphix)

        if showPlot:
            plt.show() 
        else: 
            plt.savefig(outfile,dpi=400,format='png')
    
###############################################################################
###############################################################################
###############################################################################
###############################################################################

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

def runCommand(cmd):
    """Run a command and take care of stdout

expects 'cmd' to be a string like "foo -b ar"

returns (stdout, stderr)
"""
    p = Popen(cmd.split(' '), stdout=PIPE)
    return p.communicate()

def doWork( args ):
    """ Main wrapper"""
    # call class
    HD = HitData()

    if args.printData:
        #HD.getPhylaCounts(args.hitData)
        HD.getDataWithSizeCutoff(args.hitData, args.contig_cutoff)
        # print header 
        #printHeader()
        # print out specified data
        #HD.getData(args.hitData,args.dataType,args.dataStatus)
    elif args.contigComparison:
        HD.phixContaminatedHits(args.phixFile)
        HD.buildScatterPlot(args.hitData,args.dataType,args.dataStatus)
        HD.plotScatterContigLengthvContigLength(args.showPlot, args.outfile)
        #HD.plotScatterContigLengthvContigLengthSideBySide(args.showPlot, args.outfile)
        #HD.plotScatterContigLengthvContigLengthSideBySideAndPhiXcontamination(args.showPlot, args.outfile)
    elif args.contig_cutoff:
        pass
        #HD.getDataWithSizeCutoff(args.hitData,args.contig_cutoff)
        #HD.PlotNetwork(args.showPlot, args.outfile)
    else: 
        HD.phixContaminatedHits(args.phixFile)
        HD.buildScatterPlot(args.hitData,args.dataType,args.dataStatus)
        HD.plotScatter(args.showPlot, args.outfile)

###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('hitData', help="File containing hit data in tab delimited format")
    parser.add_argument('-dt','--dataType', default= "both",help="Inter or intra phyla transfers, or both")
    parser.add_argument('-ds','--dataStatus', default = "all",help="Type of data to output e.g. all, clean or dirty")
    parser.add_argument('-p','--printData', default = False,help="Print out hit data in tab delimited format")
    parser.add_argument('-sp','--showPlot', default = False,help="Display plot on screen")
    parser.add_argument('-o','--outfile', default = False, help="outfile")
    parser.add_argument('-c','--contigComparison', default = False, help="Comparison of contig sizes between genomes")
    parser.add_argument('-phix','--phixFile', default = False, help="File containing pids_sqids on phix contaminated hits")
    parser.add_argument('-cutoff','--contig_cutoff', type=int, help="Set size cutoff for contig")
    #parser.add_argument('input_file2', help="gut_img_ids")
    #parser.add_argument('input_file3', help="oral_img_ids")
    #parser.add_argument('input_file4', help="ids_present_gut_and_oral.csv")
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
