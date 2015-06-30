#!/usr/bin/env python
###############################################################################
#
# __annotateHits__.py - Annotate LGT hits!
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
__copyright__ = "Copyright 2015"
__credits__ = ["Josh Daly"]
__license__ = "GPL3"
__version__ = "0.0.1"
__maintainer__ = "Josh Daly"
__email__ = "joshua.daly@uqconnect.edu.au"
__status__ = "Development"

###############################################################################

# system imports
import argparse
import regex
import sys
import timeit
import pyfasta
from Bio import SeqIO
import difflib

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class Fragment:

    def __init__(self):
        self.length    = 0

        self.f200 = {}
        self.f210 = {}
        self.f220 = {}
        self.f230 = {}
        self.f240 = {}
        self.f250 = {}
        self.iTAG926f = ['AAACTCAAAGGAATTGACGG',
                         'AAACTTAAAGGAATTGACGG',
                         'AAACTCAAATGAATTGACGG',
                         'AAACTTAAATGAATTGACGG',
                         'AAACTCAAAGGAATTGGCGG',
                         'AAACTTAAAGGAATTGGCGG',
                         'AAACTCAAATGAATTGGCGG',
                         'AAACTTAAATGAATTGGCGG']
        self.iTAG803f = set(['TTAGAGACCCGCGTAGTC',
                            'TTAGATACCCTGGTAGTC',
                            'TTAGATACCCCTGTAGTC',
                            'TTAGATACCCGCGTAGTC',
                            'TTAGAGACCCGTGTAGTC',
                            'TTAGAGACCCTTGTAGTC',
                            'TTAGATACCCCCGTAGTC',
                            'TTAGATACCCGTGTAGTC',
                            'TTAGAGACCCCGGTAGTC',
                            'TTAGAGACCCTGGTAGTC',
                            'TTAGAGACCCTCGTAGTC',
                            'TTAGAGACCCTAGTAGTC',
                            'TTAGATACCCTTGTAGTC',
                            'TTAGATACCCTCGTAGTC',
                            'TTAGAGACCCGGGTAGTC',
                            'TTAGAGACCCCAGTAGTC',
                            'TTAGAGACCCCCGTAGTC',
                            'TTAGATACCCCGGTAGTC',
                            'TTAGATACCCCAGTAGTC',
                            'TTAGATACCCTAGTAGTC',
                            'TTAGATACCCGAGTAGTC',
                            'TTAGATACCCGGGTAGTC',
                            'TTAGAGACCCCTGTAGTC',
                            'TTAGAGACCCGAGTAGTC'])
        self.iTAG1114f = set(['CAACGAGCGCAACCC',
                              'TAACGAGCGCAACCC'])
        self.iTAG27f = set(['AGAGTTTGATCCTGGCTCAG',
                            'AGAGTTTGATTCTGGCTCAG',
                            'AGAGTTTGATCATGGCTCAG',
                            'AGAGTTTGATTATGGCTCAG'])

    def idhit(self,
              id,
              seq,
              primerDict
              ):
        for primer in primerDict:
            #match=regex.findall("(%s){e<=5}" % (primer), str(seq))
            ## checks if more than one match! 
            #if len(match)==1:
            #ind = str(seq).index(match[0])
            
            a={id: regex.findall("(%s){e<=5}" % (primer), str(seq))}
            a={key:item for key, item in a.iteritems() if item}
                
            for key,item in a.iteritems():
                if len(item)>1:    
                    best_match=[difflib.SequenceMatcher(None, x, primer).ratio() for x in item]
                    best=item[ best_match.index(max(best_match)) ] 
                    ind = str(seq).index(best) + len(best)
                else:
                    best=item[0]
                    ind = str(seq).index(best) + len(best)
            
            if len(str(seq)[ind+20:ind+270])==250:
                self.f250[id] = str(seq)[ind+20:ind+270]
                self.f200[id] = str(seq)[ind+20:ind+220]
                self.f210[id] = str(seq)[ind+20:ind+230]
                self.f220[id] = str(seq)[ind+20:ind+240]
                self.f230[id] = str(seq)[ind+20:ind+250]
                self.f240[id] = str(seq)[ind+20:ind+260]
                return True
            
            else:
                continue
        return False

    def idhitWrapper(self,
                     id,
                     seq,
                     primer
                     ):
        primerDict = {}
        
        if primer.lower() == '803f':
            primerDict = self.iTAG803f
        
        elif primer.lower() == '926f':
            primerDict = self.iTAG926f
        
        elif primer.lower() == '27f':
            primerDict = self.iTAG27f
        
        elif primer.lower() == '1114f':
            primerDict = self.iTAG1114f
            
        return self.idhit(id, seq, primerDict)

    def read(self,
             input_sequence_path,
             primer
             ):
        print "Reading in sequences..."
        entries=pyfasta.Fasta(input_sequence_path)
        sequence_count=len(entries)
        primer_miss=0
        primer_hits=0
        print "Finished reading in sequences, scanning for primers..."
        for entry in entries.keys():
            hit = self.idhitWrapper(entry,
                                    entries[entry],
                                    primer
                                    )
            if hit:
                primer_hits+=1
            elif not hit:
                primer_miss+=1
            print "%i hits; %i misses in %i\r" % (primer_hits, primer_miss, sequence_count),
        if primer_hits!=sequence_count:
            print "WARNING: %s reads hit the primer(s) out of a file of %s. Writing misses into file: primer_misses.fa" % (str(primer_hits), str(sequence_count))
            with open("primer_misses.fa", "w") as primer_misses_out:
                for entry in primer_miss:
                    primer_misses_out.write('>%s\n' % entry.id)
                    primer_misses_out.write('%s\n' % str(entry.seq))
        elif primer_hits==sequence_count:
            print "All reads in a file of %s hit the primer(s)" % (str(sequence_count))
        print "Writing split fragments to files: f200.fa, f210.fa, f220.fa, f230.fa, f240.fa, f250.fa"
        self.writeFrags()
        return

    def writeFrags(self):
        with open('f200.fa', 'w') as w200:
            for key, item in self.f200.iteritems():
                w200.write('>%s\n' % key)
                w200.write('%s\n' % item)
        with open('f210.fa', 'w') as w210:
            for key, item in self.f210.iteritems():
                w210.write('>%s\n' % key)
                w210.write('%s\n' % item)
        with open('f220.fa', 'w') as w220:
            for key, item in self.f220.iteritems():
                w220.write('>%s\n' % key)
                w220.write('%s\n' % item)
        with open('f230.fa', 'w') as w230:
            for key, item in self.f230.iteritems():
                w230.write('>%s\n' % key)
                w230.write('%s\n' % item)
        with open('f240.fa', 'w') as w240:
            for key, item in self.f240.iteritems():
                w240.write('>%s\n' % key)
                w240.write('%s\n' % item)
        with open('f250.fa', 'w') as w250:
            for key, item in self.f250.iteritems():
                w250.write('>%s\n' % key)
                w250.write('%s\n' % item)

###############################################################################
###############################################################################
###############################################################################
###############################################################################

def runCommand(cmd):
    """Run a command and take care of stdout

expects 'cmd' to be a string like "foo -b ar"

returns (stdout, stderr)
"""
    print cmd
    args = shlex.split(cmd)
    p = subprocess.Popen(args) # shell=bash is not recommended. Only use when '>' must be in cmd. 
    return p.communicate()
    #p = Popen(cmd.split(' '), stdout=PIPE)
    #return p.communicate()

def doWork( args ):
    """ Main wrapper"""
    F = Fragment()
    F.read(args.inputSeq,
           args.primer
           )
    

###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help="--", dest='subparser_name')
    
    #------------------------------
    
    fragment_parser = subparsers.add_parser('fragment',
                                            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                            description='Fragment input sequences into amplicon sizes of 200, 210, 220, 230, 240 and 250',
                                            help='Fragment input sequences into amplicon sizes of 200, 210, 220, 230, 240 and 250',
                                            epilog='Joel Boyd')
    
    fragment_parser.add_argument('-i', '--inputSeq', help='Input query files', required=True)
    fragment_parser.add_argument('-p', '--primer',help='Primer sequence to use', required=True)

    # parse the arguments
    args = parser.parse_args()

    # do what we came here to do
    doWork(args)

###############################################################################
###############################################################################
###############################################################################
###############################################################################
