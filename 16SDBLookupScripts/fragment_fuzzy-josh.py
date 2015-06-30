#!/usr/bin/env python

import argparse
import regex
import sys
import timeit
import pyfasta
from Bio import SeqIO

def phelp():
    print """
-----------------------------------------------
               16S_db_lookup
-----------------------------------------------

  fragment    -   Build fragments from full
                  length sequences
                  -i <fasta> -t <gg_taxonomy>

  build       -   From a fragment file, build
                  a database, and return the
                  redundant sequences.
  store

Joel A. Boyd
===============================================

"""


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

    def idhit(self, id, seq):
        for primer in self.iTAG803f:
            match=regex.findall("(%s){e<=5}" % (primer), str(seq))
            if len(match)==1:
                ind = str(seq).index(match[0])
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
            else:
                continue
        return False


    def read(self, input_sequence_path):
        print "Reading in sequences..."
        entries=pyfasta.Fasta(input_sequence_path)
        sequence_count=len(entries)
        primer_miss=0
        primer_hits=0
        print "Finished reading in sequences, scanning for primers..."
        for entry in entries.keys():
            hit = self.idhit(entry, entries[entry])

            if hit:
                primer_hits+=1
            elif not hit:
                primer_miss+=1
            print "%i hits; %i misses in %i\r" % (primer_hits, primer_miss, sequence_count),
        if primer_hits!=sequence_count:
            print "WARNING: %s reads hit the primer(s) out of a file of %s. Writing misses into file: primer_misses.fa" % (str(primer_hits), str(sequence_count))
            #with open("primer_misses.fa", "w") as primer_misses_out:
            #    for entry in primer_miss:
            #        primer_misses_out.write('>%s\n' % entry.id)
            #        primer_misses_out.write('%s\n' % str(entry.seq))
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

    def main(self, args):
        print '\n--FRAGMENT--'
        self.read(args.i)
        print '\n'


class Lookup:

    def __init__(self):
        self.f = Fragment()

    def tax(self, gg_fmt_tax):
        taxonomy_hash = {}
        for line in [x.rstrip() for x in open(gg_fmt_tax).readlines()]:
            line_splt = line.split()
            taxonomy_hash[line_splt[0]] = line_splt[1]
        return taxonomy_hash

    def redundantTax(self, db_hash):

        def whatax(tax_list):
            levels={0:'Kingdom',
                    1:'Phylum',
                    2:'Class',
                    3:'Order',
                    4:'Family',
                    5:'Genus',
                    6:'Species'}
            seen=[]
            tax=[]
            seen_counter={}

            seen = list(set(tax_list))
            for item in seen:
                seen_counter[item] = tax_list.count(item)

            if len(seen)>1:
                r=[[y for y in x.split('; ') if y] for x in seen]
                ranges = [len(x) for x in r]
                if len(set(ranges)) == 1:
                    for i in range(0,sorted(ranges)[0]):
                        rank_classifications = [x[i] for x in r]
                        if len(set(rank_classifications))==1:
                            tax.append(r[0][i])
                        else:
                            divergence_dict = {}
                            divergence_counter = 1
                            for idx, rank in enumerate(rank_classifications):
                                divergence_dict["split_%s_%s" % (rank, str(seen_counter[seen[idx]]))] = '; '.join(r[idx][i:])
                            divergence_dict['2_split_rank'] = levels[i]
                            tax.append(divergence_dict)
                            return tax
                return tax

            elif len(seen)==1:
                return False
            else:
                raise Exception

        div=[]
        for key, item in db_hash.iteritems():
            if item['obs'] > 1:
                w = whatax(item['tax'])
                if w:
                    div_dict=w.pop(-1)
                    if div_dict == 'uncultured':
                        raise Exception

                    div_dict['1_split_consensus'] = '; '.join(w)
                    div_dict['0_seq']=key
                    div.append(div_dict)
                else:
                    pass

        a=[]
        for i in div:
            keys = [x for x in i.keys() if x.startswith("split") and int(x[-1]) > 5]
            if keys:
                a.append(i['2_split_rank'])

        counter=0
        for i in div:
            keys = [x for x in i.keys() if x.startswith("split")]
            for x in keys:
                if any([x for x in keys if int(x.split('_')[-1])<=5]):
                    continue
                elif i['2_split_rank']=="Kingdom":
                    counter+=1
        counter
        import IPython
        IPython.embed()

    def check(self, args):
        hashlist=[]
        mismatch_files=False
        def extd(seq, tax, lookup_tbl):

            lookup_tbl[seq] = {'tax': [tax],
                               'obs': 1,
                               'dis_tax': [tax]}
            return lookup_tbl

        def updt(seq, tax, lookup_tbl):
            if tax!=lookup_tbl[seq]['tax']:
                lookup_tbl[seq]['tax'].append(tax)
                lookup_tbl[seq]['obs']+=1
            else:
                print "Not sure what to do in this situation..."
            return lookup_tbl

        def srch(seq):
            if seq in set(self.lookup_tbl.keys()):
                return True
            elif seq not in set(self.lookup_tbl.keys()):
                return False
            else:
                raise Exception("Programming error in Lookup/check/srch")

        print 'Reading in taxonomy'
        taxonomy=self.tax(args.t) # Read in taxonomy file
        len_tax = len(taxonomy)
        print 'Finished reading taxonomy, %i entries found' % (len(taxonomy))
        start = timeit.default_timer()

        start = timeit.default_timer()
        print 'Reading in sequences'
        fragment_size=pyfasta.Fasta(args.i)
        len_seqs = len(fragment_size.keys())
        print "Done. Found %s sequences" % (len_seqs)
        stop = timeit.default_timer()



        if len_tax != len_seqs:
          print "WARNING: Taxonomy file and sequences don't match up. Are you sure you gave me the right taxonomy or sequence file?"
          print "..Tentatively continuing"
          mismatch_files=True
        self.lookup_tbl = {} # will later include import of pre-existing database
        howmany=0
        counter=1
        for item in fragment_size.keys():
            s = (True if str(fragment_size[item]) in self.lookup_tbl else False)

            try:
                tax = taxonomy[item]
            except KeyError:
                #print "Taxonomy for %s not found for some annoying reason. I can't work like this!" % str(item)
                continue
            if s:
                lookup_tbl=updt(str(fragment_size[item]), tax, self.lookup_tbl)
                howmany+=1
            elif not s:
                lookup_tbl=extd(str(fragment_size[item]), tax, self.lookup_tbl)
                howmany+=1
            counter+=1
            p_complete = round(float(counter)/float(len_seqs)*100, 1)

            if p_complete == 100.0:
                print '\n100% complete!\n'
                break
            else:
                #print 'Processed %s sequences (%s%% percent complete)\r' % (str(counter), str(p_complete)),
                continue
        print howmany
        exit()
        if p_complete != 100.0 and mismatch_files:
            print '\nWARNING: Not all sequences were processed. Probably because of the mismatched taxonomy and sequence files'
        elif p_complete != 100.0 and not mismatch_files:
            raise Exception('Programming Error')
        else:
            print "\nDone building database"
        return self.lookup_tbl

    def main(self, args):
        print "\n--BUILD--"
        result=self.check(args)
        print 'Determining redundant taxonomy'
        self.redundantTax(result)
        print '\n'

if __name__=="__main__":
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--version', action='version', version='16S_db_lookup.py v0.0.1')
    subparsers = parser.add_subparsers(help="--", dest='subparser_name')

    fragment_parser = subparsers.add_parser('fragment',
                                            description='Fragment input sequences into amplicon sizes of 200, 210, 220, 230, 240 and 250',
                                            epilog='Joel Boyd')
    fragment_parser.add_argument('-i', help='Input query files', required=True)
    fragment_parser.add_argument('-p', '--primer', help='Primer sequence to use', required=True)

    build_parser = subparsers.add_parser('build',
                                         description='Build a sequence database, and return redundant sequences',
                                         epilog='Joel Boyd')
    build_parser.add_argument('-i', help='Input query files', required=True)
    build_parser.add_argument('-t', help='Taxonomy. Required for now.', required=True)


    if(len(sys.argv) == 1 or sys.argv[1] == '-h' or sys.argv[1] == '--help'):
        phelp()
        exit(0)
    else:
        args = parser.parse_args()

    if args.subparser_name=='fragment':
        f = Fragment()
        f.main(args)
    elif args.subparser_name=='build':
        l = Lookup()
        l.main(args)
    else:
        raise Exception("Programming Error.")

