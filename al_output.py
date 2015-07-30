#!/usr/bin/env python
# -*- coding: utf-8 -*-
#alfie.py
#05/2015

__author__ = 'Igor Rodrigues da Costa'
__contact__ = 'igor.bioinfo@gmail.com'

''' This program will wrap all ALs predictor scripts in a user-friendly interface.'''
import os
import argparse
from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Align import MultipleSeqAlignment
from Bio.AlignIO import _FormatToWriter
from Bio.AlignIO.NexusIO import NexusWriter
from Bio.Nexus import Nexus


def argument_parser(hlp=False):
    '''al_output.py'''

    default_out = os.getcwd()
    parser = argparse.ArgumentParser(description = 'Concatenate and convert AL files', add_help = False,
                                     argument_default = None, fromfile_prefix_chars = '@',
                                     formatter_class = argparse.RawTextHelpFormatter)
    parser.add_argument("-h", "--help", action = "help", help = "Show this help message and exit.")
    parser.add_argument('-i', '--inpath', nargs = '?', type = str, required = True,\
                        help = 'Path to folder with the AL files.')
    parser.add_argument('-o', '--outpath', nargs = '?', type = str, default = default_out,\
                        dest = 'outpath', help = 'Path where the output will be saved.\n(default: %(default)s)')
    parser.add_argument('--concatenate', action = 'store_true', dest = 'concatenate',\
                         help = 'Concatenate all files (output will always be a single file).')
    parser.add_argument('--nexus', action = 'store_true', dest = 'nexus', help = 'Convert AL alignments to Nexus format.')
    parser.add_argument('--phylip', action = 'store_true', dest = 'phylip', help = 'Convert AL alignments to Phylip format.')
    parser.add_argument('--remove_gaps', action = 'store_true', dest = 'remove_gaps', help = 'Remove gaps from the final alignment.')
    parser.add_argument('--remove_n', action = 'store_true', dest = 'remove_n', help = 'Remove "N"s from the final alignment.')
    if hlp:
        args = parser.parse_args(['-h'])
    else:
        args = parser.parse_args().__dict__
    return args

def main():
    args = argument_parser()
    if args['inpath'][-1] != '/':
        args['inpath'] += '/'
    if args['outpath'][-1] != '/':
        args['outpath'] += '/'
    args['aligned_files'] = [args['inpath'] + aln for aln in os.listdir(args['inpath']) if '.aln' in aln]
    print args['aligned_files']
    if args['nexus']:
        make_nexus(args['aligned_files'], args['outpath'], args['remove_gaps'], args['remove_n'])
    elif args['phylip']:
        make_philip(args['aligned_files'], args['outpath'], args['remove_gaps'], args['remove_n'])
    elif args['concatenate']:
        join_fasta(args['aligned_files'], args['outpath'], args['remove_gaps'], args['remove_n'])
    else:
        print 'No option selected.'
        argument_parser(hlp=True)

def join_fasta():
    pass

def make_philip():
    pass

def make_nexus(aligned_files, outpath, remove_gaps, remove_n):
    nexus_files = []
    for infile in aligned_files:
        outfile = outpath + infile.replace('.aln', '.nexus').split('/')[-1]
        fastatonexus(infile, outfile, keep_gaps = not remove_gaps, keep_n = not remove_n)
        nexus_files.append(outfile)

def fastatonexus(infile, outfile, format_in = 'fasta', format_out = 'nexus', keep_gaps = False, keep_n = False):
    align = ''
    gaps = 0
    with open(infile, 'r') as handle:
        i = AlignIO.read(handle, format_in)
        columns = len(i[0])
        for col in range(columns):
            if keep_n or 'N' not in str(i[:,col:col+1]):
                if '-' in str(i[:,col:col+1]):  
                    gaps += 1
                    if keep_gaps:
                        if align:
                            align += i[:,col:col+1]
                        else:
                            align = i[:,col:col+1]
                else:
                    if align:
                        align += i[:,col:col+1]
                    else:
                        align = i[:,col:col+1]
#    if not align:
#        vprint('No sequences found in', infile)
#    if keep_gaps:
#        size = align.get_alignment_length()
#    else:
#        size = align.get_alignment_length() + gaps
    align.sort()
    align._alphabet = IUPAC.unambiguous_dna
    with open(outfile, 'wb') as out:
        AlignIO.write(align, out, format_out)

if __name__ == '__main__':
    main()
