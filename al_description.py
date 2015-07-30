#!/usr/bin/env python
# -*- coding: utf-8 -*-
#al_description.py
#08/2014

__author__ = 'Igor Rodrigues da Costa'
__contact__ = 'igor.bioinfo@gmail.com'

''' This program will generate a table with the sizes from every fasta entry of a multifasta file.'''

import os
import argparse
from Bio import SeqIO


def argument_parser(hlp=False):

    default_out = os.getcwd() + '/description.txt'
    parser = argparse.ArgumentParser(description = 'Generates a table with the sizes from a multifasta file.', add_help = False,
                                     argument_default = None, fromfile_prefix_chars = '@',
                                     formatter_class = argparse.RawTextHelpFormatter)
    parser.add_argument("-h", "--help", action = "help", help = "Show this help message and exit.")
    parser.add_argument('-i', '--infile', nargs = '?', type = str, required = True,\
                        help = 'Multifasta file.')
    parser.add_argument('-o', '--outfile', nargs = '?', type = str, default = default_out,\
                        dest = 'outfile', help = 'File where the table will be saved.\n(default: %(default)s)')
    if hlp:
        args = parser.print_help()
    else:
        args = parser.parse_args().__dict__
    return args

def size_parser(infile, outfile):
    '''Prints a file with fasta name -> size.'''
    with open(outfile, 'w') as out:
        for seq in SeqIO.parse(infile, 'fasta'):
            name = seq.description.split()[0]
            size = len(seq.seq)
            out.write('\t'.join([name, str(size)])+'\n')

if __name__ == '__main__':
   args = argument_parser()
   size_parser(args['infile'], args['outfile'])
