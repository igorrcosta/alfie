#!/usr/bin/env python
# -*- coding: utf-8 -*-
#uce_finder.py
#05/2016

__author__ = 'Igor Rodrigues da Costa'
__contact__ = 'igor.bioinfo@gmail.com'

''' This program will find all UCEs read from a tab separated coord fileand save them in a FASTA file.'''


import os
import csv
import argparse
from math import floor
from shlex import split as ssplit
from Bio import SeqIO


def argument_parser(hlp=False):
    '''uce_finder.py'''

    default_out = os.getcwd() + '/uce_finder.out'
    default_log = os.getcwd() + '/uce_finder.log'
    description_example = '1\t1450000\n2\t3204000\ncontig123\t400000'

    parser = argparse.ArgumentParser(description = 'Finds UCEs', add_help = False,
                                     argument_default = None, fromfile_prefix_chars = '@',
                                     formatter_class = argparse.RawTextHelpFormatter)
    parser.add_argument("-h", "--help", action = "help", help = "Show this help message and exit.")
    parser.add_argument('-g', '--genome', nargs = '?', type = str, required = True,\
                        help = 'Fasta file with the genome.')
    parser.add_argument('-o', '--outfile', nargs = '?', type = str, default = default_out,\
                        dest = 'outfile', help = 'File where the UCEs will be saved.\n(default: %(default)s)')
    parser.add_argument('-l', '--log', nargs = '?', type = str, default = default_log,\
                        dest = 'log', help = 'Log file. (default: %(default)s)')
    parser.add_argument('--locus_length', nargs = '?', type = int, default = 1000,\
                        dest = 'length', help = 'Length of the locus sequences.\n(default: %(default)s)')
    parser.add_argument('--uce_distance', nargs = '?', type = int, default = 0,\
                        dest = 'uce_dist', help = 'Distance between locus and the beggining of the UCE.\n(defaut: %(default)s)')
    parser.add_argument('-u', '--uce', nargs = '?', type = str, default = None,\
                        dest = 'uce', help = 'File with the UCE coordinates.\n')
    parser.add_argument('-v', '--verbose', action = 'store_true', dest = 'verbose', help = 'Verbose switch.')


    if hlp:
        args = parser.print_help()
    else:
        args = parser.parse_args().__dict__
    return args

def main(args):
    '''Main function'''

    #args processing:
    if args['verbose']:
        a = open(args['log'], 'w')
        a.close()
        def vprint(*a):
            # Print only when verbose
            with open(args['log'], 'a') as log:
                log.write(' '.join(a))
    else:
        def vprint(*a):
            return None
    for k in args:
        vprint(k, ' ', str(args[k]), '\n')
    locus_length = args['length'] #size of the locus.
    uce_dist = args['uce_dist'] #distance between the locus and the UCE.
    try:
        o = open(args['outfile'], 'w')
        o.close()
    except:
        vprint('Can\'t open outfile: ' + args['outfile'] + '\n')
        raise
    uce_dict = get_uce(args['uce']) #uce_dict = {chr:[(uce_start, uce_end),...]}
    n = 1 #UCE counter
    for c in SeqIO.parse(args['genome'], 'fasta'): #for every chr
        chromo = c.description.split()[0]
        print chromo
        if chromo in uce_dict:
            for uce in uce_dict[chromo]:
                uce_start = int(uce[0]) - uce_dist
                uce_end = int(uce[1]) + uce_dist
                #54321uce12345 distance 0, lenght 10
                #54321_uce_12345 distance 1, lenght 10
                locus_start = uce_start - locus_length/2
                locus_end = uce_end + locus_length/2
                locus_seq1 = c.seq[locus_start:uce_start]#Get first half seq from the chromossome
                locus_seq2 = c.seq[uce_end:locus_end] #Get second half seq from the chromossome
                locus_seq = locus_seq1 + locus_seq2
                with open(args['outfile'], 'a') as out: #write the uce to outfile
                    out.write('>UCE_{}|{}:{}:{}\n{}\n'.format(n,
                               chromo, locus_start, locus_end, locus_seq))
                    n += 1


def get_uce(uce_file):
    uce_dict = {}
    with open(uce_file, 'r') as u_file:
        for l in u_file:
            if len(l.split()) == 4:
                #chromossome = l.split()[1].replace('chr', '')
                chromossome = l.split()[1]
                start = l.split()[2]
                end = l.split()[3]
                uce_dict[chromossome] = uce_dict.get(chromossome, []) + [(start, end)]
    return uce_dict

if __name__ == "__main__":
    args = argument_parser()
    for a in args:
        if isinstance(args[a], int) and args[a] < 0:
            raise ValueError('Argument "' + a + '" can\'t receive a negative value.')
    l = main(args)
