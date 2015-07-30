#!/usr/bin/env python
# -*- coding: utf-8 -*-
#al_circos.py
#03/2014

__author__ = 'Igor Rodrigues da Costa'
__contact__ = 'igor.bioinfo@gmail.com'

import os
import argparse
import pickle

def argument_parser(hlp=False):

    parser = argparse.ArgumentParser(description = 'Prints circos list', add_help = False,
                                     argument_default = None, fromfile_prefix_chars = '@',
                                     formatter_class = argparse.RawTextHelpFormatter)
    default_out = os.getcwd() + '/'
    parser.add_argument("-h", "--help", action = "help", help = "Show this help message and exit.")
    parser.add_argument('-o', '--outpath', nargs = '?', type = str, default = default_out,\
                        dest = 'outpath', help = 'Path where the circos AL files.')
    parser.add_argument('-i', '--infile', nargs = '?', type = str, default = default_out,\
                        dest = 'infile', help = 'File where the best_hit dictionary was saved.')
    args = parser.parse_args().__dict__
    return args

def main(args):
    with open(args['infile'], 'r') as l:
        locus, dict_hits = pickle.load(l)
    als = [] #human als
    for k in dict_hits:
	if 'Homo' in k:
            homo_key = k
            with open(args['outpath']+k.split('/')[-1].replace('BLASTDB', '').split('_')[0]+'.list', 'w') as out:
                for al in dict_hits[k]:
                    als.append(al)
                    chromo = (dict_hits[k][al][0])
                    start = str(dict_hits[k][al][1])
                    end = str(dict_hits[k][al][2])
                    out.write('\t'.join(('hs' + chromo, start, end)) + '\n')
    for k in dict_hits:
        if not 'Homo' in k:
            with open(args['outpath']+k.split('/')[-1].replace('BLASTDB', '').split('_')[0]+'.list', 'w') as out:
                for al in dict_hits[k]:
                    if al in als:
                        chromo = (dict_hits[homo_key][al][0])
                        start = str(dict_hits[homo_key][al][1])
                        end = str(dict_hits[homo_key][al][2])
                        out.write('\t'.join(('hs' + chromo, start, end)) + '\n')
                    else:
                        print al, 'was not found in ', k #should never happen!

if __name__ == '__main__':
    args = argument_parser()
    main(args)
