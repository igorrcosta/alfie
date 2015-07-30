#!/usr/bin/env python
# -*- coding: utf-8 -*-
#alfie.py
#05/2015

__author__ = 'Igor Rodrigues da Costa'
__contact__ = 'igor.bioinfo@gmail.com'

''' This program will wrap all ALs predictor scripts in a user-friendly interface.'''

import os
import argparse
import al_finder
import al_formatdb
import al_blast
import al_align
from copy import deepcopy

def argument_parser(hlp=False):
    '''alfie.py'''

    default_out = os.getcwd()
    parser = argparse.ArgumentParser(description = 'Finds ALs', add_help = False,
                                     argument_default = None, fromfile_prefix_chars = '@',
                                     formatter_class = argparse.RawTextHelpFormatter)
    parser.add_argument("-h", "--help", action = "help", help = "Show this help message and exit.")
    parser.add_argument('-i', '--genomes', nargs = '*', type = str, required = True,\
                        help = 'Path to the fasta files with the genomes.')
    parser.add_argument('-o', '--outpath', nargs = '?', type = str, default = default_out,\
                        dest = 'outpath', help = 'Path where the ALs will be saved.\n(default: %(default)s)')
    parser.add_argument('-l', '--log', nargs = '?', type = str, default = 'al_finder.log',\
                        dest = 'log', help = 'Log file. (default: %(default)s)')
    parser.add_argument('-f', '--skip_formatdb', action = 'store_true', dest = 'skip_formatdb', help = 'Skip making BLAST databases, use databses from the last run.')
    parser.add_argument('--locus_length', nargs = '?', type = int, default = 1000,\
                        dest = 'length', help = 'Length of the ALs sequences.\n(default: %(default)s)')
    parser.add_argument('--max_n', nargs = '?', type = float, default = 0,\
                        dest = 'max_n', help = 'Maximum percentage of N\'s in the AL sequence.\n(default: %(default)s)')
    parser.add_argument('--inter_distance', nargs = '?', type = int, default = 200000,\
                        dest = 'idist', help = 'Minimum distance between ALs.\n(defaut: %(default)s)')
    parser.add_argument('--gene_distance', nargs = '?', type = int, default = 200000,\
                        dest = 'gdist', help = 'Minimum (or maximum, if negative) distance between ALs and genes.\n(defaut: %(default)s)')
    parser.add_argument('--gene_locus', action = 'store_true', dest = 'gene_locus',\
                        help = 'Find coding regions loci.\n(defaut: %(default)s)')
    parser.add_argument('--cds', action = 'store_true', default = False,\
                        dest = 'cds', help = 'Only considers the CDS features of GTF files. (default: %(default)s)')
    parser.add_argument('--end_distance', nargs = '?', type = int, default = 10000,\
                        dest = 'edist', help = 'Distance between ALs and the start and end of a chromosome.\n(default: %(default)s)')
    parser.add_argument('-g', '--gtf', nargs = '?', type = str, required = True,\
                        dest = 'est', help = 'GTF File with all genome features coordinates.')
    parser.add_argument('-v', '--verbose', action = 'store_true', dest = 'verbose', help = 'Verbose switch.')
    parser.add_argument('--duplication_cutoff', nargs = '?', type = int, default = 50,\
                        dest = 'dup_cut', help = 'ALs with 2 hits with identity higher than this will\
                        be considered duplicated.\n(default: %(default)s)')  
    parser.add_argument('--identity_cutoff', nargs = '?', type = int, default = 90,\
                        dest = 'id_cut', help = 'ALs with a identity higher than this will be considered homologous.\n(default: %(default)s)')
    parser.add_argument('--coverage_cutoff', nargs = '?', type = int, default = 90,\
                        dest = 'cov_cut', help = 'BLAST hits must have at least this much %%coverage to be considered hits.\n(default: %(default)s)')
    parser.add_argument('--chromossomes', nargs = '*', type = str,\
                        dest = 'excluded', help = 'Chromossomes to be excluded.')
    parser.add_argument('--min_align', nargs = '?', type = int, default = 500,\
                        dest = 'align_size', help = 'Minimum final alignment lenght.')
    parser.add_argument('--minsize', nargs = '?', type = int, default = 500,\
                        dest = 'min_size', help = 'Minimum sequence size.')
    parser.add_argument('--total_al', nargs = '?', type = int,\
                        dest = 'pick', help = 'Pick only N ALs.')
    parser.add_argument('--remove_gaps', action = 'store_true', dest = 'nogaps', help = 'Remove gaps from the final alignment.')
    if hlp:
        args = parser.print_help()
    else:
        args = parser.parse_args().__dict__
    return args

def main():
    #todo: improve log and verbose.
    args = argument_parser()
    for a in args:
        if isinstance(args[a], int) and args[a] < 0 and a not in ('idist', 'gdist'):
            raise ValueError('Argument "' + a + '" can\'t receive a negative value.')
    if args['skip_formatdb']:
        for n in range(len(args['genomes'])):
            if [str(n), 'db'] not in [f.split('.')[:2] for f in os.listdir(args['outpath']) if f.split('.')[-1] == 'nsq']:
                print 'BLAST databases not found'
                raise
    if args['outpath'][-1] != '/':
        args['outpath'] += '/'
    finder_args = deepcopy(args)
    finder_args['outfile'] = args['outpath'] + 'teste.fasta'
    finder_args['description'] = False
    finder_args['genome'] = args['genomes'][0]
    finder_args['circos'] = False
    finder_args['idist'] = 0
    al_finder.locus(finder_args)
    blast_args = deepcopy(args)
    blast_args['blast_database'] = []
    for n, infile in enumerate(args['genomes']):
        outfile = '{}.db'.format(n)
        log = args['outpath'] + '{}.formatdb.log'.format(n)
        blast_args['blast_database'].append(args['outpath'] + outfile)
        if not args['skip_formatdb']:
            al_formatdb.run_formatdb(infile, outfile, args['outpath'], log)
    blast_args['query'] = args['outpath'] + 'teste.fasta'
    blast_args['outfile'] = args['outpath'] + 'blasted.fasta'
    blast_args['blastm8'] = args['outpath'] + '*queryname**dbname*.m8'
    blast_args['log'] = args['outpath'] + '*dbname*.log' 
    blast_args['sum'] =  args['outpath'] + 'al_blast.sum'
    blast_args['border'] = 0
    for k in [k for k in blast_args.keys() if type(blast_args[k]) == str and 'queryname' in blast_args[k]]:
            blast_args[k] = blast_args[k].replace('*queryname*', blast_args['query'].split('/')[-1])
    al_blast.main(blast_args)
    align_args = deepcopy(args)
    align_args['join'] = False
    align_args['filter'] = False
    align_args['nexus'] = False
    align_args['chromo_sep'] = False
    align_args['idist'] = args['idist']
    align_args['sum'] =  args['outpath'] + 'al_blast.sum'
    align_args['dist_file'] = args['outpath'] + '/distances.txt'
    al_align.main(align_args)

if __name__ == '__main__':
    main()
