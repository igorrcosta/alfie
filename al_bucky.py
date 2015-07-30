#!/usr/bin/env python
# -*- coding: utf-8 -*-
#al_verify.py
#08/2014

__author__ = 'Igor Rodrigues da Costa'
__contact__ = 'igor.bioinfo@gmail.com'

'''Runs bucky on several nexus files.'''

import os
import argparse
from subprocess import Popen
from shlex import split as ssplit
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord 
from Bio.Alphabet import IUPAC

def argument_parser():
    '''al_verify.py'''

    default_log = os.getcwd() + '/al_verify.log'
    default_out = os.getcwd() + '/al_verify.out'

    parser = argparse.ArgumentParser(description = 'Check ALs', argument_default = None,\
                                     fromfile_prefix_chars = '@', add_help = False,\
                                     formatter_class = argparse.RawTextHelpFormatter)
    parser.add_argument("-h", "--help", action = "help", help = "Show this help message and exit.")
    parser.add_argument('-i', '--infile', nargs = '?', type = str, required = True,\
                        help = 'Sample file listing the nexus to be used.')
    parser.add_argument('--infolder', nargs = '?', type = str, required = True,\
                        help = 'Folder with the nexus files.')
    parser.add_argument('-o', '--outfolder', nargs = '?', type = str, default = default_out,\
                        dest = 'outfolder', help = 'Folder where the files will be saved.\n(default: %(default)s)')
    parser.add_argument('-l', '--log', nargs = '?', type = str, default = default_log,\
                        dest = 'log', help = 'Log file. (default: %(default)s)')
    parser.add_argument('-v', '--verbose', action = 'store_true', dest = 'verbose', help = 'Verbose switch.')

    args = parser.parse_args().__dict__
    return args

def main(args):
    '''Process the arguments and control the flow of the program.'''
    if args['infolder'][-1] != '/':
	args['infolder'] = args['infolder'] + '/'
    nexus_list = []
    with open(args['infile'], 'r') as sample:
        nexus_list = [f[:-1] for f in sample.readlines()]
    mb_command = '''

outgroup 4;
mcmc ngen=500000 nrun = 2 burnin=100000 nchain = 2 samplefreq=1000 ;
sumt;
sump;
;
'''
    outfiles = []
    for f in nexus_list:
	with open(args['infolder'] + f, 'r') as nexus:
            outfile = args['outfolder'] + 'bucky.' + f
	    with open(outfile, 'w') as out:
		out.write(nexus.read())
		out.write(mb_command)
	    outfiles.append(outfile)
    for f in outfiles:
	command = ssplit('mb ' + f)
        a = Popen(command)
        a.wait()
	
if __name__ == '__main__':
    args = argument_parser()
    main(args)
