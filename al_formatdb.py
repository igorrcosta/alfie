#!/usr/bin/env python
# -*- coding: utf-8 -*-
#al_formatdb.py
#05/2015

__author__ = 'Igor Rodrigues da Costa'
__contact__ = 'igor.bioinfo@gmail.com'

''' This program run formatdb on a fasta file.'''

import os
import argparse
from subprocess import Popen
from shlex import split as ssplit

def argument_parser():
    '''al_formatdb.py'''

    default_out = os.getcwd()
    parser = argparse.ArgumentParser(description = 'Run formatdb', argument_default = None,\
                                     fromfile_prefix_chars = '@', add_help = False,\
                                     formatter_class = argparse.RawTextHelpFormatter)
    parser.add_argument("-h", "--help", action = "help", help = "Show this help message and exit.")
    parser.add_argument('-i', '--infile', nargs = '?', type = str, required = True,\
                        dest = 'infile', help = 'Genome fasta file.')
    parser.add_argument('-o', '--outfile', nargs = '?', type = str,\
                        dest = 'outfile', help = 'Database title.\n(default: Genome name)')
    parser.add_argument('-p', '--outpath', nargs = '?', type = str, default = default_out,\
                        dest = 'outfile', help = 'Database title.\n(default: %(default)s)')
    parser.add_argument('-l', '--logfile', nargs = '?', type = str,\
                        dest = 'logfile', help = 'Log file.\n(default: /dev/null)')

    args = parser.parse_args().__dict__
    return args

def run_formatdb(fasta_file, out_file, out_path, log_file='/dev/null'):
    '''Blast query file against db using default settings.'''
    
    fasta_file = os.path.abspath(fasta_file)
    log_file = os.path.abspath(log_file)
    old_dir = os.getcwd()
    os.chdir(out_path)
    command = 'formatdb -p F -i '
    command += fasta_file
    command += ' -l ' 
    command += log_file
    if out_file:
        command += ' -n  '
        command += out_file
    command = ssplit(command)
    try:
        a = Popen(command)
        a.wait()
        os.chdir(old_dir)
    except:
        os.chdir(old_dir)
        print 'Something went wrong with formatdb!'
        raise


if __name__ == '__main__':
    args = argument_parser()
    run_formatdb(args['infile'], args['outfile'], args['outpath'])
