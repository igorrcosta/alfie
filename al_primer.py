#!/usr/bin/env python
# -*- coding: utf-8 -*-
#al_primer.py

''' This is an interface to EMBOSS' eprimer3 interface of primer3.'''

__author__ = 'Igor Rodrigues da Costa'
__contact__ = 'igor.bioinfo@gmail.com'

import os
import argparse
from subprocess import Popen
import subprocess
import shlex
from time import sleep

def argument_parser(hlp = False):
    '''al_primer.py -i seq.fasta
    Output default: "./seq.out".'''

    default_out = os.getcwd() + '/seq.out'
    parser = argparse.ArgumentParser(description = 'Finds primers',\
                                     argument_default = None, fromfile_prefix_chars = '@')
    parser.add_argument('-i', '--sequence', nargs = '?', type = str, required = True,\
                        dest = 'sequence', help = 'Fasta file with the sequence targets. (default: %(default)s)')
    parser.add_argument('-o', '--outfile', nargs = '?', type = str, default = default_out,\
                        dest = 'outfile', help = 'File where the primers will be saved. (default: %(default)s)')
    parser.add_argument('-n', '--numreturn', nargs = '?', type = int, default = 1,\
                        dest = 'numreturn', help = 'The maximum number of primer pairs to return. (default: %(default)s)')
    parser.add_argument('--mintm', nargs = '?', type = float, default = 50.0,\
                        dest = 'mintm', help = 'Minimum acceptable melting temperature for a primer oligo. (default: %(default)s)')
    parser.add_argument('--maxtm', nargs = '?', type = float, default = 65.0,\
                        dest = 'maxtm', help = 'Maximum acceptable melting temperature for a primer oligo. (default: %(default)s)')
    parser.add_argument('--mingc', nargs = '?', type = int, default = 20,\
                        dest = 'mingc', help = 'Minimum acceptable gc percentage. (default: %(default)s)')
    parser.add_argument('--maxgc', nargs = '?', type = int, default = 80,\
                        dest = 'maxgc', help = 'Maximum acceptable gc percentage. (default: %(default)s)')
    parser.add_argument('--minsize', nargs = '?', type = int, default = 18,\
                        dest = 'minsize', help = 'Minimum acceptable length of a primer. (default: %(default)s)')
    parser.add_argument('--maxsize', nargs = '?', type = int, default = 27,\
                        dest = 'maxsize', help = 'Maximum acceptable length (in bases) of ' +\
                        'a primer. Currently this parameter cannot be larger than 35. (default: %(default)s)')
    parser.add_argument('--prange', nargs = '?', type = str, default = '500-600',\
                        dest = 'prange', help = 'Product size range. Smaller range size increases speed considerably. (default: %(default)s)')
    parser.add_argument('--psizeopt', nargs = '?', type = int, default = 550,\
                        dest = 'psizeopt', help = 'The optimum size for the PCR product. '+\
                        '0 indicates that there is no optimum product size. It must be an Integer >= 0. (default: %(default)s)')
    parser.add_argument('-f', '--format', nargs = '?', choices=['fasta', 'list', 'tab'], default = 'list',\
                        dest = 'format', help = 'Print the output in fasta format (fasta), one primer sequence per line (list), one primer pair per line (tab). (default: %(default)s)')
    parser.add_argument('--ftail', nargs = '?', const = True, default = False,\
                        dest = 'ftail', help = 'Add a M13Forward primer "tail" '+\
                        'to the forward anonymous locus primer. \n'+\
                        'The M13F-21 primer sequence is:   GTTGTAAAACGACGGCCAGT. (default: %(default)s)')
    parser.add_argument('--rtail', nargs = '?', const = True, default = False,\
                        dest = 'rtail', help = 'Add a M13Reverse primer "tail" '+\
                        'to the reverse anonymous locus primer. \n'+\
                        'The M13R-29 primer sequence is:  CACAGGAAACAGCTATGACC. (default: %(default)s)')
    parser.add_argument('--targetregion', nargs = '?', type = str,\
                        dest = 'targetregion', help = 'Sequence region that must be included in the product. \n'+\
                        'Values in the form "start,end", where "start" and "end" are the position of the first and last base to be included.')
    if hlp:
        args = parser.parse_args(['-h'])
    else:
        print 'Welcome to al_primer!'
        args = parser.parse_args().__dict__
        args['numreturn'] = args['numreturn'] + 1
    return args

def primer_run(args):
    'eprimer3 shell'
    command = 'eprimer3'
    args = validate(args)
    for arg, value in args.iteritems():
        print arg, value
        if value is not False:
            command += ' -'
            command += arg
            command += ' '
            command += str(value)
        if value is True and value not in ['ftail', 'rtail']:
            command += ' -'
            command += arg
    command = shlex.split(command)
    clean = open(args['outfile'], 'w')
    clean.close()
    try:
        a = Popen(command)
        a.wait()
    except:
        print 'Something went wrong with eprimer3!'
        argument_parser(hlp = True)
        raise

def validate(args):
    'Check if sequence size is compatible with the parameters'
    min_size = int(args['prange'].split('-')[0])
    max_size = args['prange'].split('-')[1]
    with open(args['sequence'], 'r') as f:
       for line in f:
           seq = ''
           if line[0] == '>':
               line = next(f)
               while line[0] != '>':
                   seq += line[:-1]
                   try:
                       line = next(f)
                   except:
                       break
               if len(seq) < min_size:
                   args['prange'] =  str(args['maxsize']) + '-' + max_size #max primer size is the lower limit for product size in primer3
                   print 'Minimun PRANGE parameter is too big. Adjusting to fit seq lenght'
                   print 'New PRANGE is %s'%args['prange']
           break
    return args

def primer3_parser(outfile):
    '''yield [forward/primer/start/len/tm/gc/seq/n]'''
    with open(outfile, 'r') as o:
        for result in o:
            result = result.strip()
            if result.startswith('# EPRIMER3'):
                for i in range(10):
                    foward = next(o, '').split()
		    if len(foward) == 7:
			print foward
		        break
                #foward = next(o, '').split()
                next(o, '')
                reverse = next(o, '').split()
                if len(reverse) == 7:
                    yield(foward, reverse)
                else:
		    print len(foward), len(reverse)
                    yield None

def primer2ispcr(outfile, infile, ftail, rtail, recursive = False):
    # Id \t foward_sequence \t reverse_sequence
    tmp = open(outfile[:-4], 'w')
    tmp.close()
    for (p, i) in zip(primer3_parser(outfile), getfastaids(infile)):
        if not p:
            if recursive:
                with open(outfile[:-4], 'a') as fl:
                    fl.write(i + ' not found \n')
            continue
        elif p[0][0] == 'FORWARD':
            fseq = p[0][6]
            if ftail:
                fseq = 'GTTGTAAAACGACGGCCAGT' + fseq
            rseq = p[1][6]
            if rtail:
                rseq = 'CACAGGAAACAGCTATGACC' + rseq
            record = i.replace('\n', '') + '\t' + fseq + '\t' + rseq + '\n'
            with open(outfile[:-4], 'a') as fl:
                fl.write(record)
    #Popen(shlex.split('rm ' + outfile))

def primer2list(outfile, infile, ftail, rtail):
    srecord_list = [] # Id_(F/R) \t sequence
    for n, (pp, i) in enumerate(zip(primer3_parser(outfile), getfastaids(infile))):
        for p in pp: #(foward, reverse)
            seq = p[6]
            if p[0] == 'FORWARD':
                if ftail:
                    seq = 'GTTGTAAAACGACGGCCAGT' + seq
                else:
                    simple_id = i+'_F'
            elif p[0] == 'REVERSE':
                if rtail:
                    seq = 'CACAGGAAACAGCTATGACC' + seq
                simple_id = i+'_R'
            assert seq
            simple_record = simple_id + '\t' + seq + '\n'
            srecord_list.append(simple_record)
    with open(outfile[:-4], 'w') as fl:
        print len(srecord_list), 'sequences had their primers done.'
        for r in srecord_list:
            fl.write(r)
    #Popen(shlex.split('rm ' + outfile))

def primer2fasta(outfile, infile, ftail, rtail):
    record_list = [] #Fasta format: >id|n|(forward or reverse)|tm|gc|start-end \n sequence
    for n, (pp, i) in enumerate(zip(primer3_parser(outfile), getfastaids(infile))):
        for p in pp: #(foward, reverse)
            seq = p[6]
            pos = str((n+2)/2)
            start = p[2]
            end = p[3] + start
            start_end = start + '-' + end
            tm = 'tm:' + p[4]
            gc = 'gc:' + p[5]
            if p[0] == 'FORWARD':
                if ftail:
                    seq = 'GTTGTAAAACGACGGCCAGT' + seq
                fend = float(end)
                fasta_id = '|'.join([i, pos, 'forward', tm, gc, 'Start-End:'])
            elif p[0] == 'REVERSE':
                if rtail:
                    seq = 'CACAGGAAACAGCTATGACC' + seq
                product_size = float(start) - fend
                fasta_id = '|'.join([i, pos, 'reverse', tm, gc, 'Start-End:']) + start_end
                fasta_id = fasta_id + '|Product_Size:' + str(product_size)
            assert start and end and seq
            record = '>'+fasta_id+'\n'+seq+'\n'
            record_list.append(record)

    with open(outfile[:-4], 'w') as fl:
        print len(record_list), 'sequences had their primers done.'
        for r in record_list:
            fl.write(r)
    #Popen(shlex.split('rm ' + outfile))

def getfastaids(infile):
    with open(infile) as i:
        ids = []
        for seq in i.readlines():
            if seq[0] == '>':
                s = seq.split('|')
                if len(s) > 1:
                    ids.append(s[0][1:])
                else:
                    ids.append(seq[1:])
    return ids

if __name__ == "__main__":
    args = argument_parser()
    args['outfile'] = args['outfile'] + '.tmp'
    formating = args.pop('format')
    primer_run(args)
    if formating == 'fasta':
        primer2fasta(args['outfile'], args['sequence'], args['ftail'], args['rtail'])
    elif formating == 'list':
        primer2list(args['outfile'], args['sequence'], args['ftail'], args['rtail'])
    elif formating == 'tab':
        primer2ispcr(args['outfile'], args['sequence'], args['ftail'], args['rtail'], recursive = True)
