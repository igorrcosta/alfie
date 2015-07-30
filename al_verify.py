#!/usr/bin/env python
# -*- coding: utf-8 -*-
#al_verify.py
#08/2014

__author__ = 'Igor Rodrigues da Costa'
__contact__ = 'igor.bioinfo@gmail.com'

''' This program will check if ALs are really what they claim to be.

Input: Nexus file with ALs, location of the Blast databases
Output: AL-> Blast DB statistics.'''

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
    parser.add_argument('-b', '--blast_database', nargs = '*', type = str, required = True,\
                        help = 'Location of the blast databases.')
    parser.add_argument('-i', '--infile', nargs = '?', type = str, required = True,\
                        help = 'Nexus file with the ALs.')
    parser.add_argument('-o', '--outfile', nargs = '?', type = str, default = default_out,\
                        dest = 'outfile', help = 'File where the results will be saved.\n(default: %(default)s)')
    parser.add_argument('-l', '--log', nargs = '?', type = str, default = default_log,\
                        dest = 'log', help = 'Log file. (default: %(default)s)')
    parser.add_argument('-v', '--verbose', action = 'store_true', dest = 'verbose', help = 'Verbose switch.')

    args = parser.parse_args().__dict__
    return args

def main(args):
    '''Process the arguments and control the flow of the program.'''
 
    #Argument processing
    if args['verbose']:
        a = open(args['log'], 'w')
        a.close()
        def vprint(*a):
            # Print only when verbose.
            with open(args['log'], 'a') as log:
                log.write(' '.join((str(i) for i in a)))
    else:
        def vprint(*a):
            return None
    args['verbose'] = vprint
    vprint('ARGS: \n')
    for k in args:
        vprint(k, ' ', str(args[k]), '\n') #Log all arguments.
    vprint('\n')
    try:
        o = open(args['outfile'], 'w') #Statistics will be saved here.
        o.close()
    except:
        vprint('Can\'t create outfile.\n')
        raise

    results = {}
    if args['infile'].split('.')[-1] == 'nexus':
        files = nexus_parser(args['infile'])
	#files = [f[:-6] for f in os.listdir(os.getcwd()) if f.endswith('fasta')]
        for query in files:
            results[query] = blast(query + '.fasta', args) # results[query][db] = (evalue, id%)
    else:
	query = args['infile']
        results[query] = blast(query, args)
    statistics(results)

def nexus_parser(nexus_file):
    nexus = open(nexus_file, 'r')
    genes = {}
    contends = {}
    start = False
    for l in nexus:
        if 'matrix' in l:
            start = True
            l = nexus.next()
        while start:
            species = l.split()[0]
            seq = l.split()[1]
            contends[species] = seq
            l = nexus.next()
            if len(l.split()) < 2:
                start = False
        if len(l.split()) > 1:
            if l.split()[0] in contends.keys():
                contends[l.split()[0]] += l.split()[1]
            if l.split()[0] == 'CHARSET':
                genes[l.split()[1]] = (l.split()[3], l.split()[5][:-1])
    for gene in genes.keys():
        start = int(genes[gene][0]) - 1
        end = int(genes[gene][1])
        seqs = []
        for seq in contends.items():
            seqs.append(SeqRecord(Seq(seq[1][start:end], IUPAC.unambiguous_dna),
                                  id=seq[0], description=''))
        with open(gene + '.fasta', 'w') as outfile:
            SeqIO.write(seqs, outfile, 'fasta')
    return genes.keys()

def blast(query, args):
    '''Blast query against all dbs and analyse the resulting m8 files.'''
    results = {}
    for db in args['blast_database']: 
    	m8_name = query + db_name(db).split('.')[0] + '.m8' 
        run_blast(db, query, m8_name)
        results[db] = read_m8(m8_name, args, db) 
    return results

def db_name(db):
    '''Return the file name from a path.'''

    return db.split('/')[-1]
    
def run_blast(db, query, m8):
    '''Blast query file against db using default settings.'''

    command = 'blastn -task megablast -db '
    command += db
    command += ' -query '
    command += query
    command += ' -out '
    command += m8
    command += ' -outfmt 6 -num_threads 4 -evalue 0.01'
    command = ssplit(command)
    try:
        a = Popen(command)
        a.wait()
    except:
        print 'Something went wrong with blast!'
        argument_parser(hlp = True)
        raise

def read_m8(m8, args, db = ''):
    max_score = 0
    brest_hit = []
    with open(m8, 'r') as blast:
        for l in blast:
	    query = l.split()[0]
	    score = float(l.split()[11])
	    if score > max_score:
		max_score = score
		best_hit = [query]
	    elif score == max_score:
		best_hit.append(query)
    return best_hit 

def statistics(results):
    for gene, result in results.items():
        print gene
        for db, r in result.items():
	    correct = False
            for query in r:
		if query in db:
		    correct = True
	    if not correct:
		print gene, db_name(db), 'wrong predicted', r[0] 

if __name__ == '__main__':
    args = argument_parser()
    main(args)
