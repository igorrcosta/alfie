#!/usr/bin/env python
# -*- coding: utf-8 -*-
#al_blast.py
#01/2014

__author__ = 'Igor Rodrigues da Costa'
__contact__ = 'igor.bioinfo@gmail.com'

''' This program will check if ALs are unique and present
in the genomes of compared species. Blast+ required.

todo: melhorar documentação, adicionar .format as strings, melhorar erros.

Input: file with ALs, location of the Blast databases (self/non-self)
Optional input: size of the border, cutoff values (identity%, coverage%).
Output: ids or fasta (default) of the ALs that were found to be present and unique.'''


import os
import argparse
from subprocess import Popen
from shlex import split as ssplit
import uuid
import pickle
from tempfile import NamedTemporaryFile

def argument_parser():
    '''al_blast.py'''

    default_out = os.getcwd() + '/*queryname*.out'
    default_log = os.getcwd() + '/*queryname*.log'
    default_m8 = os.getcwd() + '/*queryname**dbname*.m8'
    default_sum = os.getcwd() + '/al_blast.sum'
    parser = argparse.ArgumentParser(description = 'Blasts ALs', argument_default = None,\
                                     fromfile_prefix_chars = '@', add_help = False,\
                                     formatter_class = argparse.RawTextHelpFormatter)
    parser.add_argument("-h", "--help", action = "help", help = "Show this help message and exit.")
    parser.add_argument('-b', '--blast_database', nargs = '*', type = str, required = True,\
                        help = 'Location of the blast databases.')
    parser.add_argument('-q', '--query', nargs = '?', type = str, required = True,\
                        help = 'Fasta file with the sequences to be blasted.')
    parser.add_argument('-o', '--outfile', nargs = '?', type = str, default = default_out,\
                        dest = 'outfile', help = 'File where the selected ALs fasta will be saved.\n(default: %(default)s)')
    parser.add_argument('-m', '--blastm8', nargs = '?', type = str, default = default_m8,\
                        help = 'File where the blast output will be saved.\n(default: %(default)s)')
    parser.add_argument('-l', '--log', nargs = '?', type = str, default = default_log,\
                        dest = 'log', help = 'Log file. (default: %(default)s)')
    parser.add_argument('-s', '--sum', nargs = '?', type = str, default = default_sum,\
                        dest = 'sum', help = 'File to save a tabular sumary of the result. (default: %(default)s)')
    parser.add_argument('--duplication_cutoff', nargs = '?', type = int, default = 50,\
                        dest = 'dup_cut', help = 'ALs with 2 hits with identity higher than this will be considered duplicated.\n(default: %(default)s)')
    parser.add_argument('--identity_cutoff', nargs = '?', type = int, default = 90,\
                        dest = 'id_cut', help = 'ALs with a identity higher than this will be considered homologous.\n(default: %(default)s)')
    parser.add_argument('--coverage_cutoff', nargs = '?', type = int, default = 90,\
                        dest = 'cov_cut', help = 'BLAST hits must have at least this much %%coverage to be considered hits.\n(default: %(default)s)')
    parser.add_argument('--border_size', nargs = '?', type = int, default = 0,\
                        dest = 'border', help = 'Size of the start and end margin to be excluded from the blast search.\n(default: %(default)s)')
    parser.add_argument('--megablast', action = 'store_true', dest = 'megablast', help = 'Run blastn with the "-task" flag set to megablast.')
    parser.add_argument('-v', '--verbose', action = 'store_true', dest = 'verbose', help = 'Verbose switch.')

    args = parser.parse_args().__dict__
    return args

def main(args):
    '''Process the arguments and control the flow of the program.'''
 
    #Argument processing
    args['log'] = os.path.abspath(args['log'])
    args['outfile'] = os.path.abspath(args['outfile'])
    #args['query'] = os.path.abspath(args['query'])
    args['blast_database'] = [os.path.abspath(arg) for arg in args['blast_database']]
    args['sum'] = os.path.abspath(args['sum'])
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
    vprint('ARGS: \n')
    for k in args:
        vprint(k, ' ', str(args[k]), '\n') #Log all arguments.
    vprint('\n')
    try:
        o = open(args['outfile'], 'w') #Final ALs will be saved here (fasta).
        o.close()
    except:
        vprint('Can\'t create outfile.\n')
        raise
    #Initiating variables
    excluded_als = []
    excluded_last_run = 0
    best_hits = {} #best_hits = {genome_db:{al_id:(subj_name, subj_start, subj_end)}}
    for db in args['blast_database']: #Main loop, run for every genome.
        #Save all ALs, except those in excluded_als list, in the temp file.
        all_ids, query_files = create_query(args['query'], args['border'], excluded_als, vprint=vprint) 
        total_ids = (len(all_ids) - 1) * 1000 + len(all_ids[-1])
        vprint(str(total_ids) + ' putative ALs.\n Runing blast against ' + db + '\n')
        #Run megablast using the temp query against "db".
        for splice_ids, query in zip(all_ids, query_files):
            run_blast(db, query.name, args['blastm8'], megablast=False, vprint=vprint)
            #Update excluded_als with all ALs that failed id, coverage or duplication filters. Save best hits.
            excluded_als, best_hits[db] = read_m8(args, splice_ids, excluded_als, vprint, db) 
            os.remove(query.name) #Deletes temporary query file.
        vprint(str(len(excluded_als) - excluded_last_run) + ' ALs have been excluded against ' + db.split('/')[-1] + ' database.\n')
        excluded_last_run = len(excluded_als)
    #Writing output
    final_als = create_output(args, excluded_als) #Save the final filtered set of ALs in fasta format.
    create_sum(args, best_hits, final_als) #save a tabular file with a sumary of the blast searchs.
    return final_als, best_hits

def create_sum(args, best_hits, final_als):
    '''write the start and end of als found in a .sum file. 
    example for 2 genomes (db_name1 and db_name2):
    
    al_id   bd_name1    subj_name1   subj_start  subj_end
    al_id   bd_name2    subj_name1   subj_start  subj_end
    al_id2  bd_name1    subj_name2   subj_start  subj_end
    al_id2  bd_name2    subj_name2   subj_start  subj_end'''
    
    sumfile = args['sum']
    with open(sumfile, 'w') as sumf:
        for al in final_als:
            for db in sorted(best_hits.keys()):
                sumf.write('\t'.join((str(al), db_name(db)) + best_hits[db][al]) + '\n')
    
def db_name(db):
    '''return the file name from a path.'''

    return db.split('/')[-1]
    
def create_output(args, excluded):
    '''copy the contents of query file to the output file, excludind all ids that match the excluded list.
     Returns a list with all selected al ids.'''

    final_als = []
    with open(args['query'], 'r') as ff:
        with open(args['outfile'], 'w') as out_file:
            for name, seq, none in readfq(ff):
                al = name.split('|')[0].split('_')[1]
                if al in excluded:
                    continue
                out_file.write('>' + name + '\n' + seq + '\n')
                final_als.append(al)
    return final_als

def create_query(fasta_file,  border = 0, exclude = [], vprint = lambda x: None):
    '''Creates a temporary query fasta file from fasta_file seqs, return the list of query ids.'''

    all_ids = [[]]
    tmp_files = []
    try:
        tmp_file = NamedTemporaryFile(delete=False) #Temporary query with unique filename, will be deleted at the end.
        tmp_file.close()
    except:
        vprint('Can\'t create temporary file.\n')
        raise
    saved_fasta = 0 # Change files every 1k seq.
    with open(fasta_file, 'r') as ff:
        for name, seq, none in readfq(ff):
            try:
                n = name.split('|')[0].split('_')[1]
            except:
                print(name + ' is not well formated.')
                continue
            if n not in exclude:
                with open(tmp_file.name, 'a') as query:
                    query.write('>' + name + '\n')
                    if border > 0:
                        query.write(seq[border:-border] + '\n')
                    else:
                        query.write(seq + '\n')
                    all_ids[-1].append(n)
                    saved_fasta += 1
                if not saved_fasta % 1000:
                    tmp_files.append(tmp_file)
                    tmp_file = NamedTemporaryFile(delete=False) #Temporary query with unique filename, will be deleted at the end.
                    tmp_file.close()
                    all_ids.append([])
    if tmp_file not in tmp_files:
        tmp_files.append(tmp_file)
    assert len(all_ids) == len(tmp_files)
    if len(all_ids) > 1:
        assert len(all_ids[0]) == 1000
    assert len(all_ids[-1]) <= 1000
    return all_ids, tmp_files
    
def readfq(fp):
    '''This is fasta/q parser generator function.'''

    last = None
    while True:
        if not last:
            for line in fp:
                if line[0] in '>@': # Fasta/q header line
                    last = line[:-1] # Save this line
                    break
        if not last: break # No more headers
        name, seqs, last = last[1:], [], None
        for l in fp: # Read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # This is a fasta record
            yield name, ''.join(seqs), None # Yield last fasta record
            if not last: break
        else: # This is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # Read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # Have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # Yield a fastq record
                    break
            if last: # Reach EOF before reading enough quality
                yield name, seq, None # Yield a fasta record instead
                break
            
def run_blast(db, query, m8, megablast=True, vprint=lambda x: None):
    '''Blast query file against db using default settings.'''

    if '*dbname*' in m8:
        m8 = m8.replace('*dbname*', db_name(db))
    command = 'blastn'
    if megablast:
        command += ' -task megablast'
    command += ' -db '
    command += db
    command += ' -query '
    command += query
    command += ' -out '
    command += m8
    command += ' -outfmt 6 -num_threads 6 -evalue 0.01'
    vprint('blast command: ', command)
    command = ssplit(command)
    try:
        a = Popen(command)
        a.wait()
    except OSError:
        print('blastn not found!')
        raise

def read_m8(args, all_ids, excluded = [], vprint = lambda x: None, db = ''):
    '''Return a list of queries with lower than cutoff identity or duplicated and a dictionary with the best hits.'''
    
    found = []
    duplicated = []
    low_id = []
    low_coverage = []
    total = 0
    m8 = args['blastm8']
    best_hits = {} #best_hits[id] = (subj_name, subj_start, subj_end)
    if '*dbname*' in m8 and db:
        m8 = m8.replace('*dbname*', db_name(db))
    with open(m8, 'r') as blast:
        for l in blast:
            total += 1
            n = l.split('|')[0].split('_')[1] #if "l = >AL_17|1:2221075:2223075", n = 17 or 17.2 if UCE
            identity = float(l.split()[2])
            subj_start = int(l.split()[8])
            subj_end = int(l.split()[9]) 
            coverage = abs(subj_end - subj_start)
            subj_name = str(l.split()[1])
            query_size = int(l.split('|')[1].split(':')[2].split()[0]) -  int(l.split('|')[1].split(':')[1])
            #print("n: ", n, "identity: ", identity, "start: ", subj_start, "coverage: ", coverage, "name: ", subj_name, "q size: ", query_size)
            #print(n not in excluded, query_size*100.0/coverage >= args['cov_cut'])
            if n not in excluded and query_size*100.0/coverage >= args['cov_cut']: # not already excluded and coverage above cutoff
                if identity > args['id_cut']: #above identity cutoff, may be duplicated...
                    if n in low_coverage:
                        low_coverage.remove(n) #found a hit with high coverage
                    if n in found: #excluded due to duplication (both above identity cutoff)
                        excluded.append(n)
                        duplicated.append(n)
                        best_hits.pop(n)
                    elif n in low_id: #excluded due to another hit above duplication cutoff
                        excluded.append(n)
                        duplicated.append(n)
                        low_id.remove(n)
                    else:
                        found.append(n)
                        best_hits[n] = (subj_name, str(subj_start), str(subj_end))
                elif identity > args['dup_cut']: #above duplication cutoff, might exclude another AL 
                    if n in low_coverage:
                        low_coverage.remove(n) #found a hit with high coverage
                    if n in found: #exclude the high id AL found due to duplication
                        duplicated.append(n)
                        excluded.append(n)
                        best_hits.pop(n)
                    elif n not in low_id:
                        low_id.append(n) #there is a hit above duplication cutoff
                elif query_size*100.0/coverage < args['cov_cut']: #not in excluded, but below coverage cutoff
                    low_coverage.append(n)
        
    not_found = []
    for i in all_ids:
        if i not in found:
            excluded.append(i)
            not_found.append(i)
    vprint(str(len(low_coverage)) + ' excluded due to low coverage.\n')
    vprint(str(len(low_id)) + ' excluded due to low identity.\n')
    vprint(str(len(duplicated)) + ' excluded due to duplication.\n')
    vprint(str(len(not_found)) + ' excluded due to not being found.\n')
    vprint(str(len(excluded)) + ' total excluded. \n')
    return excluded, best_hits

if __name__ == '__main__':
    args = argument_parser()
    for k in [k for k in args.keys() if type(args[k]) == str and 'queryname' in args[k]]:
        args[k] = args[k].replace('*queryname*', args['query'].split('/')[-1])
    locus, dict_hits = main(args)
#    with open('list.txt', 'w') as l:
#        pickle.dump((locus, dict_hits), l)
