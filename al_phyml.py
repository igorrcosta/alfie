#!/usr/bin/env python
# -*- coding: utf-8 -*-
# al_phyml.py

__author__ = 'Igor Rodrigues da Costa'
__contact__ = 'igor.bioinfo@gmail.com'

import os
import shlex
import argparse
from subprocess import Popen
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC, generic_dna
from Bio import Phylo

def argument_parser(hlp = False):

    default_out = os.getcwd() + '/'
    parser = argparse.ArgumentParser(description = 'Runs phyml for all .aln files in a folder',\
                                     argument_default = None, fromfile_prefix_chars = '@')
    parser.add_argument('-i', '--inpath', nargs = '?', type = str, required = True,\
                        dest = 'inpath', help = 'Path to the folder with genbank sequences. (default: %(default)s)')
    parser.add_argument('-m', '--model', nargs = '?', type = str, default = 'HKY85',\
                        dest = 'model', help = 'Model to use in phyml. (default: %(default)s)')
    parser.add_argument('-p', '--protein', nargs = '?', const = True, default = False,\
                        dest = 'protein', help = 'Set this flag for protein sequences alignment and phylogeny. (default: %(default)s)')
    parser.add_argument('-t', '--topologies', nargs = '?', const = True, default = False,\
                        dest = 'topologies', help = 'Do topologies analysis. (default: %(default)s)')
    parser.add_argument('-c', '--concatenated', nargs = '?', const = True, default = False,\
                        dest = 'concatenated', help = 'Do a concatenated phylogeny (super-matrix). (default: %(default)s)')
    parser.add_argument('-g', '--gene_tree', nargs = '?', const = False, default = True,\
                        dest = 'gene_tree', help = 'Default individual AL phylogeny. (default: %(default)s)')
    if hlp:
        args = parser.parse_args(['-h'])
    else:
        args = parser.parse_args().__dict__
    return args

def main(args):
    if args['gene_tree']:
    	gene_trees(args)
    elif args['concatenated']:
        concatenated_tree(args)
    else:
	gene_trees(args)
    if args['topologies']:
        topologies(args)

def topologies(args):
    path = args['inpath']
    tree_files = [f for f in os.listdir(path) if ('tree.txt' in f) and 'all' not in f]
    hits = {'Pan':0, 'Gorilla':0, 'Pongo*': 0, 'Out':0}
    soma = len(tree_files)
    for f in sorted(tree_files):
        topology = 'Out'
        try:
	    tree = Phylo.read(path + f, 'newick')
	except:
	    print f
	    raise
        homo_clade = tree.find_clades('Homo').next()
        brother_clade = [c for c in tree.get_path(homo_clade)[-2] if c.name != 'Homo'][0].name
        if brother_clade == 'Pan':
            topology = '((H,C),G)'
        elif brother_clade == 'Gorilla':
            topology = '((H,G),C)'
        elif brother_clade == None:
            topology = '((C,G),H)'
	print f.split('.')[0] + '\t' + topology
        try:
            hits[brother_clade] += 1
        except:
            assert brother_clade == None
            hits['Out'] += 1
    print soma, 'trees'
    for k in hits.keys():
        print k + ':', hits[k], '('+str((hits[k]*100.0)/soma) + '%)'

def concatenated_tree(args):
    outpath = args['inpath']
    protein = args['protein']
    if args['model'] == 'HKY85' and protein:
        model = 'JTT'
    else:
        model = args['model']
    skip_phyml = False
    junta(outpath, protein)
    if protein:
        fastatophy(outpath + 'all_aa.aln', outpath + 'all_aa.phy')
        fastatophy(outpath + 'all_aa.aln', outpath + 'all_aa.nex', 'fasta', 'nexus')
    else:
        fastatophy(outpath + 'all_nuc.aln', outpath + 'all_nuc.phy')
        fastatophy(outpath + 'all_nuc.aln', outpath + 'all_nuc.nex', 'fasta', 'nexus', protein = False)
    command_nuc = 'phyml -m ' + model + ' -b 100 -v 0.0 -c 4 -a 4 -f m -i ' + outpath + 'all_nuc.phy'
    command_aa = 'phyml -d aa -m ' + model + ' -b 100 -v 0.0 -c 4 -a 4 -f m -i '+ outpath + 'all_aa.phy'
    try:
        a = open(outpath + 'log_phyml.txt', 'w')
        a.close()
    except:
       print 'Was not able to open log_phyml.txt. Check your permissions.'
       skip_phyml = True
    if skip_phyml:
       print 'Skiping phylogeny.'
       return None
    else:
        with open(outpath + 'log_phyml.txt', 'a') as log:
            if protein:
                log.write(command_aa + '\n')
                a = Popen(shlex.split(command_aa), stdout=log, stderr=log)
                a.wait()
            else:
                log.write(command_nuc + '\n')
                a = Popen(shlex.split(command_nuc), stdout=log, stderr=log)
                a.wait()

def gene_trees(args):
    path = args['inpath']
    if args['model'] == 'HKY85' and args['protein']:
        model = 'JTT'
    else:
        model = args['model']
    aln_genes = [f for f in os.listdir(path) if ('.nexus' in f and 'all' not in f and '.phy' not in f and '.aln' not in f)]
    for f in aln_genes:
	alnfile = path + f
	phyfile = path + f + '.phy'
        print f, alnfile
        if args['protein']:
            command = 'nohup phyml -d aa -m %s -b 100 -v 0.0 -c 4 -a 4 -f m -i %s'%(model, phyfile)
        else:
            command = 'nohup phyml -m %s -b 100 -v 0.0 -c 4 -a 4 -f m -i %s'%(model, phyfile)
        fastatophy(alnfile, phyfile, format_in = 'nexus')
        with open(path + 'log_phyml.txt', 'a') as log:
            log.write(command + '\n')
            a = Popen(shlex.split(command), stdout=log, stderr=log)
            a.wait()

def fastatophy(infile, outfile, format_in = 'fasta', format_out = 'phylip', protein = True):
    seq_records = []
    with open(infile, 'r') as handle:
        i = SeqIO.parse(handle, format_in)
        for seq in i:
            if format_out == 'nexus':
                if protein:
                    seq.seq.alphabet = IUPAC.protein
                else:
                    seq.seq.alphabet = IUPAC.unambiguous_dna
            if seq.id == 'Pongo':
                seq.id += '*'
            seq_records.append(seq)
    with open(outfile, 'wb') as out:
    	SeqIO.write(seq_records, out, format_out)

def junta(path, protein = False):
    end = '.aln'
    spec_dic = {}
    for f in os.listdir(path):
        if f.endswith(end):
            for seq in SeqIO.parse(path + f, 'fasta'):
                spec = seq.description.split('_')[0]
                if spec in spec_dic.keys():
                    spec_dic[spec].seq = spec_dic[spec].seq + seq.seq
                else:
                    spec_dic[spec] = SeqRecord(seq = Seq(str(seq.seq)), id = spec, description = '') #spec_dic = {especie1:str(gene1)+str(gene2), especie2: str(gene1)+str(gene2), ...}
    if protein:
        a = open(path + 'all_aa.aln', 'w')
        a.close()
        SeqIO.write(spec_dic.values(), path + 'all_aa.aln', 'fasta')
    else:
        a = open(path + 'all_nuc.aln', 'w')
        a.close()
        SeqIO.write(spec_dic.values(), path + 'all_nuc.aln', 'fasta')

if __name__ == '__main__':
    args = argument_parser()
    main(args)
