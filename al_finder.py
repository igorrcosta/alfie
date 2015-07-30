#!/usr/bin/env python
# -*- coding: utf-8 -*-
#al_finder.py
#12/2013

__author__ = 'Igor Rodrigues da Costa'
__contact__ = 'igor.bioinfo@gmail.com'

''' This program will find all Anonymous Regions from the genome in a sequential way,
and will split them in 1kb fragments.'''


import os
import csv
import argparse
from math import floor
from shlex import split as ssplit
from Bio import SeqIO


def argument_parser(hlp=False):
    '''al_finder.py'''

    default_out = os.getcwd() + '/al_finder.out'
    description_example = '1\t1450000\n2\t3204000\ncontig123\t400000'
    parser = argparse.ArgumentParser(description = 'Finds ALs', add_help = False,
                                     argument_default = None, fromfile_prefix_chars = '@',
                                     formatter_class = argparse.RawTextHelpFormatter)
    parser.add_argument("-h", "--help", action = "help", help = "Show this help message and exit.")
    parser.add_argument('-g', '--genome', nargs = '?', type = str, required = True,\
                        help = 'Fasta file with the genome.')
    parser.add_argument('-o', '--outfile', nargs = '?', type = str, default = default_out,\
                        dest = 'outfile', help = 'File where the ALs will be saved.\n(default: %(default)s)')
    parser.add_argument('-l', '--log', nargs = '?', type = str, default = 'al_finder.log',\
                        dest = 'log', help = 'Log file. (default: %(default)s)')
    parser.add_argument('--locus_length', nargs = '?', type = int, default = 2000,\
                        dest = 'length', help = 'Length of the ALs sequences.\n(default: %(default)s)')
    parser.add_argument('--max_n', nargs = '?', type = float, default = 0,\
                        dest = 'max_n', help = 'Maximum percentage of N\'s in the AL sequence.\n(default: %(default)s)')
    parser.add_argument('--inter_distance', nargs = '?', type = int, default = 10000,\
                        dest = 'idist', help = 'Minimum distance between ALs (negative means superposition).\n(defaut: %(default)s)')
    parser.add_argument('--gene_distance', nargs = '?', type = int, default = 10000,\
                        dest = 'gdist', help = 'Minimum (or maximum, if negative) distance between ALs and genes.\n(defaut: %(default)s)')
    parser.add_argument('--gene_locus', action = 'store_true', dest = 'gene_locus',\
                        help = 'Find coding regions loci.\n(defaut: %(default)s)')
    parser.add_argument('--cds', action = 'store_true', default = False,\
                        dest = 'cds', help = 'Only considers the CDS features of GTF files. (default: %(default)s)')
    parser.add_argument('--circos', action = 'store_true', default = False,\
                        dest = 'circos', help = 'Print anonymous regions: chr\tstart\tstop.\n(default: %(default)s)')
    parser.add_argument('--end_distance', nargs = '?', type = int, default = 10000,\
                        dest = 'edist', help = 'Distance between ALs and the start and end of a chromosome.\n(default: %(default)s)')
    parser.add_argument('-d', '--description', nargs = '?', type = str, default = None,\
                        dest = 'description', help = 'File with the id of all contigs to be analised\nand optionally their size in base pairs.\n' +
                        'Eg.: \n' + description_example)
    parser.add_argument('-e', '--est', nargs = '?', type = str, required = True,\
                        dest = 'est', help = 'File with all genome features coordinates. (GTF or INFO)')
    parser.add_argument('-v', '--verbose', action = 'store_true', dest = 'verbose', help = 'Verbose switch.')

    if hlp:
        args = parser.print_help()
    else:
        args = parser.parse_args().__dict__
    return args

def locus(args):
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
    locus_length = args['length'] #size of the putative ALs.
    idist = args['idist'] #distance between putative ALs. !!!
    max_n = floor(locus_length * args['max_n']) #args['max_n'] is the max % of N's
    try:
        o = open(args['outfile'], 'w')
        o.close()
    except:
        vprint('Can\'t open outfile: ' + args['outfile'] + '\n')
        raise
    #find the size of the chromosomes only if there is no description file.
    if args['description']:
        chromosome_dict = description_parser(args['genome'], args['description'], vprint) #fast
    else:
        chromosome_dict = chromosome_parser(args['genome'], vprint) #slow
    chromosomes = chromosome_dict.keys() #chromosome_dict = {chr:chr_size}
    chromosomes.sort() #ordered by name
    #find the coordinates of all GFF featues.
    try:
        est_list = gtf_parser(args['est'], args['cds'])
        vprint('\nGTF file found and parsed. \n')
    except:
	raise
        try:
            est_list = info_parser(args['est'])
            vprint('\nInfo file found and parsed.\n')
        except:
            vprint('EST file must be in GTF or info format.\n')
            raise argparse.ArgumentTypeError('File with genomic features must be of gtf or info format')
    #make a dictionary with the gene coordinates for every chromosome
    est_dict = est_info(est_list, chromosomes, vprint)
    #make a dictionary with the anonymous regions of every chromosome
    locus_dict = {}
    regions_total = 0
    for crom in chromosomes:
        est_pairs = est_dict[crom]
	if args['gene_locus']:
	    locus_dict[crom] = intersect(est_pairs)
	else:
            chromosome_coords = (args['edist'], chromosome_dict[crom] - args['edist'])
            if len(est_pairs) > 0:
                locus_dict[crom] = locus_finder(est_pairs, chromosome_coords, args['gdist'])
            elif args['gdist'] >= 0:
                locus_dict[crom] = [chromosome_coords]
	    else:
	        vprint('no genes in cromossome %s'%crom)
	        locus_dict[crom] = []
        regions_found = len(locus_dict[crom])
        regions_total += regions_found
        if not args['gene_locus']:
            vprint(str(regions_found) + ' anonymous regions found in chromosome ' + crom + '.\n')
        else:
	    vprint(str(regions_found) + ' coding regions found in chromosome ' + crom + '.\n')
    if not args['gene_locus']:
        vprint(str(regions_total) + ' total anonymous regions found.\n')
    else:
        vprint(str(regions_total) + ' total coding regions found.\n')
    if args['circos']:
        for chromo in locus_dict:
            for lr in locus_dict[chromo]:
                print '\t'.join(('hs' + str(chromo), str(lr[0]), str(lr[1])))
    #split the anonymous regions in several putative anonymous loci
    n = 1 #AL counter
    for c in SeqIO.parse(args['genome'], 'fasta'): #for every chr
        chromo = c.description.split()[0]
        try:
            locus_dict[chromo]
        except KeyError:
            if args['description']:
                vprint('Skipping chromosome "' + chromo + '":\nNot found in description file.')
                continue
            else:
                raise
        for locus_coords in locus_dict[chromo]: #for every anonymous region:
            count = 0
            start = locus_coords[0]
            end = start + locus_length
            while end <= locus_coords[1]:
                count += 1
                al = str(c.seq)[start:end] #grab a slice
                al = al.replace('n', 'N')
                if al.count('N') > max_n: #check number of Ns
                    if max_n == 0:
                        start += recursive_find(al, 'N', -1) + 1
                    else: #look for another slice in the next possible index
                        start += recursive_find(al, 'N', -max_n)
                    end = start + locus_length #look for another slice in the next possible index
                else:
                    with open(args['outfile'], 'a') as out: #save a pAL
                        out.write('>AL_' + str(n) + '|' + chromo + ':' + str(start) + ':' + str(end) + '\n' + al + '\n')
                        n += 1
                    start = end + idist #save the next slice index
                    end = start + locus_length

def est_info(est_list, chromosome_names, verbose = lambda *a, **b : None):
    '''Separate ESTs based on chromosome location.
    Returns a dictionary with chromosome name -> list of ESTs
    (start, stop) on that chromosome
    '''

    est_dict = {}
    est_not_in_genome = []
    for c in chromosome_names:
        est_dict[c] = []
    for est in est_list:
        if est[0] not in est_not_in_genome:
            try:
                est_dict[est[0]].append((int(est[1]), int(est[2])))
            except KeyError:
                est_not_in_genome.append(est[0])
                verbose('est: ' + str(est[0]) + ' is not in the genome.\n')
            except ValueError:
                verbose('est: ' + str(est) + ' is not int.\n')
                raise 
    return est_dict

def chromosome_parser(genome, verbose = lambda *a, **b : None):
    '''Returns a dictionary with chromosome name -> size.'''

    chromosome_dict = {}
    verbose('\nReading chromosomes\' sizes and names.\n')
    for chromo in SeqIO.parse(genome, 'fasta'):
        name = chromo.description.split()[0]
        size = len(chromo.seq)
        chromosome_dict[name] = size
        verbose('\t'.join([name, str(size)])+'\n')
    return chromosome_dict

def description_parser(genome, description, verbose = lambda *a, **b : None):
    chromosome_dict = {}
    failed_chromosomes = []
    verbose('\nReading chromosomes\' sizes and names.\n')
    with open(description, 'r') as desc:
        for l in desc:
            items = l.split()
            if len(items) > 1:
                try:
                    chromosome_dict[items[0]] = int(items[1])
                    verbose('\t'.join(items[:2])+'\n')
                except ValueError:
                    verbose('Description file needs an integer after the name of the chromosome.\n')
                    failed_chromosomes.append(items[0])
            elif len(items) == 1:
                failed_chromosomes.append(items[0])
    if len(failed_chromosomes) > 0:
        for chromo in SeqIO.parse(genome, 'fasta'):
            name = chromo.description.split()[0]
            if name in failed_chromosomes:
                size = len(chromo.seq)
                chromosome_dict[name] = size
                verbose('\t'.join([name, str(size)])+'\n')
                failed_chromosomes.pop(name)
                if len(failed_chromosomes) == 0:
                    break
    if len(failed_chromosomes) > 0:
        verbose('The following chromosomes were not found in the genome:\n')
        for f in failed_chromosomes:
            verbose(f + '\n')
        verbose('\n')
                
    return chromosome_dict
            
def info_parser(arq):
    '''Read INFO files.
    INFO is a tabular file with:
    <contig name> <initial pos> <end pos> <strand>'''

    with open(arq, 'rb') as f:
        spamReader = csv.reader(f, delimiter='\t')
        l = []
        for row in spamReader:
            l.append(row)
    return l

def gtf_parser(arq, gene = False):
    '''Read GTF files (ENSEMBL's General Transfer Format, with ESTs)'''

    with open(arq, 'rb') as f:
        handle = csv.reader(f, delimiter='\t')
        l = []
	if not gene:
            for row in handle:
                try:
                    l.append([row[0], row[3], row[4], row[6]])
                except:
                    if len(row) < 7:
                        continue
	else:
	    for row in handle:
		try:
	            if row[1] == 'protein_coding' and row[2] == 'gene':
		        l.append([row[0], row[3], row[4], row[6]])
                except:
                    if len(row) < 7:
                        continue
		    else:
			raise
    if not l:
       raise RuntimeWarning('No genes found on GTF file!') 
    return l

def gtf_to_info(arq):
    '''Transforms a gft file into an INFO file'''
    with open(arq, 'rb') as f:
        with open('INFO', 'wb') as w:
            read = csv.reader(f, delimiter='\t')
            write = csv.writer(w, delimiter='\t')
            for row in read:
                write.writerow([row[0], row[3], row[4], row[6]])

def locus_finder(est_pairs, genome_coords, distance):
    ''' Returns all anonimous loci regions coordinates.

    est_pairs is a list of (start, end) of all the ESTs
    genome_coords has the (start, end) of the contig/chromosome
    dist is the minimum (or maximum, if negative) distance for a locus to be considered anonimous

    >>>al = locus_finder( [(10,100), (5000,5010)], (0,10000), 1000)    
    >>>print al
    [(1101, 3999), (6011, 10000)]
    >>>for l in al:
        print l

    (1100, 4000)
    (6010, 10000)
    >>>
    '''

    exclude = [] # regions to be excluded from the genome
    dist = abs(distance)
    if distance < 0:
	dist += 1
    for est in est_pairs:
        est = tuple_sort(est)
        exclude.append((max(genome_coords[0], est[0] - dist), est[1] + dist)) 
    exclude.sort()
    exclude = intersect(exclude) #remove intersections
    locus = remove(exclude, genome_coords) #remove exclude coords from genome
    if distance < 0:
	exclude = locus[:] #exclude all regions too far from genes
        for est in est_pairs:
	    est = tuple_sort(est)
	    exclude.append(est) #exclude genes
	exclude.sort()
	exclude = intersect(exclude)
	locus = remove(exclude, genome_coords)
    if locus == [(0,0)]:
        return []
    else:
        return locus

def remove(set_list, big_set):
    '''Returns the set operation BIG_SET - SET_LIST,
    thinking of sets as coordinates, eg:
    (10,20) - [(9,11), (15,18), (25,32)] = [(12,14), (19,20)]'''

    set_list.sort()
    tuple_sort(big_set)
    ns = [big_set]
    for i in set_list:
        n = 0
        while n < len(ns):
            big_set = ns[n]
            if i[1] < big_set[0]:
                n += 1
                continue
            elif i[0] > big_set[1]:
                n += 1
                continue
            elif (big_set[0] <= i[0] <= big_set[1]) or \
                 (big_set[0] <= i[1] <= big_set[1]):
                ns.remove(big_set)
                ns += set_minus(big_set, i)
                n += 1
    return ns

def set_minus(a, b):
    '''Return a set (a - b),
    with sets representing (x,y) coordinates.
    eg: (1,5) - (2,3) = [(1,1), (4,5)]
    Sets must be ordered, run tuple_sort.
    Only to be used with integers.'''

    if b[1] < a[0] or b[0] > a[1]:
        return [a]
    elif b[0] == a[0]:
        if b[1] >= a[1]:
            return []
        else:
            return [(b[1] + 1, a[1])]
    else: #b[0] > a[0]
        if b[1] >= a[1]:
            return [(a[0], b[0] - 1)]
        else:
            return [(a[0], b[0] - 1), (b[1] + 1, a[1])]

def intersect(set_list):
    '''Returns a list of sets, duplicated elements removed.
    >>>intersect([(0,3), (2,4), (2,4), (5,7), (6,6)])
    [(0,4), (5,7)]
    '''
    sl = set_list
    nl = []
    for n in xrange(len(set_list)):
        if len(nl) == 0:
            nl.append(sl[n])
        elif sl[n][1] <= nl[-1][1]:
            continue
        elif sl[n][0] > nl[-1][1]:
            nl.append(sl[n])
            continue
        elif sl[n][1] > nl[-1][1]: #[(0,1) (0,2)]
            old = nl.pop()
            new = (old[0], sl[n][1])
            nl.append(new)

    return nl

def tuple_sort(s):
    '''Sort-a-tuple, smaller element first.'''

    if s[0] > s[1]:
        return(s[1], s[0])
    return s

def recursive_find(string, element, n):
    '''find index of the 'n'th occurence of 'element' in 'string'.'''
    start = 0
    end = len(string)
    if string.find(element, start, end) == -1:
        return -1
    if n > 0:
        while n > 0:
            last = string.find(element, start, end) #index of the first occurence
            start = last + 1
            n = n - 1
    elif n < 0:
        string = string[::-1]
        while n < 0:
            last = string.find(element, start, end) #index of the last occurence
            start = last + 1
            n = n + 1
        if last != -1:
            last = end - last - 1
    elif n == 0:
        raise ValueError('Nth can\'t be 0.')
    return last  

if __name__ == "__main__":
    args = argument_parser()
    for a in args:
        if isinstance(args[a], int) and args[a] < 0 and a not in ('idist', 'gdist'):
            raise ValueError('Argument "' + a + '" can\'t receive a negative value.')
    l = locus(args)
