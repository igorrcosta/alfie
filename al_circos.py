#!/usr/bin/env python
# -*- coding: utf-8 -*-
#al_circos.py

''' Return anonymous regions coordinates for circos.'''

__author__ = 'Igor Rodrigues da Costa'
__contact__ = 'igor.bioinfo@gmail.com'

from sys import argv

def get_regions(al_file):
    #Concatenate all contiguous putative loci from fasta and write in regions.txt file.
    circos = []
    chromo = '0'
    start = '0'
    end = '0'
    with open(al_file, 'r') as a:
        for l in a:
            if l[0] == '>':
                chromo_new = l.split(':')[0].split('|')[-1]
                start_new = l.split(':')[1]
                end_new = l.split(':')[2][:-1]
                #Concatenate contiguous regions. 
                if chromo == chromo_new and start_new == end:
                    print l
                    print circos.pop()
                    end = end_new
                else:
                    start = start_new
                    end = end_new
                    chromo = chromo_new
                circos.append('\t'.join(['hs'+chromo, start, end]))
    with open('regions.txt', 'w') as regions_file:
        for c in circos:
            regions_file.write(c + '\n')

def circos(al_coords):
    #Write AL coordinates in al_circos.txt file.
    with open(al_coords, 'r') as coords, open('al_circos.txt', 'w') as outfile:
        for l in coords:
            al_id = l.split()[0]
            al_chromo = l.split()[1]
            al_start = l.split()[2]
            al_end = l.split()[3]
            outfile.write('\t'.join(['hs'+al_chromo, al_start, al_end])+'\n')

if len(argv) > 1:
    circos(argv[1])
    get_regions(argv[2])
else:
    print './al_circos.py al_coords pl_coords.fasta'

