#!/usr/bin/env python
# -*- coding: utf-8 -*-
#al_chromosome.py

''' Return chromosome location AL statistics.'''

__author__ = 'Igor Rodrigues da Costa'
__contact__ = 'igor.bioinfo@gmail.com'

from sys import argv

def get_chromo(al):
    chromo = {}
    with open(al, 'r') as a:
        for l in a:
            if l[0] == '>':
                c = l.split(':')[0].replace('>', '')
                #c = l.split('|')[1].split(':')[0]
                try:
                    chromo[c] += 1
                except KeyError:
                    chromo[c] = 1
    keys = []
    for k in chromo.keys():
        try:
            keys.append(int(k))
        except:
            keys.append(k)
    keys.sort()
    for k in keys:
        print k, chromo[str(k)]


if len(argv) > 1:
    get_chromo(argv[1])
else:
    print './al_chromosome.py al_file'
