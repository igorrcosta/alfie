#!/usr/bin/env python
# -*- coding: utf-8 -*-
#nexus2phylip.py
#3/2014

__author__ = 'Igor Rodrigues da Costa'
__contact__ = 'igor.bioinfo@gmail.com'

''' This program will get ALs from a nexus file and make a phylip file out of them.

Though I'm not sure what Bryan sent me was a real phylip format file... whatever, lets do this!'''

import os
import argparse

def argument_parser(hlp = False):
    '''parse my args!'''
    
    parser = argparse.ArgumentParser(description = 'nexus2phylip.py',\
                                     argument_default = None, fromfile_prefix_chars = '@')
    parser.add_argument('-o', '--outfile', nargs = '?', type = str, default = os.getcwd() + '/AL.phylip',\
                        help = 'Path where the phylip file with results will be saved. (default: %(default)s)')
    parser.add_argument('-i', '--infile', nargs = '?', type = str, \
                        help = 'Nexus to be parsed into a nice phylip.')
    parser.add_argument('-a', '--al_run', action = 'store_true', dest = 'al_run', help = 'Run for all distances (v13_2).')
    if hlp:
        args = parser.parse_args(['-h'])
    else:
        args = parser.parse_args().__dict__
    return args

def folder_run(args):
    root = '/'.join(args['infile'].split('/')[:-1]) + '/'
    for nexus_file in [f for f in os.listdir(args['infile']) if f.endswith('.nexus')]:
        outfile = '.'.join(nexus_file.split('.')[:-1]) + '.phylip'
        print root+nexus_file, root+outfile
	nexus2phylip(root+nexus_file, root+outfile)

def nexus2phylip(infile, outfile): 
    n = infile.split('.')[0].split('/')[-1]
    with open(infile, 'r') as nexus:
        c = open(outfile, 'w')
        c.close()
        l = nexus.next()
        seq = ['']*4
        while l.split()[0] != 'Gorilla': #I feel this is going to be fun!
            l = nexus.next()
        while l[0] != ';':
            if len(l) > 1:
                if l.split()[0] == 'Gorilla':
                    seq[0] += l.split()[1]
                elif l.split()[0] == 'Homo':
                    seq[1] += l.split()[1]
                elif l.split()[0] == 'Pan':
                    seq[2] += l.split()[1]
                elif l.split()[0] == 'Pongo':
                    seq[3] += l.split()[1]
            l = nexus.next()

    with open(outfile, 'w') as out:
        out.write('\n 4     ' + str(len(seq[0])) + '\n\n')
        out.write('G-AL' + n + '^1\n')
        out.write(seq[0] + '\n')
        out.write('H-AL' + n + '^2\n')
        out.write(seq[1] + '\n')
        out.write('C-AL' + n + '^3\n')
        out.write(seq[2] + '\n')
        out.write('O-AL' + n + '^4\n')
        out.write(seq[3] + '\n\n')

def al_run():
    dist_list = ['200X200/', '200X10/', '10X200/', '10X10/']
    for dist in dist_list:
        print '##################'
        print 'directory :', dist
        infolder = '/home/igor/AL/final_finder/v13_2/nogap300/' + dist
        with open(infolder + 'sample.log', 'r') as sample:
            file_list = sample.readlines()
        folder = '/home/igor/AL/final_finder/v13_2/phylip/nogap300/' + dist
        if not os.path.isdir(folder):
            os.makedirs(folder)
        for f in file_list:
            print 'converting file: ', f[:-1]
            nexus2phylip({'infile': infolder + f[:-1], 'outfile': folder + f.split('.')[0]+'.phylip'})
    for dist in dist_list:
        print '##################'
        print 'directory :', dist
        infolder = '/home/igor/AL/final_finder/v13_2/gap300/' + dist
        with open(infolder + 'sample.log', 'r') as sample:
            file_list = sample.readlines()
        folder = '/home/igor/AL/final_finder/v13_2/phylip/gap300/' + dist
        if not os.path.isdir(folder):
            os.makedirs(folder)
        for f in file_list:
            print 'converting file: ', f[:-1]
            nexus2phylip({'infile': infolder + f[:-1], 'outfile': folder + f.split('.')[0]+'.phylip'})

if __name__ == '__main__':
    args = argument_parser()
    if args['al_run']:
	al_run()
    elif args['infile']:
        if os.path.isdir(args['infile']):
	    folder_run(args)
        else:
            nexus2phylip(args['infile'], args['outfile'])
    else:
	argument_parser(hlp=True)
