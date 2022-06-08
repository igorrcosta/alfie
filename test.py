#!/usr/bin/env python
# -*- coding: utf-8 -*-
#test.py
#11/2015

from subprocess import Popen
from shlex import split as ssplit

def test():
    command = 'python alfie.py -r test/homo_y.fasta -i test/pan_y.fasta -g test/y.gtf -o tmp/ -v -l log'
    command = ssplit(command)
    a = Popen(command)
    a.wait()

test()
