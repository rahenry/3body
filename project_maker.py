#!/usr/bin/env python

#input generator
#takes parameters from file argument and creates a bunch of input files!
import subprocess
import os
import re
import time
import sys
import math

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('name', nargs='?', default=None)
parser.add_argument('--kappa_min', nargs='+', default=None)
parser.add_argument('--kappa_max', nargs='+', default='1')
parser.add_argument('--kappa_log', nargs='?', default=None, const='1')
parser.add_argument('--kappa_number', nargs='+', default='1')
parser.add_argument('--kappa_reciprocal', nargs='?', default=None, const='1')
parser.add_argument('--etas', nargs='+', default='1')
parser.add_argument('--kappas', nargs='+', default='1')
parser.add_argument('--Ls', nargs='+', default='1')
parser.add_argument('--sf_channels', nargs='+', default=['3', '2', '0'])

args = parser.parse_args(sys.argv[1:])

result = '#'
for a in sys.argv:
    result += a + ' '
result += '\n'

def kfun(x):
    if (args.kappa_log):
        return(math.log(x))
    else:
        return x

def kfun_inv(x):
    if (args.kappa_log):
        return(math.exp(x))
    else:
        return x


if (len(args.kappa_max) != len(args.kappa_number)):
    print "WHAT IS WRONG WITH YOU?"
    exit()
if (args.kappa_min):
    if (len(args.kappa_max) != len(args.kappa_min) and not args.kappa_reciprocal):
        print "WHAT IS WRONG WITH YOU?"
        exit()

def get_min(i):
    if (not args.kappa_min):
        return 1.0
    if (args.kappa_reciprocal and args.kappa_log and args.kappa_min[i] == 'auto'):
        return(-kfun(float(args.kappa_max[i])))
    else:
        return(kfun(float(args.kappa_min[i])))

kappas = []
if (args.kappas != '1'):
    kappas = args.kappas
else: 
    for i in range(0, len(args.kappa_max)):
        kappas_new = []
        kmax = kfun(float(args.kappa_max[i]))
        kmin = get_min(i)
        nk = int(args.kappa_number[i])

        if (nk == 1):
            kappas_new = [kmax]
        else:
            kappas_new = []
            kgrad = (kmax - kmin) / (nk - 1)
            for i in range(0, nk):
                kappas_new.append(kmin + i * kgrad)

        kappas += map(kfun_inv, kappas_new)

result += 'kappas '
for k in kappas:
    result += str(k) + ' '
result += '\n'

result += 'etas '
for e in args.etas:
    result += str(e) + ' '
result += '\n'

result += 'Ls '
for l in args.Ls:
    result += str(l) + ' '
result += '\n'

result += 'sf_channels '
for s in args.sf_channels:
    result += str(s) + ' '
result += '\n'

if (args.name):
    #print(result)
    f = open(args.name, 'w+')
    f.write(result)
else:
    print(result)
    
    

