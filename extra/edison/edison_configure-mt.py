#!/usr/bin/python
#
# desciption: Use this script to permute configurations.
#             Run in directory with p3dfft.
# Options: stride1, useeven, single

import os
import sys
import math
from subprocess import call

cwd = os.getcwd()
opts = sys.argv[1:]

if len(opts) < 2:
  print('usage: ./configure ARG0 ARG1... ARGN')
else:
  for i in range(pow(2,len(opts))):
    d = cwd + '/' + 'p3dfft' + `i` + '-mt'
    try:
      os.mkdir(d)
    except:
      pass
    call('cp -r ' + cwd + '/p3dfft-mt/* ' + d, shell=True)
    os.chdir(d)
    b = list(bin(i))[2:]
    b = map(int,['0']*(len(opts)-len(b)) + b)
    c = './configure --enable-intel --enable-fftw --with-fftw=$FFTW_DIR/.. FC=ftn CC=cc FCFLAGS=-openmp CFLAGS=-openmp'
    for i in range(len(opts)):
      if b[i]:
        c += ' --enable-' + opts[i]
    call(c, shell=True)
    call('make', shell=True)
