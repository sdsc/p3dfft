#!/usr/bin/python
#
# desciption: Use this script to permute configurations.
#             Run in directory with p3dfft.
import os
import sys
import math
from subprocess import call

cwd = os.getcwd()
opts = sys.argv[1:]

if len(opts) < 2:
  print('usage: ./configure useeven stride1 single ... ARGN')
else:
  for i in range(pow(2,len(opts))):
    d = cwd + '/' + 'p3dfft' + `i`
    try:
      os.mkdir(d)
    except:
      pass
    call('cp -r ' + cwd + '/p3dfft/* ' + d, shell=True)
    os.chdir(d)
    b = list(bin(i))[2:]
    b = map(int,['0']*(len(opts)-len(b)) + b)
    c = './configure --enable-intel --enable-fftw'
    for i in range(len(opts)):
      if b[i]:
        c += ' --enable-' + opts[i]
    call(c + ' --with-fftw=$FFTWHOME FC=mpif90 CC=mpicc', shell=True)
    call('make', shell=True)

