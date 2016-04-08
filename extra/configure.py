#!/usr/bin/python2

import getopt
import sys
import os
import math
from subprocess import call

platforms = ["comet", "gordon", "edison", "cori", "stampede"]   # todo: add mira
compilers = ["intel", "gnu", "pgi", "cray", "ibm"]
options   = ["useeven", "stride1", "single"]
configs   = { ("comet", False, False)   : './configure --enable-fftw --with-fftw=$FFTWHOME FC=mpif90 CC=mpicc',
              ("comet", True, False)    : './configure --enable-fftw --with-fftw=$FFTWHOME FC=mpif90 CC=mpicc --enable-openmp',
              ("comet", False, True)    : './configure --enable-fftw --with-fftw=$FFTWHOME FC=mpif90 CC=mpicc --enable-stride1 CFLAGS="-O3"',
              ("comet", True, True)     : './configure --enable-fftw --with-fftw=$FFTWHOME FC=mpif90 CC=mpicc --enable-stride1 CFLAGS="-O3" --enable-openmp',
              ("gordon", False, False)  : './configure --enable-fftw --with-fftw=$FFTWHOME FC=mpif90 CC=mpicc',
              ("gordon", True, False)   : './configure --enable-fftw --with-fftw=$FFTWHOME FC=mpif90 CC=mpicc --enable-openmp',
              ("gordon", False, True)   : './configure --enable-fftw --with-fftw=$FFTWHOME FC=mpif90 CC=mpicc --enable-stride1 CFLAGS="-O3"',
              ("gordon", True, True)    : './configure --enable-fftw --with-fftw=$FFTWHOME FC=mpif90 CC=mpicc --enable-stride1 CFLAGS="-O3" --enable-openmp',
              ("edison", False, False)  : './configure --enable-fftw --with-fftw=$FFTW_DIR/.. FC=ftn CC=cc',
              ("edison", True, False)   : './configure --enable-fftw --with-fftw=$FFTW_DIR/.. FC=ftn CC=cc --enable-openmp',
              ("edison", False, True)   : './configure --enable-fftw --with-fftw=$FFTW_DIR/.. FC=ftn CC=cc --enable-stride1 CFLAGS="-O3"',
              ("edison", True, True)    : './configure --enable-fftw --with-fftw=$FFTW_DIR/.. FC=ftn CC=cc --enable-stride1 CFLAGS="-O3 --enable-openmp"',
              ("cori", False, False)    : './configure --enable-fftw --with-fftw=$FFTW_DIR/.. FC=ftn CC=cc',
              ("cori", True, False)     : './configure --enable-fftw --with-fftw=$FFTW_DIR/.. FC=ftn CC=cc --enable-openmp',
              ("cori", False, True)     : './configure --enable-fftw --with-fftw=$FFTW_DIR/.. FC=ftn CC=cc --enable-stride1 CFLAGS="-O3"',
              ("cori", True, True)      : './configure --enable-fftw --with-fftw=$FFTW_DIR/.. FC=ftn CC=cc --enable-stride1 CFLAGS="-O3" --enable-openmp',
              ("stampede", False, False): './configure --enable-fftw --with-fftw=$TACC_FFTW3_DIR FC=mpif90 CC=mpicc',
              ("stampede", True, False) : './configure --enable-fftw --with-fftw=$TACC_FFTW3_DIR FC=mpif90 CC=mpicc --enable-openmp',
              ("stampede", False, True) : './configure --enable-fftw --with-fftw=$TACC_FFTW3_DIR FC=mpif90 CC=mpicc --enable-stride1 CFLAGS="-O3"',
              ("stampede", True, True)  : './configure --enable-fftw --with-fftw=$TACC_FFTW3_DIR FC=mpif90 CC=mpicc --enable-stride1 CFLAGS="-O3" --enable-openmp',
            }
sourcedir = "p3dfft"
destdir   = { (False, False) : "p3dfft", (True, False) : "p3dfft-mt",
              (False, True)  : "p3dfft-p", (True, True) : "p3dfft-mt-p"
            }

def usage_exit(msg):
    print msg
    print "USAGE: ./configure.py -s comet|gordon|edison|cori|stampede [-c intel|gnu|pgi|cray|ibm] [-m] [-p] [-f extra flags]"
    print "-s  specifies which platform"
    print "-c  to specify non-default compiler"
    print "-m  to build -mt branch"
    print "-p  to build performance test"
    print "-f  extra configure flags"
    sys.exit(1)

def main():
    platform = None
    comp = None
    extra = None
    mt = False
    perf = False

    # parse command line options
    try:
        opts = getopt.getopt(sys.argv[1:], 's:c:mpf:')
    except getopt.GetoptError as err:
        usage_exit(str(err))
    for o, a in opts[0]:
        if o == '-s':
            platform = a
        elif o == '-c':
            comp = a
        elif o == '-m':
            mt = True
        elif o == '-p':
            perf = True
        elif o == '-f':
            extra = a
        else:
            assert False, "unhandled option"
    if platform == None:
        usage_exit("no platform specified")
    elif platform not in platforms:
        usage_exit("invalid platform specified")

    # make configline according to compiler
    configline = configs[(platform, mt, perf)]
    if comp == None:
        configline += " --enable-intel"
    else:
        if comp not in compilers:
            usage_exit("invalid compiler specified")
        configline += " --enable-" + comp
    if extra != None:
        configline += " " + extra
    if comp == "gnu":
        configline += " CFLAGS=-lm"

    # ensure that the source dir exists
    source = sourcedir
    dest = destdir[(mt, perf)]
    cwd = os.getcwd()
    if not os.path.isdir(cwd + '/' + source):
        usage_exit(source + " dir does not exist. Check it.")

    # start build
    print configline
    print "Source Directory: " + source
    print "Destination Directory: " + dest
    print "********** Starting build... **********"

    if perf:
        d = cwd + '/' + dest
        try:
            os.mkdir(d)
        except:
            pass
        call('cp -r ' + cwd + '/' + source + '/* ' + d, shell=True)
        os.chdir(d)
        c = configline
        call(c, shell=True)
        call('make', shell=True)
    else:
        for i in range(pow(2,len(options))):
            d = cwd + '/' + dest + str(i)
            try:
                os.mkdir(d)
            except:
                pass
            call('cp -r ' + cwd + '/' + source + '/* ' + d, shell=True)
            os.chdir(d)
            b = list(bin(i))[2:]
            b = map(int,['0']*(len(options)-len(b)) + b)
            c = configline
            for i in range(len(options)):
                if b[i]:
                    c += ' --enable-' + options[i]
            print "Configuring " + d + " with "
            print "\t" + c
            c += " > config_output"
            ret = call(c, shell=True)
            if ret != 0:
                usage_exit("CONFIG FAILED! CHECK config_output for log")
            print "Configured " + d
            print "Making " + d
            ret = call('make > make_output', shell=True)
            if ret != 0:
                usage_exit("MAKE FAILED! CHECK make_output for log")
            print "Built " + d

    print "********** Done. **********"

if __name__ == '__main__':
    main()
