#!/usr/bin/python

# This script will generate a permutative batch job script on Mira.
#
# More specifically, the generated batch script will test for correctness
# of all tests (under sample/FORTRAN and sample/C) for each configuration
# of p3dfft available (see below). Because we only care about correctness,
# 1 node with 16 cores will be used with four dimensions: 4x4, 1x16, 16x1, 1x1.
#
# The p3dfft directories need to be in the current working directory from which
# this script is executed, named 'p3dfft0' to 'p3dfftX' (use the configure.py
# script for that).
#
# The jobs/ directory needs to exist in the current working directory. This
# is where all batch job files are written to.


# Make another script to run the resultant script:
# #!/bin/bash
# qsub -A P3DFFT -t 30 -n 1 --mode script mira_permutative.sh

import os
import re

NUMCORES = 16

# Factorisation helper function
def factors(n):
    return set(reduce(list.__add__,
        ([i, n//i] for i in range(1, int(n**0.5) + 1) if n % i == 0)))

    # Open batch job file to be written to.
batchf = open('jobs/mira_permutative.sh', 'w')

# Write SBATCH header commands.
batchf.write('#!/bin/bash\n')
batchf.write('\n')

basedir = os.getcwd()

# Get all p3dfft config directories
p3dfft_dirs = next(os.walk('.'))[1]
pattern = re.compile('p3dfft\d+$')
p3dfft_dirs = sorted(filter(pattern.match, p3dfft_dirs))

# Get all test names using first directory
f_dir = os.path.join(p3dfft_dirs[0], 'sample/FORTRAN')
c_dir = os.path.join(p3dfft_dirs[0], 'sample/C')

pattern = re.compile('test_\S+_[cf].x')
f_test_files = filter(pattern.match, next(os.walk(f_dir))[2])
c_test_files = filter(pattern.match, next(os.walk(c_dir))[2])

# Get full paths to tests in all dirs
all_tests = []

for d in p3dfft_dirs:
    f_dir = os.path.join(d, 'sample/FORTRAN')
    c_dir = os.path.join(d, 'sample/C')
    for test in f_test_files:
        all_tests.append(os.path.join(f_dir, test))
    for test in c_test_files:
        all_tests.append(os.path.join(c_dir, test))

# Calculate dims
all_dims = []
facs = sorted(factors(NUMCORES))
if (len(facs) % 2 == 0):
    # take the two factors in the middle
    all_dims.append("'" + str(facs[len(facs)/2-1]) + " " + str(facs[len(facs)/2]) + "'")
else:
    # perfect square, take the factor in the middle
    all_dims.append("'" + str(facs[len(facs)/2]) + " " + str(facs[len(facs)/2]) + "'")
all_dims.append("'" + str(facs[len(facs)-1]) + " " + str(facs[0]) + "'")
all_dims.append("'" + str(facs[0]) + " " + str(facs[len(facs)-1]) + "'")
#all_dims = ["'4 4'", "'16 1'", "'1 16'", "'1 1'"]

# Run all tests for all dims
for test in all_tests:
    if "cheby" in test:
        batchf.write("echo '128 128 129 2 1' > stdin\n")
    elif "many" in test:
        batchf.write("echo '128 128 128 2 5 1' > stdin\n")
    else:
        batchf.write("echo '128 128 128 2 1' > stdin\n")
    for dims in all_dims:
        # write dims
        batchf.write("echo " + dims + " > dims\n")
        # run test
        batchf.write("runjob -p " + str(NUMCORES) + " --np " + str(NUMCORES) + " " + basedir + "/" + test + "\n")
    # 1x1 dims test
    batchf.write("echo '1 1' > dims\n")
    batchf.write("runjob -p 1 --np 1 " + basedir + "/" + test + "\n")

# Truncate previous content if any existed.
#batchf.truncate()

# Close the file. Done.
batchf.close()
