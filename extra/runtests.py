#!/usr/bin/python
#description: Use this script to submit test jobs
from fractions import Fraction
from subprocess import call
import os
import math
PWD=os.getcwd()
NUMCORES=128
batchscripts = []
while NUMCORES <= 16384:
    path0 = PWD + '/test' + `NUMCORES` + '/16' + `NUMCORES/16`
    if not os.path.exists(path0): 
        os.makedirs(path0)
    path1 = PWD + '/test' + `NUMCORES` + '/8' + `NUMCORES/8`
    if not os.path.exists(path1):
        os.makedirs(path1)
    path2 = PWD + '/test' + `NUMCORES` + '/4' + `NUMCORES/4`
    if not os.path.exists(path2):
        os.makedirs(path2)
    paths = [path0,path1,path2]
    NMAX=int(math.floor(math.pow(NUMCORES*4*math.pow(10,9)/48,Fraction(1,3))))
    for n in range(127, NMAX+1):
        if ((n != 0) and not (n & (n-1))):
            for i in range(0,3):
                path = paths[i] + '/' + `n`
                if not os.path.exists(path):
                    os.makedirs(path)
                    os.chdir(path)
                    file=open("stdin","w")
                    file.write("%d %d %d 2 10\n" % (n, n, n))
                    file.close()
                    file=open("dims","w")
                    if i == 0:
                        file.write("16 %d\n" % (NUMCORES/16))
                    elif i == 1:
                        file.write("8 %d\n" % (NUMCORES/8))
                    else:
                        file.write("4 %d\n" % (NUMCORES/4))
                    file.close()
            os.chdir(PWD)
            filename = "batch" + `NUMCORES` + `n`
            batchscripts.append(filename)
            file=open(filename,"w")
            file.write("#!/bin/bash\n#PBS -q normal\n#PBS -l nodes=%d:ppn=16:native\n" % (NUMCORES/16))
            file.write("#PBS -l walltime=00:10:00\n#PBS -N P3DFFT\n#PBS -o out.$PBS_JOBID\n#PBS -e err.$PBS_JOBID\n")
            file.write("#PBS -A sds53\n#PBS -M m1samuel@ucsd.edu\n#PBS -V\n\n")
            file.write("cd ${SCRATCH}\n")
            file.write("cd " + path0 + '/' + `n` + '\n')
            file.write("mpirun_rsh -export -np "+`NUMCORES`+" -hostfile $PBS_NODEFILE ~/p3dfft-multivar/sample/FORTRAN/test_sine_f.x\n")
            file.write("cd " + path1 + '/' + `n` + '\n')
            file.write("mpirun_rsh -export -np "+`NUMCORES`+" -hostfile $PBS_NODEFILE ~/p3dfft-multivar/sample/FORTRAN/test_sine_f.x\n")
            file.write("cd " + path2 + '/' + `n` + '\n')
            file.write("mpirun_rsh -export -np "+`NUMCORES`+" -hostfile $PBS_NODEFILE ~/p3dfft-multivar/sample/FORTRAN/test_sine_f.x\n")
            file.close()
    NUMCORES *= 2
os.chdir(PWD)
for i in range(0,len(batchscripts)):
    call("qsub " + batchscripts[i],shell = True)
