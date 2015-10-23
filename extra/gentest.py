#!/usr/bin/python

NUMCORES=64

file=open("tstfnctnlty","w")

def wrt1(fh,i):
    fh.write("mpirun_rsh -export -np " + `NUMCORES` + " -hostfile $PBS_NODEFILE ~/p3dfft" + `i` + "/sample/FORTRAN/test_rand_many_f.x\n")
    fh.write("mpirun_rsh -export -np " + `NUMCORES` + " -hostfile $PBS_NODEFILE ~/p3dfft" + `i` + "/sample/FORTRAN/test_sine_many_f.x\n")
    fh.write("mpirun_rsh -export -np " + `NUMCORES` + " -hostfile $PBS_NODEFILE ~/p3dfft" + `i` + "/sample/FORTRAN/test_sine_inplace_many_f.x\n")

def wrt2(fh,i):
    fh.write("mpirun_rsh -export -np " + `NUMCORES` + " -hostfile $PBS_NODEFILE ~/p3dfft" + `i` + "/sample/FORTRAN/test_inverse_f.x\n")
    fh.write("mpirun_rsh -export -np " + `NUMCORES` + " -hostfile $PBS_NODEFILE ~/p3dfft" + `i` + "/sample/FORTRAN/test_noop_f.x\n")
    fh.write("mpirun_rsh -export -np " + `NUMCORES` + " -hostfile $PBS_NODEFILE ~/p3dfft" + `i` + "/sample/FORTRAN/test_sine_f.x\n")
    fh.write("mpirun_rsh -export -np " + `NUMCORES` + " -hostfile $PBS_NODEFILE ~/p3dfft" + `i` + "/sample/FORTRAN/test_rand_f.x\n")
    fh.write("mpirun_rsh -export -np " + `NUMCORES` + " -hostfile $PBS_NODEFILE ~/p3dfft" + `i` + "/sample/FORTRAN/test_sine_inplace_f.x\n")
    fh.write("mpirun_rsh -export -np " + `NUMCORES` + " -hostfile $PBS_NODEFILE ~/p3dfft" + `i` + "/sample/FORTRAN/test_spec_f.x\n")
    fh.write("mpirun_rsh -export -np " + `NUMCORES` + " -hostfile $PBS_NODEFILE ~/p3dfft" + `i` + "/sample/C/test_inverse_c.x\n")
    fh.write("mpirun_rsh -export -np " + `NUMCORES` + " -hostfile $PBS_NODEFILE ~/p3dfft" + `i` + "/sample/C/test_sine_c.x\n")
    fh.write("mpirun_rsh -export -np " + `NUMCORES` + " -hostfile $PBS_NODEFILE ~/p3dfft" + `i` + "/sample/C/test_rand_c.x\n")
    fh.write("mpirun_rsh -export -np " + `NUMCORES` + " -hostfile $PBS_NODEFILE ~/p3dfft" + `i` + "/sample/C/test_sine_inplace_c.x\n")
    fh.write("mpirun_rsh -export -np " + `NUMCORES` + " -hostfile $PBS_NODEFILE ~/p3dfft" + `i` + "/sample/C/test_spec_c.x\n")

file.write("#!/bin/bash\n#PBS -q normal\n#PBS -l nodes=%d:ppn=16:native\n" % (NUMCORES/16))
file.write("#PBS -l walltime=01:00:00\n#PBS -N P3DFFT\n#PBS -o out.$PBS_JOBID\n#PBS -e err.$PBS_JOBID\n")
file.write("#PBS -A sds153\n#PBS -M m1samuel@ucsd.edu\n#PBS -V\n\n")
file.write("cd ${SCRATCH}\n")

file.write("cd ${SCRATCH}/many/88\n")
for i in range(16):
    wrt1(file,i)


file.write("cd ${SCRATCH}/many/164\n")
for i in range(16):
    wrt1(file,i)

file.write("cd ${SCRATCH}/many/641\n")
for i in range(16):
    wrt1(file,i)

file.write("cd ${SCRATCH}/many/11\n")
for i in range(16):
    file.write("mpirun_rsh -export -np 1 -hostfile $PBS_NODEFILE ~/p3dfft" + `i`+ "/sample/FORTRAN/test_rand_many_f.x\n")
    file.write("mpirun_rsh -export -np 1 -hostfile $PBS_NODEFILE ~/p3dfft" + `i`+ "/sample/FORTRAN/test_sine_many_f.x\n")
    file.write("mpirun_rsh -export -np 1 -hostfile $PBS_NODEFILE ~/p3dfft" + `i`+ "/sample/FORTRAN/test_sine_inplace_many_f.x\n")

file.write("cd ${SCRATCH}/cheby/88\n")
for i in range(16):
    file.write("mpirun_rsh -export -np " + `NUMCORES` + " -hostfile $PBS_NODEFILE ~/p3dfft" + `i` + "/sample/FORTRAN/test_cheby_f.x\n")

file.write("cd ${SCRATCH}/cheby/164\n")
for i in range(16):
    file.write("mpirun_rsh -export -np " + `NUMCORES` + " -hostfile $PBS_NODEFILE ~/p3dfft" + `i` + "/sample/FORTRAN/test_cheby_f.x\n")

file.write("cd ${SCRATCH}/cheby/641\n")
for i in range(16):
    file.write("mpirun_rsh -export -np " + `NUMCORES` + " -hostfile $PBS_NODEFILE ~/p3dfft" + `i` + "/sample/FORTRAN/test_cheby_f.x\n")

file.write("cd ${SCRATCH}/cheby/11\n")
for i in range(16):
    file.write("mpirun_rsh -export -np 1 -hostfile $PBS_NODEFILE ~/p3dfft" + `i` + "/sample/FORTRAN/test_cheby_f.x\n")

file.write("cd ${SCRATCH}/other/88\n")
for i in range(16):
    wrt2(file,i)

file.write("cd ${SCRATCH}/other/164\n")
for i in range(16):
    wrt2(file,i)

file.write("cd ${SCRATCH}/other/641\n")
for i in range(16):
    wrt2(file,i)

file.write("cd ${SCRATCH}/other/11\n")
for i in range(16):
    file.write("mpirun_rsh -export -np 1 -hostfile $PBS_NODEFILE ~/p3dfft" + `i` + "/sample/FORTRAN/test_inverse_f.x\n")
    file.write("mpirun_rsh -export -np 1 -hostfile $PBS_NODEFILE ~/p3dfft" + `i` + "/sample/FORTRAN/test_noop_f.x\n")
    file.write("mpirun_rsh -export -np 1 -hostfile $PBS_NODEFILE ~/p3dfft" + `i` + "/sample/FORTRAN/test_sine_f.x\n")
    file.write("mpirun_rsh -export -np 1 -hostfile $PBS_NODEFILE ~/p3dfft" + `i` + "/sample/FORTRAN/test_rand_f.x\n")
    file.write("mpirun_rsh -export -np 1 -hostfile $PBS_NODEFILE ~/p3dfft" + `i` + "/sample/FORTRAN/test_sine_inplace_f.x\n")
    file.write("mpirun_rsh -export -np 1 -hostfile $PBS_NODEFILE ~/p3dfft" + `i` + "/sample/FORTRAN/test_spec_f.x\n")
    file.write("mpirun_rsh -export -np 1 -hostfile $PBS_NODEFILE ~/p3dfft" + `i` + "/sample/C/test_inverse_c.x\n")
    file.write("mpirun_rsh -export -np 1 -hostfile $PBS_NODEFILE ~/p3dfft" + `i` + "/sample/C/test_sine_c.x\n")
    file.write("mpirun_rsh -export -np 1 -hostfile $PBS_NODEFILE ~/p3dfft" + `i` + "/sample/C/test_rand_c.x\n")
    file.write("mpirun_rsh -export -np 1 -hostfile $PBS_NODEFILE ~/p3dfft" + `i` + "/sample/C/test_sine_inplace_c.x\n")
    file.write("mpirun_rsh -export -np 1 -hostfile $PBS_NODEFILE ~/p3dfft" + `i` + "/sample/C/test_spec_c.x\n") 
