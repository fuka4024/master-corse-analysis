#! /bin/bash
#PBS -l select=1:ncpus=64:mpiprocs=1:ompthreads=64:jobtype=largemem
#PBS -l walltime=4:00:00

echo "start: " `date "+%Y-%m-%d %H:%M:%S"`

source /apl/bio/etc/bio.sh
module load salmon
cd task_0023/output/

salmon quantmerge --colum numreads -o /home/users/fch/task_0023/output/FB_EP.txt --quants {1..18}

echo "end: " `date "+%Y-%m-%d %H:%M:%S"`



