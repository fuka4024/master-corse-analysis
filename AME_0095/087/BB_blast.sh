#! /bin/bash
#PBS -l select=1:ncpus=64:mpiprocs=1:ompthreads=12:jobtype=largemem
#PBS -l walltime=8:00:00



echo "start: " `date "+%Y-%m-%d %H:%M:%S"`

source /apl/bio/etc/bio.sh
module load blast+



makeblastdb -in /home/users/fch/task_0017/GCF_014905235.1_Bmori_2016v1.0_rna.fna -out Bmori -dbtype nucl -parse_seqids
blastn -db Bmori -query /home/fuuka40/rna.fna -outfmt 6  -num_threads=8 > /home/fuuka40/task_0024/BB_blast_result.txt


echo "end: " `date "+%Y-%m-%d %H:%M:%S"`