#! /bin/bash
#PBS -l select=1:ncpus=64:mpiprocs=1:ompthreads=12:jobtype=largemem
#PBS -l walltime=4:00:00



echo "start: " `date "+%Y-%m-%d %H:%M:%S"`

source /apl/bio/etc/bio.sh
module load blast+



makeblastdb -in /home/users/fch/task_GCF_000001215/GCF_000001215.4_Release_6_plus_ISO1_MT_protein.faa -out Dm -dbtype prot -parse_seqids
blastx -db Dm -query /home/users/fch/task_0025/BmToll.fna -outfmt 6  -num_threads=8 > /home/users/fch/task_0025/BmToll_Dm_blast_result.txt


echo "end: " `date "+%Y-%m-%d %H:%M:%S"`
