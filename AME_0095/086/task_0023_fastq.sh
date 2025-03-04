#! /bin/bash
#PBS -l select=1:ncpus=64:mpiprocs=1:ompthreads=12:jobtype=largemem
#PBS -l walltime=4:00:00

echo "start: " `date "+%Y-%m-%d %H:%M:%S"`

# PATHの設定（既存のPATHを保持する）
export PATH=$PATH:/fck@ccfep.ims.ac.jp:/home/users/fck/sratoolkit.3.1.1-centos_linux64/bin

cd /home/users/fch/task_0023

# SRAファイルのリストを処理
while read sra; do
  echo "Processing /home/users/fch/task_0023/${sra}.sra"

  # fasterq-dumpを実行
  fasterq-dump "/home/users/fch/task_0023/${sra}.sra" -O /home/users/fch/task_0023/fastq -e 4 --verbose
done < /home/users/fch/task_0023/SRA_Acc_list.txt

echo "end: " `date "+%Y-%m-%d %H:%M:%S"`



