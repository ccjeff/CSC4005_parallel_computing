#!/bin/bash
#PBS -l nodes=1:ppn=5,mem=1g,walltime=00:30:00
#PBS -q batch
#PBS -m abe
#PBS -V

for ((pro=20; pro<=20; pro++))
#do
#for ((arr_size=20; arr_size<=1000; arr_size=arr_size+100))
  do
    echo "processor: $pro arr_size: 20000"
    mpiexec -n $pro /code/117010008/a.out 20000 1
    echo -e "\n"
# done
done