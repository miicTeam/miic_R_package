#!/usr/bin/sh

for i in `seq 10 10 100`
do
  for j in `seq 50 10 100`
  do
    for k in `seq 80 5 95`
    do
      Rscript benchmark.R -e 0 -k $i -c $j -f $k -n 1000 -s 2019 -d "insurance" -o "/home/mribeirodantas/dev/miic_r_package/tmp/figs/benchmark/edge_filtering"
    done
  done
done
