#!/usr/bin/sh

for i in `seq 50 10 100`
do
  for j in `seq 10 10 100`
  do
    for k in `seq 80 5 95`
    do
      Rscript robustness_evaluation.R -e 2 -n $j -c $i -f $k -s 2019 -d "/home/mribeirodantas/dev/miic_r_package/tmp/figs/loops/edge_filtering_only_fulldataset"
    done
  done
done
