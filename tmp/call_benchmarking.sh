#!/usr/bin/sh

for j in `seq 50 10 100`
do
  for k in `seq 80 5 95`
  do
    Rscript benchmark.R -e 0 -k 100 -c $j -f $k -n 1000 -s 2019 -d "barley" -o "/home/mribeirodantas/dev/miic_r_package/tmp/data/output/barley_1000_samples"
  done
done

