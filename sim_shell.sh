#!/bin/bash

param_file=$1
echo $param_file

pass_parameters() {
  param="$1"
  Rscript sim.R "$param"
}

export -f pass_parameters

cat $param_file|parallel --will-cite -j 32 --retry-failed --use-cpus-instead-of-cores pass_parameters

exit
