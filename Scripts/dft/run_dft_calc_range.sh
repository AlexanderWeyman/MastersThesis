#!/bin/bash

# usage: ./run_dft_calc_range.sh Nmin Nmax Nstep delta_s eta l_Bjerrum q_min q_max q_step

Nmin=$1
Nmax=$2
Nstep=$3
deltaS=$4
eta=$5
l_b=$6
qmin=$7
qmax=$8
qstep=$9

for N in $(seq $Nmin $Nstep $Nmax)
do
    ./calcPartStructurefactors.py $N $deltaS $eta $l_b $qmin $qmax $qstep
done
