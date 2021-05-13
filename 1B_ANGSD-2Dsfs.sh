#!/bin/sh

realSFS=$1 # Path to realSFS in angsd program directory
dir=$2 # Path to working directory

Pop1=$3 # Path to Pop1_HiRise_rh_HiC_nosites.saf.idx
Pop2=$4 # Path to Pop2_HiRise_rh_HiC_nosites.saf.idx

${realSFS} ${Pop1} ${Pop2} -nSites 10000000 -P 8 > ${dir}2Dsfs.Pop1_Pop2.sfs 2> ${dir}2Dsfs.Pop1_Pop2.log