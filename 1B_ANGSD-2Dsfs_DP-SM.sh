#!/bin/sh

realSFS=/home/ubuntu/USS/caitlin/Programs/angsd/misc/realSFS
dir=/home/ubuntu/USS/caitlin/PPM-demography/

SM_PPM=/home/ubuntu/USS/conservation-genetics/PPM/Angsd_HiC/SM_PPM_HiRise_rh_HiC_nosites.saf.idx
DP_PPM=/home/ubuntu/USS/conservation-genetics/PPM/Angsd_HiC/DP_PPM_HiRise_rh_HiC_nosites.saf.idx

${realSFS} ${DP_PPM} ${SM_PPM} -nSites 10000000 -P 8 > ${dir}2Dsfs.DP_SM.sfs 2> ${dir}2Dsfs.DP_SM.log