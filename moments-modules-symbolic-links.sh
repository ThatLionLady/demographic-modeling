#!/bin/sh

MPDIR=$1 # path to moments_pipeline directory
WKDIR=$2 # path to working directory

ln -s ${MPDIR}/Two_Population_Pipeline/Models_2D.py ${WKDIR}/Models_2D.py
ln -s ${MPDIR}/Two_Population_Pipeline/Summarize_Outputs.py ${WKDIR}/Summarize_Outputs.py
ln -s ${MPDIR}/Optimize_Functions.py Optimize_${WKDIR}/Functions.py
ln -s ${MPDIR}/Goodness_of_Fit/Optimize_Functions_GOF.py ${WKDIR}/Optimize_Functions_GOF.py