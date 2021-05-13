#!/bin/sh

DIR=$1 # path to moments_pipeline

ln -s ${DIR}/Two_Population_Pipeline/Models_2D.py Models_2D.py
ln -s ${DIR}/Two_Population_Pipeline/Summarize_Outputs.py Summarize_Outputs.py
ln -s ${DIR}/Optimize_Functions.py Optimize_Functions.py
ln -s ${DIR}/Goodness_of_Fit/Optimize_Functions_GOF.py Optimize_Functions_GOF.py