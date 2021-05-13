# This script finds the number of SNPs for all possible projections for two populations with equal sample size.

import sys
import os
import numpy
import moments
import pylab
import pandas as pd
from datetime import datetime

mylist=[]

for x in range(65,1,-1): # range starts with number of total alleles-1
	fs = moments.Spectrum.from_file("DP_SS.2dsfs")
	proj = [x, x]
	fs=fs.project(proj)
	mylist.append(format(numpy.around(fs.S(), 2)))

# create a pandas dataframe using pd.DataFrame function with a dictionary with two lists as input.
df = pd.DataFrame({'proj':range(65,1,-1), 'SNPs':mylist})

print(df)
df.to_csv(r'find-projection.csv', index = False)