import sys
import os
import numpy
import moments
import pylab
from datetime import datetime
import Optimize_Functions
import Models_2D

#************** EDIT THESE LINES

#=========================================================================================
# Part 1A: Import data to create joint-site frequency spectrum (comment out if not using)
#=========================================================================================

snps = "DP-SM_MAF2dadiSNPinfile" #**************
dd = moments.Misc.make_data_dict(snps)
pop_ids=["DP", "SS"] #**************
proj = [33, 33] #**************
fs = moments.Spectrum.from_data_dict(dd, pop_ids=pop_ids, projections = proj, polarized = False)
fs.to_file("DP-SS_MAF2dadiSNP2sfs.fs") #**************
fs=fs.project(proj)
moments.Plotting.plot_single_2d_sfs(fs, vmin=5)
pylab.show()

#==========================================================================================================================
# Part 2: Import site frequency spectrumn from angsd OR reading in previously made file above (comment out if not using)
#==========================================================================================================================

pop_ids=["DP", "SS"] #**************
proj = [33, 33] #**************
fs = moments.Spectrum.from_file("DP_SS.2dsfs") #**************
fs=fs.project(proj)
fs=fs.fold()
moments.Plotting.plot_single_2d_sfs(fs, vmin=5)
pylab.show()

#===========================================================================
# Print some useful information
#===========================================================================

print("Projection: {}".format(proj))
print("Sample sizes: {}".format(fs.sample_sizes))
print("Sum of SFS: {}".format(numpy.around(fs.S(), 2)))
