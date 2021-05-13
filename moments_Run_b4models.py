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

snps = "Pop1-Pop2_MAF2dadiSNPinfile" #**************
dd = moments.Misc.make_data_dict(snps)
pop_ids=["Pop1", "Pop2"] #**************
proj = [33, 33] #**************
fs = moments.Spectrum.from_data_dict(dd, pop_ids=pop_ids, projections = proj, polarized = False)
fs.to_file("Pop1_Pop2.2dsfs") #**************
fs=fs.project(proj)
moments.Plotting.plot_single_2d_sfs(fs, vmin=5)
pylab.show()

#==========================================================================================================================
# Part 2: Import site frequency spectrumn from angsd OR reading in previously made file above (comment out if not using)
#==========================================================================================================================

pop_ids=["Pop1", "Pop2"] #**************
proj = [33, 33] #**************
fs = moments.Spectrum.from_file("Pop1_Pop2.2dsfs") #**************
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
