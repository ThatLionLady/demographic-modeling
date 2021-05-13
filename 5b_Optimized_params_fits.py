import sys
import os
import numpy
import moments
import pylab
from datetime import datetime
import Optimize_Functions_GOF
import Models_2D
from Models_2D import no_mig
from Models_2D import sym_mig
from Models_2D import asym_mig
from Models_2D import anc_sym_mig
from Models_2D import anc_asym_mig
from Models_2D import sec_contact_sym_mig
from Models_2D import sec_contact_asym_mig
from Models_2D import no_mig_size
from Models_2D import sym_mig_size
from Models_2D import asym_mig_size
from Models_2D import anc_sym_mig_size
from Models_2D import anc_asym_mig_size
from Models_2D import sec_contact_sym_mig_size
from Models_2D import sec_contact_asym_mig_size

# Import SFS

fs = moments.Spectrum.from_file("Pop1_Pop2.2dsfs") #**************
print("Sum of SFS: {}".format(numpy.around(fs.S(), 2)))

# Find projection will use the largest number of SNPs

fs = moments.Spectrum.from_file("Pop1_Pop2.2dsfs") #**************
proj = [10, 10] #**************
fs=fs.project(proj)
print("Sum of SFS: {}".format(numpy.around(fs.S(), 2)))

# Project to largest number of SNPs

pop_ids=["Pop1", "Pop2"] #**************
proj = [33, 33] #**************
fs=fs.project(proj)
fs=fs.fold()
fs_folded = True

#these are the optimized parameters in the best-fitting model
#the things that need to be changed for each model are commented after the line 
emp_params=[0.2417,20.6599,0.0917,0.1644,4.8612,0.2341,0.0406] #copy from Pop1Pop2_compareMomentsModels.R output
#generate sfs from optimized params
scaled_fs = Optimize_Functions_GOF.Optimize_Empirical(fs, "Pop1Pop2_bestmod","sec_contact_sym_mig_size",sec_contact_sym_mig_size, emp_params, fs_folded=fs_folded) #(fs, "outfile", "model_name", model_name, emp_params, fs_folded=True)
model = Models_2D.sec_contact_sym_mig_size(emp_params, fs.sample_sizes) 
#plot the observed vs. fit SFS with residuals between observed and fit (need to adjust vmin and resid-range for your data)
moments.Plotting.plot_2d_comp_multinom(model,fs, vmin=0.4,resid_range=5) #can fiddle with vmin & resid_range

# All models

emp_params=[0.2417,20.6599,0.0917,0.1644,4.8612,0.2341,0.0406]
scaled_fs = Optimize_Functions_GOF.Optimize_Empirical(fs, "Pop1Pop2_comp","sec_contact_sym_mig_size",sec_contact_sym_mig_size, emp_params, fs_folded=fs_folded) #(fs, "outfile", "model_name", model_name, emp_params, fs_folded=True)
model = Models_2D.sec_contact_sym_mig_size(emp_params, fs.sample_sizes) 
moments.Plotting.plot_2d_comp_multinom(model,fs, vmin=0.4,resid_range=5)

emp_params=[0.1713,4.9709,0.0913,0.229,3.4997,0.1978,0.0733]
scaled_fs = Optimize_Functions_GOF.Optimize_Empirical(fs, "Pop1Pop2_comp","sym_mig_size",sym_mig_size, emp_params, fs_folded=fs_folded)
model = Models_2D.sym_mig_size(emp_params, fs.sample_sizes)
moments.Plotting.plot_2d_comp_multinom(model,fs, vmin=0.4,resid_range=5)

emp_params=[0.0869,0.2795,3.1256,0.0967]
scaled_fs = Optimize_Functions_GOF.Optimize_Empirical(fs, "Pop1Pop2_comp","sym_mig",sym_mig, emp_params, fs_folded=fs_folded)
model = Models_2D.sym_mig(emp_params, fs.sample_sizes)
moments.Plotting.plot_2d_comp_multinom(model,fs, vmin=0.4,resid_range=5)

emp_params=[7.1512,1.5002,0.0878,0.3801,4.5254,2.0757,5.4081,1.1072]
scaled_fs = Optimize_Functions_GOF.Optimize_Empirical(fs, "Pop1Pop2_comp", "asym_mig_size", asym_mig_size, emp_params, fs_folded=True)
model = Models_2D.asym_mig_size(emp_params, fs.sample_sizes) 
moments.Plotting.plot_2d_comp_multinom(model,fs, vmin=0.4,resid_range=5)

emp_params=[0.0434,0.1443,8.3089,5.886,0.1297,0.3704]
scaled_fs = Optimize_Functions_GOF.Optimize_Empirical(fs, "Pop1Pop2_comp", "sec_contact_asym_mig", sec_contact_asym_mig, emp_params, fs_folded=True)
model = Models_2D.sec_contact_asym_mig(emp_params, fs.sample_sizes) 
moments.Plotting.plot_2d_comp_multinom(model,fs, vmin=0.4,resid_range=5)

emp_params=[0.4681,1.391,0.7351,12.0598,3.1693]
scaled_fs = Optimize_Functions_GOF.Optimize_Empirical(fs, "Pop1Pop2_comp", "sec_contact_sym_mig", sec_contact_sym_mig, emp_params, fs_folded=True)
model = Models_2D.sec_contact_sym_mig(emp_params, fs.sample_sizes) 
moments.Plotting.plot_2d_comp_multinom(model,fs, vmin=0.4,resid_range=5)

emp_params=[0.1564,0.3588,1.6215,2.8602,0.2312]
scaled_fs = Optimize_Functions_GOF.Optimize_Empirical(fs, "Pop1Pop2_comp", "asym_mig", asym_mig, emp_params, fs_folded=True)
model = Models_2D.asym_mig(emp_params, fs.sample_sizes) 
moments.Plotting.plot_2d_comp_multinom(model,fs, vmin=0.4,resid_range=5)

emp_params=[0.5495,2.0355,0.47,1.2997,0.749,0.4676,5.7103,3.3224]
scaled_fs = Optimize_Functions_GOF.Optimize_Empirical(fs, "Pop1Pop2_comp", "sec_contact_asym_mig_size", sec_contact_asym_mig_size, emp_params, fs_folded=True)
model = Models_2D.sec_contact_asym_mig_size(emp_params, fs.sample_sizes) 
moments.Plotting.plot_2d_comp_multinom(model,fs, vmin=0.4,resid_range=5)

emp_params=[0.0737,4.1364,12.6723,0.1255,4.2493,0.0841,0.0128]
scaled_fs = Optimize_Functions_GOF.Optimize_Empirical(fs, "Pop1Pop2_comp", "anc_sym_mig_size", anc_sym_mig_size, emp_params, fs_folded=True)
model = Models_2D.anc_sym_mig_size(emp_params, fs.sample_sizes) 
moments.Plotting.plot_2d_comp_multinom(model,fs, vmin=0.4,resid_range=5)

emp_params=[0.2576,0.4794,1.6468,1.9547,0.712,0.0101]
scaled_fs = Optimize_Functions_GOF.Optimize_Empirical(fs, "Pop1Pop2_comp", "anc_asym_mig", anc_asym_mig, emp_params, fs_folded=True)
model = Models_2D.anc_asym_mig(emp_params, fs.sample_sizes) 
moments.Plotting.plot_2d_comp_multinom(model,fs, vmin=0.4,resid_range=5)

emp_params=[0.0377,0.172,0.0193]
scaled_fs = Optimize_Functions_GOF.Optimize_Empirical(fs, "Pop1Pop2_comp", "no_mig", no_mig, emp_params, fs_folded=True)
model = Models_2D.no_mig(emp_params, fs.sample_sizes) 
moments.Plotting.plot_2d_comp_multinom(model,fs, vmin=0.4,resid_range=5)

emp_params=[1.8084,9.9994,0.6358,2.0317,2.9654,3.0663,22.1418,0.3006]
scaled_fs = Optimize_Functions_GOF.Optimize_Empirical(fs, "Pop1Pop2_comp", "anc_asym_mig_size", anc_asym_mig_size, emp_params, fs_folded=True)
model = Models_2D.anc_asym_mig_size(emp_params, fs.sample_sizes) 
moments.Plotting.plot_2d_comp_multinom(model,fs, vmin=0.4,resid_range=5)

emp_params=[4.5506,14.4216,0.1023,0.3513,0.0553,0.0518]
scaled_fs = Optimize_Functions_GOF.Optimize_Empirical(fs, "Pop1Pop2_comp", "no_mig_size", no_mig_size, emp_params, fs_folded=True)
model = Models_2D.no_mig_size(emp_params, fs.sample_sizes) 
moments.Plotting.plot_2d_comp_multinom(model,fs, vmin=0.4,resid_range=5)

emp_params=[0.5259,1.2377,1.4422,9.8968,0.0737]
scaled_fs = Optimize_Functions_GOF.Optimize_Empirical(fs, "Pop1Pop2_comp", "anc_sym_mig", anc_sym_mig, emp_params, fs_folded=True)
model = Models_2D.anc_sym_mig(emp_params, fs.sample_sizes) 
moments.Plotting.plot_2d_comp_multinom(model,fs, vmin=0.4,resid_range=5)