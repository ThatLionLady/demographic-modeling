'''
Usage: python moments_Run_2D_Set.py

This is a modified version of the 'moments_Run_Optimizations.py' script in which
we run optimizations for 2D comparisons for a large set of models that have been
made available as part of published work. These models are stored in the
Models_2D.py script, and will be called directly here. The user can delete or
comment out models in order to analyze a subset of the models available. 

This script must be in the same working directory as Optimize_Functions.py, which
contains all the functions necessary, as well as the  Models_2D.py script, which
has all the model definitions.


General workflow:
 The optimization routine runs a user-defined number of rounds, each with a user-defined
 or predefined number of replicates. The starting parameters are initially random, but after
 each round is complete the parameters of the best scoring replicate from that round are
 used to generate perturbed starting parameters for the replicates of the subsequent round.
 The arguments controlling steps of the optimization algorithm (maxiter) and perturbation
 of starting parameters (fold) can be supplied by the user for more control across rounds.
 The user can also supply their own set of initial parameters, or set custom bounds on the
 parameters (upper_bound and lower_bound) to meet specific model needs. This flexibility
 should allow these scripts to be generally useful for model-fitting with any data set.

 
Outputs:
 For each model run, there will be a log file showing the optimization steps per replicate
 and a summary file that has all the important information. Here is an example of the output
 from a summary file, which will be in tab-delimited format:
 
 Model	Replicate	log-likelihood	AIC	chi-squared	theta	optimized_params(nu1, nu2, m, T)
 sym_mig	Round_1_Replicate_1	-1684.99	3377.98	14628.4	383.04	0.2356,0.5311,0.8302,0.182
 sym_mig	Round_1_Replicate_2	-2255.47	4518.94	68948.93	478.71	0.3972,0.2322,2.6093,0.611
 sym_mig	Round_1_Replicate_3	-2837.96	5683.92	231032.51	718.25	0.1078,0.3932,4.2544,2.9936
 sym_mig	Round_1_Replicate_4	-4262.29	8532.58	8907386.55	288.05	0.3689,0.8892,3.0951,2.8496
 sym_mig	Round_1_Replicate_5	-4474.86	8957.72	13029301.84	188.94	2.9248,1.9986,0.2484,0.3688


Notes/Caveats:
 The likelihood and AIC returned represent the true likelihood only if the SNPs are
 unlinked across loci. For ddRADseq data where a single SNP is selected per locus, this
 is true, but if SNPs are linked across loci then the likelihood is actually a composite
 likelihood and using something like AIC is no longer appropriate for model comparisons.
 See the discussion group for more information on this subject. 


Citations:
 If you use these moments-based scripts, please
 cite the following publication:
    Leache, A.D., Portik, D.M., Rivera, D., Rodel, M.-O., Penner, J.,
    Gvozdik, V., Greenbaum, E., Jongsma, G.F.M., Ofori-Boateng, C., 
    Burger, M., Eniang, E.A., Bell, R.C., and M.K. Fujita. 2019.
    Exploring rain forest diversification using demographic model testing
    in the African foam-nest tree frog Chiromantis rufescens. Journal of
    Biogeography, Early View. doi: 10.1111/jbi.13716

 If you use the diversification models for your work, please
 cite the following publications:
    Portik, D.M., Leache, A.D., Rivera, D., Blackburn, D.C., Rodel, M.-O.,
    Barej, M.F., Hirschfeld, M., Burger, M., and M.K. Fujita. 2017.
    Evaluating mechanisms of diversification in a Guineo-Congolian forest
    frog using demographic model selection. Molecular Ecology 26: 5245-5263.
    doi: 10.1111/mec.14266

    Charles, K.C., Bell, R.C., Blackburn, D.C., Burger, M., Fujita, M.K., 
    Gvozdik, V., Jongsma, G.F.M., Leache, A.D., and D.M. Portik. 2018. 
    Sky, sea, and forest islands: diversification in the African leaf-folding 
    frog Afrixalus paradorsalis (Order: Anura, Family: Hyperoliidae). Journal 
    of Biogeography 45: 1781-1794. doi: 10.1111/jbi.13365     
   
 If you are interesting in contributing your models to this workflow, please email me!

-------------------------
Written for Python 2.7 <-- BEEN USING WITH PYTHON 3.8.5 JUST FINE (CJC 2021)
Python modules required:
-Numpy
-Scipy
-dadi
-------------------------

Daniel Portik
daniel.portik@gmail.com
https://github.com/dportik
Updated August 2019
'''

import sys
import os
import numpy
import moments
import pylab
from datetime import datetime
import Optimize_Functions
import Models_2D


#=========================================================================================
# Part 1A: Import data to create joint-site frequency spectrum (comment out if not using)
#=========================================================================================

#**************
snps = "DP-SM_MAF2dadiSNPinfile"

#Create python dictionary from snps file
dd = moments.Misc.make_data_dict(snps)

#**************
#pop_ids is a list which should match the populations headers of your SNPs file columns
pop_ids=["DP", "SS"]

#**************
#projection sizes, in ALLELES not individuals
proj = [33, 33]

#Convert this dictionary into folded AFS object
#[polarized = False] creates folded spectrum object
fs = moments.Spectrum.from_data_dict(dd, pop_ids=pop_ids, projections = proj, polarized = False)

#**************
#write to file for later use
fs.to_file("DP_SS.2dsfs")

#project to different sample size
fs=fs.project(proj)


#==========================================================================================================================
# Part 2: Import site frequency spectrumn from angsd OR reading in previously made file above (comment out if not using)
#==========================================================================================================================


#**************
#pop_ids is a list which should match the populations headers of your SNPs file columns
pop_ids=["DP", "SS"]

#**************
#projection sizes, in ALLELES not individuals
proj = [33, 33]

#**************
#•Read in SFS file generated by angsd
fs = moments.Spectrum.from_file("DP_SS.2dsfs")

#project to different sample size
fs=fs.project(proj)

#fold unfolded SFS
fs=fs.fold()


#===========================================================================
# Print some useful information
#===========================================================================

print("\n\n============================================================================")
print("\nData for site frequency spectrum\n")
print("Projection: {}".format(proj))
print("Sample sizes: {}".format(fs.sample_sizes))
print("Sum of SFS: {}".format(numpy.around(fs.S(), 2)))
print("\n============================================================================\n")


#===========================================================================
# Plot the data before running models
#===========================================================================

moments.Plotting.plot_single_2d_sfs(fs, vmin=5)
pylab.show()

#If running as a script, this will pause the run as long as the plot is open. Close plot to continue.


#=======================================================================================
# Part 3: Calling external 2D models from the Models_2D.py script with Altered Bounds
#=======================================================================================
'''
 We will use a function from the Optimize_Functions.py script for our optimization routines:
  
 Optimize_Routine(fs, outfile, model_name, func, rounds, param_number, fs_folded=True,
                       reps=None, maxiters=None, folds=None, in_params=None, 
                       in_upper=None, in_lower=None, param_labels=" "):

   Mandatory Arguments =
    fs:  spectrum object name
    outfile:  prefix for output naming
    model_name: a label to slap on the output files; ex. "no_mig"
    func: access the model function from within 'moments_Run_Optimizations.py' or 
          from a separate python model script, ex. after importing Models_2D, calling Models_2D.no_mig
    rounds: number of optimization rounds to perform
    param_number: number of parameters in the model selected (can count in params line for the model)
    fs_folded: A Boolean value (True or False) indicating whether the empirical fs is folded (True) or not (False).

   Optional Arguments =
     reps: a list of integers controlling the number of replicates in each of the optimization rounds
     maxiters: a list of integers controlling the maxiter argument in each of the optimization rounds
     folds: a list of integers controlling the fold argument when perturbing input parameter values
     in_params: a list of parameter values 
     in_upper: a list of upper bound values
     in_lower: a list of lower bound values
     param_labels: list of labels for parameters that will be written to the output file to keep track of their order

Below, I give all the necessary information to call each model available in the
Models_2D.py script. I have set the optimization routine to be the same for each
model using the optional lists below, which are included as optional arguments for
each model. This particular configuration will run 4 rounds as follows:
Round1 - 10 replicates, maxiter = 3, fold = 3
Round2 - 20 replicates, maxiter = 5, fold = 2
Round3 - 30 replicates, maxiter = 10, fold = 2
Round4 - 40 replicates, maxiter = 15, fold = 1

If this script was run as is, each model would be called and optimized sequentially;
this could take a very long time. For your actual analyses, I strongly recommend
creating multiple scripts with only a few models each and running them
independently.

'''

#create a prefix based on the population names to label the output files
#ex. Pop1_Pop2
prefix = "_".join(pop_ids)


#**************
#Set the number of rounds here
rounds = 4

#define the lists for optional arguments
#you can change these to alter the settings of the optimization routine
reps = [10,20,30,40]
maxiters = [3,5,10,15]
folds = [3,2,2,1]


#**************
#Indicate whether your frequency spectrum object is folded (True) or unfolded (False)
fs_folded = True


#**************
# To change bounds, edit upper and lower. 
# To use default, comment out upper and lower and remove "in_upper=upper, in_lower=lower," from Optimize_Routine.


# Split into two populations, no migration.
upper = [30, 30, 30]
lower = [1e-05, 1e-05, 1e-05]
Optimize_Functions.Optimize_Routine(fs, prefix, "no_mig", Models_2D.no_mig, rounds, 3, fs_folded=fs_folded,
                                        reps=reps, maxiters=maxiters, folds=folds, in_upper=upper, in_lower=lower, param_labels = "nu1, nu2, T")


# Split into two populations, with continuous symmetric migration.
upper = [30, 30, 30, 30]
lower = [1e-05, 1e-05, 0, 1e-05]
Optimize_Functions.Optimize_Routine(fs, prefix, "sym_mig", Models_2D.sym_mig, rounds, 4, fs_folded=fs_folded,
                                        reps=reps, maxiters=maxiters, folds=folds, in_upper=upper, in_lower=lower, param_labels = "nu1, nu2, m, T")


# Split into two populations, with continuous asymmetric migration.
upper = [30, 30, 30, 30, 30]
lower = [1e-05, 1e-05, 0, 0, 1e-05]
Optimize_Functions.Optimize_Routine(fs, prefix, "asym_mig", Models_2D.asym_mig, rounds, 5, fs_folded=fs_folded,
                                        reps=reps, maxiters=maxiters, folds=folds, in_upper=upper, in_lower=lower, param_labels = "nu1, nu2, m12, m21, T")


# Split with continuous symmetric migration, followed by isolation.
upper = [30, 30, 30, 30, 30]
lower = [1e-05, 1e-05, 0, 1e-05, 1e-05]
Optimize_Functions.Optimize_Routine(fs, prefix, "anc_sym_mig", Models_2D.anc_sym_mig, rounds, 5, fs_folded=fs_folded,
                                        reps=reps, maxiters=maxiters, folds=folds, in_upper=upper, in_lower=lower, param_labels = "nu1, nu2, m, T1, T2")


# Split with continuous asymmetric migration, followed by isolation.
upper = [30, 30, 30, 30, 30, 30]
lower = [1e-05, 1e-05, 0, 0, 1e-05, 1e-05]
Optimize_Functions.Optimize_Routine(fs, prefix, "anc_asym_mig", Models_2D.anc_asym_mig, rounds, 6, fs_folded=fs_folded,
                                        reps=reps, maxiters=maxiters, folds=folds, in_upper=upper, in_lower=lower, param_labels = "nu1, nu2, m12, m21, T1, T2")


# Split with no gene flow, followed by period of continuous symmetrical gene flow.
upper = [30, 30, 30, 30, 30]
lower = [1e-05, 1e-05, 0, 1e-05, 1e-05]
Optimize_Functions.Optimize_Routine(fs, prefix, "sec_contact_sym_mig", Models_2D.sec_contact_sym_mig, rounds, 5, fs_folded=fs_folded,
                                        reps=reps, maxiters=maxiters, folds=folds, in_upper=upper, in_lower=lower, param_labels = "nu1, nu2, m, T1, T2")


# Split with no gene flow, followed by period of continuous asymmetrical gene flow.
upper = [30, 30, 30, 30, 30, 30]
lower = [1e-05, 1e-05, 0, 0, 1e-05, 1e-05]
Optimize_Functions.Optimize_Routine(fs, prefix, "sec_contact_asym_mig", Models_2D.sec_contact_asym_mig, rounds, 6, fs_folded=fs_folded,
                                        reps=reps, maxiters=maxiters, folds=folds, in_upper=upper, in_lower=lower, param_labels = "nu1, nu2, m12, m21, T1, T2")


# Split with no migration, then instantaneous size change with no migration.
upper = [30, 30, 30, 30, 30, 30]
lower = [1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05]
Optimize_Functions.Optimize_Routine(fs, prefix, "no_mig_size", Models_2D.no_mig_size, rounds, 6, fs_folded=fs_folded,
                                        reps=reps, maxiters=maxiters, folds=folds, in_upper=upper, in_lower=lower, param_labels = "nu1a, nu2a, nu1b, nu2b, T1, T2")


# Split with symmetric migration, then instantaneous size change with continuous symmetric migration.
upper = [30, 30, 30, 30, 30, 30, 30]
lower = [1e-05, 1e-05, 1e-05, 1e-05, 0, 1e-05, 1e-05]
Optimize_Functions.Optimize_Routine(fs, prefix, "sym_mig_size", Models_2D.sym_mig_size, rounds, 7, fs_folded=fs_folded,
                                        reps=reps, maxiters=maxiters, folds=folds, in_upper=upper, in_lower=lower, param_labels = "nu1a, nu2a, nu1b, nu2b, m, T1, T2")


# Split with different migration rates, then instantaneous size change with continuous asymmetric migration.
upper = [30, 30, 30, 30, 30, 30, 30, 30]
lower = [1e-05, 1e-05, 1e-05, 1e-05, 0, 0, 1e-05, 1e-05]
Optimize_Functions.Optimize_Routine(fs, prefix, "asym_mig_size", Models_2D.asym_mig_size, rounds, 8, fs_folded=fs_folded,
                                        reps=reps, maxiters=maxiters, folds=folds, in_upper=upper, in_lower=lower, param_labels = "nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2")


# Split with continuous symmetrical gene flow, followed by instantaneous size change with no migration.
upper = [30, 30, 30, 30, 30, 30, 30]
lower = [1e-05, 1e-05, 1e-05, 1e-05, 0, 1e-05, 1e-05]
Optimize_Functions.Optimize_Routine(fs, prefix, "anc_sym_mig_size", Models_2D.anc_sym_mig_size, rounds, 7, fs_folded=fs_folded,
                                        reps=reps, maxiters=maxiters, folds=folds, in_upper=upper, in_lower=lower, param_labels = "nu1a, nu2a, nu1b, nu2b, m, T1, T2")


# Split with continuous asymmetrical gene flow, followed by instantaneous size change with no migration.
upper = [30, 30, 30, 30, 30, 30, 30, 30]
lower = [1e-05, 1e-05, 1e-05, 1e-05, 0, 0, 1e-05, 1e-05]
Optimize_Functions.Optimize_Routine(fs, prefix, "anc_asym_mig_size", Models_2D.anc_asym_mig_size, rounds, 8, fs_folded=fs_folded,
                                        reps=reps, maxiters=maxiters, folds=folds, in_upper=upper, in_lower=lower, param_labels = "nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2")


# Split with no gene flow, followed by instantaneous size change with continuous symmetrical migration.
upper = [30, 30, 30, 30, 30, 30, 30]
lower = [1e-05, 1e-05, 1e-05, 1e-05, 0, 1e-05, 1e-05]
Optimize_Functions.Optimize_Routine(fs, prefix, "sec_contact_sym_mig_size", Models_2D.sec_contact_sym_mig_size, rounds, 7, fs_folded=fs_folded,
                                        reps=reps, maxiters=maxiters, folds=folds, in_upper=upper, in_lower=lower, param_labels = "nu1a, nu2a, nu1b, nu2b, m, T1, T2")


# Split with no gene flow, followed by instantaneous size change with continuous asymmetrical migration.
upper = [30, 30, 30, 30, 30, 30, 30, 30]
lower = [1e-05, 1e-05, 1e-05, 1e-05, 0, 0, 1e-05, 1e-05]
Optimize_Functions.Optimize_Routine(fs, prefix, "sec_contact_asym_mig_size", Models_2D.sec_contact_asym_mig_size, rounds, 8, fs_folded=fs_folded,
                                        reps=reps, maxiters=maxiters, folds=folds, in_upper=upper, in_lower=lower, param_labels = "nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2")


# These are the more complex models:

# Split into two populations, with continuous symmetric migration, rate varying across two epochs.
upper = [30, 30, 30, 30, 30, 30]
lower = [1e-05, 1e-05, 0, 0, 1e-05, 1e-05]
Optimize_Functions.Optimize_Routine(fs, prefix, "sym_mig_twoepoch", Models_2D.sym_mig_twoepoch, rounds, 6, fs_folded=fs_folded,
                                        reps=reps, maxiters=maxiters, folds=folds, in_upper=upper, in_lower=lower, param_labels = "nu1, nu2, m1, m2, T1, T2")


# Split into two populations, with continuous asymmetric migration, rate varying across two epochs.
upper = [30, 30, 30, 30, 30, 30, 30, 30]
lower = [1e-05, 1e-05, 0, 0, 0, 0, 1e-05, 1e-05]
Optimize_Functions.Optimize_Routine(fs, prefix, "asym_mig_twoepoch", Models_2D.asym_mig_twoepoch, rounds, 8, fs_folded=fs_folded,
                                        reps=reps, maxiters=maxiters, folds=folds, in_upper=upper, in_lower=lower, param_labels = "nu1, nu2, m12a, m21a, m12b, m21b, T1, T2")


# Split with no gene flow, followed by period of continuous symmetrical migration, then isolation.
upper = [30, 30, 30, 30, 30, 30]
lower = [1e-05, 1e-05, 0, 1e-05, 1e-05, 1e-05]
Optimize_Functions.Optimize_Routine(fs, prefix, "sec_contact_sym_mig_three_epoch", Models_2D.sec_contact_sym_mig_three_epoch, rounds, 6, fs_folded=fs_folded,
                                        reps=reps, maxiters=maxiters, folds=folds, in_upper=upper, in_lower=lower, param_labels = "nu1, nu2, m, T1, T2, T3")


# Split with no gene flow, followed by period of continuous asymmetrical migration, then isolation.
upper = [30, 30, 30, 30, 30, 30, 30]
lower = [1e-05, 1e-05, 0, 0, 1e-05, 1e-05, 1e-05]
Optimize_Functions.Optimize_Routine(fs, prefix, "sec_contact_asym_mig_three_epoch", Models_2D.sec_contact_asym_mig_three_epoch, rounds, 7, fs_folded=fs_folded,
                                        reps=reps, maxiters=maxiters, folds=folds, in_upper=upper, in_lower=lower, param_labels = "nu1, nu2, m12, m21, T1, T2, T3")


# Split with no gene flow, followed by instantaneous size change with continuous symmetrical migration, then isolation.
upper = [30, 30, 30, 30, 30, 30, 30, 30]
lower = [1e-05, 1e-05, 1e-05, 1e-05, 0, 1e-05, 1e-05, 1e-05]
Optimize_Functions.Optimize_Routine(fs, prefix, "sec_contact_sym_mig_size_three_epoch", Models_2D.sec_contact_sym_mig_size_three_epoch, rounds, 8, fs_folded=fs_folded,
                                        reps=reps, maxiters=maxiters, folds=folds, in_upper=upper, in_lower=lower, param_labels = "nu1a, nu2a, nu1b, nu2b, m, T1, T2, T3")


# Split with no gene flow, followed by instantaneous size change with continuous asymmetrical migration, then isolation.
upper = [30, 30, 30, 30, 30, 30, 30, 30, 30]
lower = [1e-05, 1e-05, 1e-05, 1e-05, 0, 0, 1e-05, 1e-05, 1e-05]
Optimize_Functions.Optimize_Routine(fs, prefix, "sec_contact_asym_mig_size_three_epoch", Models_2D.sec_contact_asym_mig_size_three_epoch, rounds, 9, fs_folded=fs_folded,
                                        reps=reps, maxiters=maxiters, folds=folds, in_upper=upper, in_lower=lower, param_labels = "nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2, T3")