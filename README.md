# Demographic Modeling

## Demographic modeling using site frequency spectrum of two populations.


- Requires `moments` to be installed (https://bitbucket.org/simongravel/moments/src/master/):
```
pip install git+https://bitbucket.org/simongravel/moments.git
```

- Clone `moments_pipeline`:
```
git clone https://github.com/dportik/moments_pipeline.git
```
# Steps: 

## Part 1A. Allele Frequencies:
Reformat `angsd` allele frequencies to `dadi` SNP format for use with the `moments_pipeline` to calculate SFS and perform demographic model optimizations and comparisons.

- `R`: reformat `angsd` *GlobalSNP_HiC_100DP670_MinInd100_MinQ20.mafs_fst.vep.txt* output to `dadi` SNP format
    - scripts = `1A_MAF2dadiSNPinfile.R`
    - output = *Pop1-Pop2_MAF2dadiSNPinfile*
- `Python`: calculate 2D SFS with `moments`
    - edit scripts = `moments_Run_2D_Set_CJC.py` or `moments_Run_b4models.py`
    - output1 = site frequency spectrum
        - *Pop1_Pop2.2dsfs*

## Part 1B. `angsd` 2D SFS:
Calculate 2D site frequency spectum in `angsd` then read in data to `moments_pipeline`.

- `bash`: calculate 2D SFS in `angsd` from indexed indiviudal SFS *HiRise_rh_HiC_nosites.saf.idx*
    - scripts = `1B_ANGSD-2Dsfs_Pop1-Pop2.sh`
    - output = *2Dsfs.Pop1_Pop2.sfs*
    - edited output = *Pop1_Pop2.2dsfs*

## Perform demographic model optimizations and do model comparisons:

2. `Python`: read in *2Dsfs.Pop1_Pop2.sfs* through `moments` models via `moments_pipeline` scripts
    - edit scripts = `moments_Run_2D_Set_CJC.py` or `moments_Run_b4models.py`
3. `Python`: run *Pop1_Pop2.2dsfs* through `moments` models via `moments_pipeline` scripts
    - edit scripts = `moments_Run_2D_Set_CJC.py`
    - output = all runs and optimized runs for each model
        - *Pop1-Pop2.model.log.txt*
        - *Pop1-Pop2.model.optimized.txt*
4. `R`: compare and find best model
    - scripts = `4_CompareMomentsModels.R`
    - output = a table of the best AIC for each model, the best model, and mutation rate relative to N<small>E</small>
        - *RESULTS.Pop1_Pop2.AIC_output*
        - *RESULTS.Pop1_Pop2.bestmod*
        - *RESULTS.Pop1_Pop2.mr*
5. `Python`: look at how well the optimized parameters fit the data
    - scripts = `5b_Optimized_params_fits.py`
    - output = 3 heat maps (data, model, residuals) and a residuals bar graph

# Part 1A. Allele Frequencies:
### `dadi` SNP format

- ***column 1***: the in-group reference sequence at that SNP, including flanking bases (unknown bases denoted by -). 
    - The header label is arbitrary.
- ***column 2***: the aligned outgroup reference sequence at that SNP, including flanking bases (unknown bases denoted by -).
    - The header label is arbitrary.
- ***Allele1 column***: first segregating allele. 
    - The header must be exactly *Allele1*.
- subsequent columns: one for each population, each giving the
number of times Allele1 was **observed** in that population. 
    - The header for each column should be the population identifer.
- ***Allele2 column***: second segregating allele. 
    - The column header must be exactly *Allele2*.
- subsequent columns: one for each population, each giving the
number of times Allele2 was **observed** in that population. 
    - The header for each column should be the population identifer and the columns should be in the same order as for the Allele1 entries.
- subsequent columns concatenated with _ to assign a label for each SNP.
    -The Allele1 and Allele2 headers must be exactly those values because the number of columns between those two is used to infer the number of populations in the file.
    
```
Maj     Min     Allele1 DP      SS      Allele2 DP.1    SS.1    Chr     Pos
-G-     -T-     G       26      25      T       0       3       HiC_scaffold_1  6065
-G-     -T-     G       30      32      T       0       10      HiC_scaffold_1  12608
-C-     -T-     C       27      30      T       7       12      HiC_scaffold_1  12612
-T-     -C-     T       22      32      C       8       8       HiC_scaffold_1  12707
-A-     -C-     A       21      31      C       7       5       HiC_scaffold_1  12810
-G-     -A-     G       35      20      A       1       0       HiC_scaffold_1  13227
-T-     -A-     T       33      23      A       1       5       HiC_scaffold_1  13329
-G-     -A-     G       35      23      A       3       7       HiC_scaffold_1  13336
-G-     -A-     G       20      22      A       0       4       HiC_scaffold_1  13402
```

- R script to reformat merged MAF output from `angsd` (ex. below).
```
#!/bin/bash
BAMDIR=$1 # Path to directory of bamfiles
BAMLIST=$2 # Path to textfile listing bamfiles to include in population-level MAF 
REFERENCE=$3 # Path to reference genome 
SNPLIST=$4 #Path to Global SNPList, e.g. 'Global_SNPList_'$OUTBASE'.txt' from global SNP calling step 
OUTDIR=$5 # Path to output files 
OUTBASE=$6 # Basename for output files, e.g. GlobalSNP_100DP646_MinInd100_MinQ20
POPULATION=$7 # Population name
MINDP=$8 # min depth across population
MININD=$9 # min individuals across population

cd ${BAMDIR}

angsd -b ${BAMLIST} -anc ${REFERENCE} -out ${OUTDIR}${POPULATION}'_'${OUTBASE} \
-dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doPost 1 -doVcf 1 \
-doCounts 1 -doDepth 1 -dumpCounts 1 -P 12 -setMinDepth ${MINDP} \
-minInd ${MININD} -sites ${SNPLIST} \
>& ${OUTDIR}${POPULATION}'_'${OUTBASE}'.log'
```

```R
library(stringr)

mafs=read.table('GlobalSNP_HiC_100DP670_MinInd100_MinQ20.mafs_fst.vep.txt',header=T,sep='\t')

#multiply the minor allele frequency by number of alleles (individuals*2) to get observed Allele1
mafs$Pop1=round(mafs$Pop1_MAF*mafs$Pop1_nInd*2)
mafs$Pop2=round(mafs$Pop2_MAF*mafs$Pop2_nInd*2)

#1 minus the minor allele frequency multiplied by number of alleles (individuals*2) to get observed Allele2
mafs$Pop1mj=round((1-mafs$Pop1_MAF)*mafs$Pop1_nInd*2)
mafs$Pop2mj=round((1-mafs$Pop2_MAF)*mafs$Pop2_nInd*2)

#Create and Write file for Pop1-Pop2

sfs=mafs[,c('Chr','Pos','major','minor','Pop1','Pop2','Pop1mj','Pop2mj')]
sfs$Maj=apply(sfs,1,function(x) paste(c('-',x[3],'-'),collapse=''))
sfs$Min=apply(sfs,1,function(x) paste(c('-',x[4],'-'),collapse=''))
names(sfs)[3]='Allele1'
names(sfs)[4]='Allele2'
sfs=sfs[,c(9:10,3,7:8,4:6,1:2)] #Maj, Min,Allele1,pop1maj,pop2maj,Allele2,pop1min,pop2min,contig,pos
names(sfs)[4:5]=c('Pop1','Pop2')
sfs=sfs[complete.cases(sfs),] #eliminates lines with NA alleles

write.table(sfs,file='Pop1-Pop2_MAF2dadiSNPinfile',row.names=F,sep='\t',quote=F)
```

### `moments_pipeline` before models

#### To run moments_pipeline, edit lines marked by `#**************` in `moments_Run_2D_Set.py` in *moments_pipeline/Two_Population_Pipeline/* directory
- Written for Python 2.7 but I have been using 3.8.5 and it seems to work just fine.
- `moments_Run_b4models.py` has only the code needed to run through check.
- Ex:
```python
import sys
import os
import numpy
import moments
import pylab
from datetime import datetime
import Optimize_Functions
import Models_2D

#===========================================================================
# Import data to create joint-site frequency spectrum
#===========================================================================

#**************
snps = "Pop1-Pop2_MAF2dadiSNPinfile"

#Create python dictionary from snps file
dd = moments.Misc.make_data_dict(snps)

#**************
#pop_ids is a list which should match the populations headers of your SNPs file columns
pop_ids=["Pop1", "Pop2"]

#**************
#projection sizes, in ALLELES not individuals
proj = [33, 88]

#Convert this dictionary into folded AFS object
#[polarized = False] creates folded spectrum object
fs = moments.Spectrum.from_data_dict(dd, pop_ids=pop_ids, projections = proj, polarized = False)
```

- print some useful information about the afs or jsfs
```python
print("\n\n============================================================================")
print("\nData for site frequency spectrum\n")
print("Projection: {}".format(proj))
print("Sample sizes: {}".format(fs.sample_sizes))
print("Sum of SFS: {}".format(numpy.around(fs.S(), 2)))
print("\n============================================================================\n")
```

- To write SFS to a file for easy upload later:
```python
fs.to_file("Pop1_Pop2.2dsfs")
```

# Part 1B. 2D SFS:

### `angsd`
- calculate 2D SFS with:
    - nSites = number of sites == 10,000,000 sites
    - P = number of threads  
    
***This takes a LONG time!*** Pop1-SS took 1 week and Pop1-Pop2 took 3+ weeks on 8 threads.
```
#!/bin/sh

realSFS=$1 # Path to realSFS in angsd program directory
dir=$2 # Path to working directory

Pop1=$3 # Path to Pop1_HiRise_rh_HiC_nosites.saf.idx
Pop2=$4 # Path to Pop2_HiRise_rh_HiC_nosites.saf.idx

${realSFS} ${Pop1} ${Pop2} -nSites 10000000 -P 8 > ${dir}2Dsfs.Pop1_Pop2.sfs 2> ${dir}2Dsfs.Pop1_Pop2.log
```
# Perform model comparisons:

- all modules need to be in the working directory
    - instead of copying files into the currect directory or always running it in the `moments_pipeline` directory, make symbolic links directly to the files with `moments-modules-symbolic-links.sh`.
```
#!/bin/sh

DIR=$1 # path to moments_pipeline

ln -s ${DIR}/Two_Population_Pipeline/Models_2D.py Models_2D.py
ln -s ${DIR}/Two_Population_Pipeline/Summarize_Outputs.py Summarize_Outputs.py
ln -s ${DIR}/Optimize_Functions.py Optimize_Functions.py
ln -s ${DIR}/Goodness_of_Fit/Optimize_Functions_GOF.py Optimize_Functions_GOF.py
```

## Part 2. `moments_pipeline` before models

### To run moments_pipeline, edit lines marked by `#**************` in `moments_Run_2D_Set.py` in *moments_pipeline/Two_Population_Pipeline/* directory
- `moments_Run_b4models.py` has only the code needed to run through check.
```python
import sys
import os
import numpy
import moments
import pylab
from datetime import datetime
import Optimize_Functions
import Models_2D
```

- Assign population names:
```python
pop_ids=["Pop1", "Pop2"]
```

- Read in SFS file:
```python
fs = moments.Spectrum.from_file("Pop1_Pop2.2dsfs")
```

- Project to different sample size:
```python
proj = [33, 33]
fs=fs.project(proj)
```

- Fold unfolded SFS:
```python
fs=fs.fold()
```

- print some useful information about the afs or jsfs
```python
print("\nData for site frequency spectrum\n")
print("Projection: {}".format(proj))
print("Sample sizes: {}".format(fs.sample_sizes))
print("Sum of SFS: {}".format(numpy.around(fs.S(), 2)))
```

- print a heatmap of SFS data
```python
moments.Plotting.plot_single_2d_sfs(fs, vmin=5)
pylab.show()
```

## 3. run models with `moments_pipeline`

### run live on command line or as script:

- to run as script:
```
python edited_moments_Run_2D_Set.py
```

- `moments_Run_2D_Set_CJC.py` includes all models of the Diversification Model Set with bounds changed from default.
    - To use default, comment out or remove `upper = []` and `lower = []` and remove `in_upper=upper, in_lower=lower,` from `Optimize_Functions.Optimize_Routine`.
```python
#================================================================================
# Calling external 2D models from the Models_2D.py script with Altered Bounds
#================================================================================

#create a prefix based on the population names to label the output files
#calls population names assigned above "pop_ids="
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
#True if "fs=fs.fold()" is run earlier
fs_folded = True

#**************
# Split into two populations, no migration.
upper = [30, 30, 30]
lower = [1e-05, 1e-05, 1e-05]
Optimize_Functions.Optimize_Routine(fs, prefix, "no_mig", Models_2D.no_mig, rounds, 3, fs_folded=fs_folded,
                                        reps=reps, maxiters=maxiters, folds=folds, in_upper=upper, in_lower=lower, param_labels = "nu1, nu2, T")

```

## 4. What's the best model?

- Use `4_CompareMomentsModels.R`

- Set working directory
```R
setwd('~/File/Path/to/Directory/Pop1_Pop2/')
```
- assign list of models
```R
models=c('no_mig','sym_mig','asym_mig','anc_sym_mig','anc_asym_mig','sec_contact_sym_mig','sec_contact_asym_mig','no_mig_size','sym_mig_size','asym_mig_size','anc_sym_mig_size','anc_asym_mig_size','sec_contact_sym_mig_size','sec_contact_asym_mig_size','sym_mig_twoepoch','asym_mig_twoepoch','sec_contact_sym_mig_three_epoch','sec_contact_asym_mig_three_epoch','sec_contact_sym_mig_size_three_epoch','sec_contact_asym_mig_size_three_epoch')
```

- This loop reads in all the optim files for the models above
```R
prefix='Pop1_Pop2'
fits=NULL

for (model in models){
	file=paste(c(prefix,model,'optimized.txt'),collapse='.')
	optim=read.table(file,sep='\t',header=T)
	optim=optim[which(optim$Model!='Model'),]#this line was added to remove multiple headers in the files
	optim[,3:6]=apply(optim[,3:6],2,as.character)
	optim[,3:6]=apply(optim[,3:6],2,as.numeric)
	optim$params=names(optim)[7]
	names(optim)[7]='values'
	optim$AIC=as.numeric(as.character(optim$AIC))
	fits=rbind(fits,optim[which(optim$AIC==min(optim$AIC,na.rm=T)),])
	}
fits=fits[order(fits$AIC),]
```

- write the ordered AIC output to a file
```R
write.table(fits,'RESULTS.Pop1_Pop2.AIC_output', quote = FALSE, sep = "\t",row.names = FALSE, col.names = TRUE)
```

- this part will estimate the parameters in demographic units for the best model
```R
bestmod=fits[1,]
write.table(bestmod,'RESULTS.Pop1_Pop2.bestmod', quote = FALSE, sep = "\t",row.names = FALSE, col.names = TRUE)

params=as.numeric(unlist(strsplit(as.character(bestmod$values),",",)))
paramnames=unlist(strsplit(as.character(bestmod$params),"\\.\\."))
paramnames[1]=unlist(strsplit(paramnames[1],"\\."))[2]
paramnames[length(paramnames)]=unlist(strsplit(paramnames[length(paramnames)],"\\."))

type=sapply(strsplit(paramnames,split=''),function(x) x[1])
theta=bestmod$theta
```

- assign the appropriate value in the next three lines
```R
mu=2.96e-09 #mutation rate
s=5319629.82 #sum of SFS (from the line in moments)
seqlgth=2212099196*(s/18381387) #(total length of genome)*(s/(SNPs in genome))
nref=theta/(4*mu*seqlgth)
```
- convert values into regular units 
```R
tparams=params
tparams[which(type=='n')]=params[which(type=='n')]*nref
tparams[which(type=='T')]=params[which(type=='T')]*nref*2
mr=params[4]/(nref*2)#migration rate
paramnames #names of parameters
tparams #values in regular units
mr #migration rate (relative to Ne)
```

## 5. Does it fit the data?

- Find and use projection with the largest amount of SNPs.
    - `5a_FindProjection.py`.
```python
import sys
import os
import numpy
import moments
import pylab
import pandas as pd
from datetime import datetime

mylist=[]

for x in range(65,1,-1): # range starts with number of total alleles-1
	fs = moments.Spectrum.from_file("Pop1_Pop2.2dsfs")
	proj = [x, x]
	fs=fs.project(proj)
	mylist.append(format(numpy.around(fs.S(), 2)))

# create a pandas dataframe using pd.DataFrame function with a dictionary with two lists as input.
df = pd.DataFrame({'proj':range(65,1,-1), 'SNPs':mylist})

print(df)
df.to_csv(r'find-projection.csv', index = False)

``` 

### Best run live from `5b_Optimized_params_fits.py`, edit fields marked by `#**************`
- Ex.
```python
import sys
import os
import numpy
import moments
import pylab
from datetime import datetime
import Optimize_Functions_GOF
import Models_2D
from Models_2D import sec_contact_sym_mig_size

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

#**************
emp_params=[0.1129,0.2103,2.9732,1.4448,0.1397]
scaled_fs = Optimize_Functions_GOF.Optimize_Empirical(fs, "bestmod","asym_mig",asym_mig, emp_params, fs_folded=fs_folded)
model = Models_2D.asym_mig(emp_params, fs.sample_sizes)
moments.Plotting.plot_2d_comp_multinom(model,fs, vmin=0.4,resid_range=5)#,resid_range=3
```

### things that need to be changed in last block for each model
- copy parameters from `CompareMomentsModels.R` output
```python
emp_params=[0.1129,0.2103,2.9732,1.4448,0.1397]
```
- generate sfs from optimized params
    - designate prefix for output naming
    - "model_name" is to help label the output files
    - model_name (no "") will access the model function so it should match the model name in `Models_2D.py`
    - don't change fs, emp_params, fs_folded=True
```python
scaled_fs = Optimize_Functions_GOF.Optimize_Empirical(fs, "outfile", "model_name", model_name, emp_params, fs_folded=True)
```
- change model_name
    - may need to import the model from the module
        - `from Models_2D import model_name`
```python
model = Models_2D.model_name(emp_params, fs.sample_sizes) 
```
- plot the observed vs. fit SFS with residuals between observed and fit
    - need to adjust vmin and resid-range for your data
```python
moments.Plotting.plot_2d_comp_multinom(model,fs, vmin=0.4,resid_range=5)
```

# Models
1.	no_mig 
    - Split into two populations, no migration.
```python
# Call external 2D model 
Optimize_Functions.Optimize_Routine(fs, prefix, "no_mig", Models_2D.no_mig, rounds, 3, fs_folded=fs_folded,
                                        reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, T")
# Optimizated params
emp_params=[0.0377,0.172,0.0193]
scaled_fs = Optimize_Functions_GOF.Optimize_Empirical(fs, "Pop1Pop2_comp", "no_mig", no_mig, emp_params, fs_folded=True)
model = Models_2D.no_mig(emp_params, fs.sample_sizes) 
moments.Plotting.plot_2d_comp_multinom(model,fs, vmin=0.4,resid_range=5)
                                      
```
2.	sym_mig 
    - Split into two populations, with continuous symmetric migration.
```python
# Call external 2D model 
Optimize_Functions.Optimize_Routine(fs, prefix, "sym_mig", Models_2D.sym_mig, rounds, 4, fs_folded=fs_folded,
                                        reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m, T")
# Optimizated params
emp_params=[0.0869,0.2795,3.1256,0.0967]
scaled_fs = Optimize_Functions_GOF.Optimize_Empirical(fs, "Pop1Pop2_comp","sym_mig",sym_mig, emp_params, fs_folded=fs_folded)
model = Models_2D.sym_mig(emp_params, fs.sample_sizes)
moments.Plotting.plot_2d_comp_multinom(model,fs, vmin=0.4,resid_range=5)
```
3.	asym_mig 
    - Split into two populations, with continuous asymmetric migration.
```python
# Call external 2D model 
Optimize_Functions.Optimize_Routine(fs, prefix, "asym_mig", Models_2D.asym_mig, rounds, 5, fs_folded=fs_folded,
                                        reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m12, m21, T")
# Optimizated params
emp_params=[0.1564,0.3588,1.6215,2.8602,0.2312]
scaled_fs = Optimize_Functions_GOF.Optimize_Empirical(fs, "Pop1Pop2_comp", "asym_mig", asym_mig, emp_params, fs_folded=True)
model = Models_2D.asym_mig(emp_params, fs.sample_sizes) 
moments.Plotting.plot_2d_comp_multinom(model,fs, vmin=0.4,resid_range=5)
```
4.	anc_sym_mig 
    - Split with continuous symmetric migration, followed by isolation.
```python
# Call external 2D model 
Optimize_Functions.Optimize_Routine(fs, prefix, "anc_sym_mig", Models_2D.anc_sym_mig, rounds, 5, fs_folded=fs_folded,
                                        reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m, T1, T2")
# Optimizated params
emp_params=[0.5259,1.2377,1.4422,9.8968,0.0737]
scaled_fs = Optimize_Functions_GOF.Optimize_Empirical(fs, "Pop1Pop2_comp", "anc_sym_mig", anc_sym_mig, emp_params, fs_folded=True)
model = Models_2D.anc_sym_mig(emp_params, fs.sample_sizes) 
moments.Plotting.plot_2d_comp_multinom(model,fs, vmin=0.4,resid_range=5)
```
5.	anc_asym_mig 
    - Split with continuous asymmetric migration, followed by isolation.
```python
# Call external 2D model 
Optimize_Functions.Optimize_Routine(fs, prefix, "anc_asym_mig", Models_2D.anc_asym_mig, rounds, 6, fs_folded=fs_folded,
                                        reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m12, m21, T1, T2")
# Optimizated params
emp_params=[0.2576,0.4794,1.6468,1.9547,0.712,0.0101]
scaled_fs = Optimize_Functions_GOF.Optimize_Empirical(fs, "Pop1Pop2_comp", "anc_asym_mig", anc_asym_mig, emp_params, fs_folded=True)
model = Models_2D.anc_asym_mig(emp_params, fs.sample_sizes) 
moments.Plotting.plot_2d_comp_multinom(model,fs, vmin=0.4,resid_range=5)
```
6.	sec_contact_sym_mig 
    - Split with no gene flow, followed by period of continuous symmetrical gene flow.
```python
# Call external 2D model 
Optimize_Functions.Optimize_Routine(fs, prefix, "sec_contact_sym_mig", Models_2D.sec_contact_sym_mig, rounds, 5, fs_folded=fs_folded,
                                        reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m, T1, T2")
# Optimizated params
emp_params=[0.4681,1.391,0.7351,12.0598,3.1693]
scaled_fs = Optimize_Functions_GOF.Optimize_Empirical(fs, "Pop1Pop2_comp", "sec_contact_sym_mig", sec_contact_sym_mig, emp_params, fs_folded=True)
model = Models_2D.sec_contact_sym_mig(emp_params, fs.sample_sizes) 
moments.Plotting.plot_2d_comp_multinom(model,fs, vmin=0.4,resid_range=5)
```
7.	sec_contact_asym_mig 
    - Split with no gene flow, followed by period of continuous asymmetrical gene flow.
```python
# Call external 2D model 
Optimize_Functions.Optimize_Routine(fs, prefix, "sec_contact_asym_mig", Models_2D.sec_contact_asym_mig, rounds, 6, fs_folded=fs_folded,
                                        reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m12, m21, T1, T2")
# Optimizated params
emp_params=[0.0434,0.1443,8.3089,5.886,0.1297,0.3704]
scaled_fs = Optimize_Functions_GOF.Optimize_Empirical(fs, "Pop1Pop2_comp", "sec_contact_asym_mig", sec_contact_asym_mig, emp_params, fs_folded=True)
model = Models_2D.sec_contact_asym_mig(emp_params, fs.sample_sizes) 
moments.Plotting.plot_2d_comp_multinom(model,fs, vmin=0.4,resid_range=5)
```
8.	no_mig_size 
    - Split with no migration, then instantaneous size change with no migration.
```python
# Call external 2D model 
Optimize_Functions.Optimize_Routine(fs, prefix, "no_mig_size", Models_2D.no_mig_size, rounds, 6, fs_folded=fs_folded,
                                        reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1a, nu2a, nu1b, nu2b, T1, T2")
# Optimizated params
emp_params=[4.5506,14.4216,0.1023,0.3513,0.0553,0.0518]
scaled_fs = Optimize_Functions_GOF.Optimize_Empirical(fs, "Pop1Pop2_comp", "no_mig_size", no_mig_size, emp_params, fs_folded=True)
model = Models_2D.no_mig_size(emp_params, fs.sample_sizes) 
moments.Plotting.plot_2d_comp_multinom(model,fs, vmin=0.4,resid_range=5)
```
9.	sym_mig_size 
    - Split with symmetric migration, then instantaneous size change with continuous symmetric migration.
```python
# Call external 2D model 
Optimize_Functions.Optimize_Routine(fs, prefix, "sym_mig_size", Models_2D.sym_mig_size, rounds, 7, fs_folded=fs_folded,
                                        reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1a, nu2a, nu1b, nu2b, m, T1, T2")
# Optimizated params
emp_params=[0.1713,4.9709,0.0913,0.229,3.4997,0.1978,0.0733]
scaled_fs = Optimize_Functions_GOF.Optimize_Empirical(fs, "Pop1Pop2_comp","sym_mig_size",sym_mig_size, emp_params, fs_folded=fs_folded)
model = Models_2D.sym_mig_size(emp_params, fs.sample_sizes)
moments.Plotting.plot_2d_comp_multinom(model,fs, vmin=0.4,resid_range=5)
```
10.	asym_mig_size 
    - Split with different migration rates, then instantaneous size change with continuous asymmetric migration.
```python
# Call external 2D model 
Optimize_Functions.Optimize_Routine(fs, prefix, "asym_mig_size", Models_2D.asym_mig_size, rounds, 8, fs_folded=fs_folded,
                                        reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2")
# Optimizated params
emp_params=[7.1512,1.5002,0.0878,0.3801,4.5254,2.0757,5.4081,1.1072]
scaled_fs = Optimize_Functions_GOF.Optimize_Empirical(fs, "Pop1Pop2_comp", "asym_mig_size", asym_mig_size, emp_params, fs_folded=True)
model = Models_2D.asym_mig_size(emp_params, fs.sample_sizes) 
moments.Plotting.plot_2d_comp_multinom(model,fs, vmin=0.4,resid_range=5)
```
11.	anc_sym_mig_size 
    - Split with continuous symmetrical gene flow, followed by instantaneous size change with no migration.
```python
# Call external 2D model 
Optimize_Functions.Optimize_Routine(fs, prefix, "anc_sym_mig_size", Models_2D.anc_sym_mig_size, rounds, 7, fs_folded=fs_folded,
                                        reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1a, nu2a, nu1b, nu2b, m, T1, T2")
# Optimizated params
emp_params=[0.0737,4.1364,12.6723,0.1255,4.2493,0.0841,0.0128]
scaled_fs = Optimize_Functions_GOF.Optimize_Empirical(fs, "Pop1Pop2_comp", "anc_sym_mig_size", anc_sym_mig_size, emp_params, fs_folded=True)
model = Models_2D.anc_sym_mig_size(emp_params, fs.sample_sizes) 
moments.Plotting.plot_2d_comp_multinom(model,fs, vmin=0.4,resid_range=5)
```
12.	anc_asym_mig_size 
    - Split with continuous asymmetrical gene flow, followed by instantaneous size change with no migration.
```python
# Call external 2D model 
Optimize_Functions.Optimize_Routine(fs, prefix, "anc_asym_mig_size", Models_2D.anc_asym_mig_size, rounds, 8, fs_folded=fs_folded,
                                        reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2")
# Optimizated params
emp_params=[1.8084,9.9994,0.6358,2.0317,2.9654,3.0663,22.1418,0.3006]
scaled_fs = Optimize_Functions_GOF.Optimize_Empirical(fs, "Pop1Pop2_comp", "anc_asym_mig_size", anc_asym_mig_size, emp_params, fs_folded=True)
model = Models_2D.anc_asym_mig_size(emp_params, fs.sample_sizes) 
moments.Plotting.plot_2d_comp_multinom(model,fs, vmin=0.4,resid_range=5)
```
13.	sec_contact_sym_mig_size 
    - Split with no gene flow, followed by instantaneous size change with continuous symmetrical migration.
```python
# Call external 2D model 
Optimize_Functions.Optimize_Routine(fs, prefix, "sec_contact_sym_mig_size", Models_2D.sec_contact_sym_mig_size, rounds, 7, fs_folded=fs_folded,
                                        reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1a, nu2a, nu1b, nu2b, m, T1, T2")
# Optimizated params
emp_params=[0.2417,20.6599,0.0917,0.1644,4.8612,0.2341,0.0406]
scaled_fs = Optimize_Functions_GOF.Optimize_Empirical(fs, "Pop1Pop2_comp","sec_contact_sym_mig_size",sec_contact_sym_mig_size, emp_params, fs_folded=fs_folded) #(fs, "outfile", "model_name", model_name, emp_params, fs_folded=True)
model = Models_2D.sec_contact_sym_mig_size(emp_params, fs.sample_sizes) 
moments.Plotting.plot_2d_comp_multinom(model,fs, vmin=0.4,resid_range=5)
```
14.	sec_contact_asym_mig_size 
    - Split with no gene flow, followed by instantaneous size change with continuous asymmetrical migration.
```python
# Call external 2D model 
Optimize_Functions.Optimize_Routine(fs, prefix, "sec_contact_asym_mig_size", Models_2D.sec_contact_asym_mig_size, rounds, 8, fs_folded=fs_folded,
                                        reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2")
# Optimizated params
emp_params=[0.5495,2.0355,0.47,1.2997,0.749,0.4676,5.7103,3.3224]
scaled_fs = Optimize_Functions_GOF.Optimize_Empirical(fs, "Pop1Pop2_comp", "sec_contact_asym_mig_size", sec_contact_asym_mig_size, emp_params, fs_folded=True)
model = Models_2D.sec_contact_asym_mig_size(emp_params, fs.sample_sizes) 
moments.Plotting.plot_2d_comp_multinom(model,fs, vmin=0.4,resid_range=5)
```
15.	sym_mig_twoepoch 
    - Split into two populations, with continuous symmetric migration, rate varying across two epochs.
```python
# Call external 2D model 
Optimize_Functions.Optimize_Routine(fs, prefix, "sym_mig_twoepoch", Models_2D.sym_mig_twoepoch, rounds, 6, fs_folded=fs_folded,
                                        reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m1, m2, T1, T2")
```
16.	asym_mig_twoepoch 
    - Split into two populations, with continuous asymmetric migration, rate varying across two epochs.
```python
# Call external 2D model 
Optimize_Functions.Optimize_Routine(fs, prefix, "asym_mig_twoepoch", Models_2D.asym_mig_twoepoch, rounds, 8, fs_folded=fs_folded,
                                        reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m12a, m21a, m12b, m21b, T1, T2")
```
17.	sec_contact_sym_mig_three_epoch 
    - Split with no gene flow, followed by period of continuous symmetrical migration, then isolation.
```python
# Call external 2D model 
Optimize_Functions.Optimize_Routine(fs, prefix, "sec_contact_sym_mig_three_epoch", Models_2D.sec_contact_sym_mig_three_epoch, rounds, 6, fs_folded=fs_folded,
                                        reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m, T1, T2, T3")
```
18.	sec_contact_asym_mig_three_epoch 
    - Split with no gene flow, followed by period of continuous asymmetrical migration, then isolation.
```python
# Call external 2D model 
Optimize_Functions.Optimize_Routine(fs, prefix, "sec_contact_asym_mig_three_epoch", Models_2D.sec_contact_asym_mig_three_epoch, rounds, 7, fs_folded=fs_folded,
                                        reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m12, m21, T1, T2, T3")
```
19.	sec_contact_sym_mig_size_three_epoch 
    - Split with no gene flow, followed by instantaneous size change with continuous symmetrical migration, then isolation.
```python
# Call external 2D model 
Optimize_Functions.Optimize_Routine(fs, prefix, "sec_contact_sym_mig_size_three_epoch", Models_2D.sec_contact_sym_mig_size_three_epoch, rounds, 8, fs_folded=fs_folded,
                                        reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1a, nu2a, nu1b, nu2b, m, T1, T2, T3")
```
20.	sec_contact_asym_mig_size_three_epoch
    - Split with no gene flow, followed by instantaneous size change with continuous asymmetrical migration, then isolation.
```python
# Call external 2D model
Optimize_Functions.Optimize_Routine(fs, prefix, "sec_contact_asym_mig_size_three_epoch", Models_2D.sec_contact_asym_mig_size_three_epoch, rounds, 9, fs_folded=fs_folded,
                                        reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2, T3")
```

# Interpretation

The final output is a plot that shows a heatmap of simulated data and the actual data as well as a heatmap and bar graph of the residuals.

For a well fit model:
- the data and the model heatmaps will be identical.