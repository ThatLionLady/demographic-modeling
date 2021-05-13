library(stringr)

mafs=read.table('PPM_GlobalSNP_HiC_100DP670_MinInd100_MinQ20.mafs_fst.vep.txt',header=T,sep='\t')

#multiply the minor allele frequency by number of alleles (individuals*2) to get observed Allele1
mafs$DP=round(mafs$DP_MAF*mafs$DP_nInd*2)
mafs$SS=round(mafs$SS_MAF*mafs$SS_nInd*2)
mafs$SM=round(mafs$SM_MAF*mafs$SM_nInd*2) #direct observations for SM
mafs$SMn=round(mafs$SM_MAF*mafs$DP_nInd*2*(33/88)) #"normalized" observations for SM


#1 minus the minor allele frequency multiplied by number of alleles (individuals*2) to get observed Allele2
mafs$DPmj=round((1-mafs$DP_MAF)*mafs$DP_nInd*2)
mafs$SSmj=round((1-mafs$SS_MAF)*mafs$SS_nInd*2)
mafs$SMmj=round((1-mafs$SM_MAF)*mafs$SM_nInd*2) #direct observations for SM
mafs$SMnmj=round((1-mafs$SM_MAF)*mafs$DP_nInd*2*(33/88)) #normalized observations for SM by multiplying by ratio of DP/SM

#Create and Write file for DP-SS

sfs=mafs[,c('Chr','Pos','major','minor','DP','SS','DPmj','SSmj')]
sfs$Maj=apply(sfs,1,function(x) paste(c('-',x[3],'-'),collapse=''))
sfs$Min=apply(sfs,1,function(x) paste(c('-',x[4],'-'),collapse=''))
names(sfs)[3]='Allele1'
names(sfs)[4]='Allele2'
sfs=sfs[,c(9:10,3,7:8,4:6,1:2)]#Maj, Min,Allele1,pop1maj,pop2maj,Allele2,pop1min,pop2min,contig,pos
names(sfs)[4:5]=c('DP','SS')

sfs=sfs[complete.cases(sfs),]

write.table(sfs,file='DP-SS_MAF2dadiSNPinfile',row.names=F,sep='\t',quote=F)

#Create and Write file for DP-SM

sfs=mafs[,c('Chr','Pos','major','minor','DP','SM','DPmj','SMmj')]
sfs$Maj=apply(sfs,1,function(x) paste(c('-',x[3],'-'),collapse=''))
sfs$Min=apply(sfs,1,function(x) paste(c('-',x[4],'-'),collapse=''))
names(sfs)[3]='Allele1'
names(sfs)[4]='Allele2'
sfs=sfs[,c(9:10,3,7:8,4:6,1:2)]#Maj, Min,Allele1,pop1maj,pop2maj,Allele2,pop1min,pop2min,contig,pos
names(sfs)[4:5]=c('DP','SM')

sfs=sfs[complete.cases(sfs),]

write.table(sfs,file='DP-SM_MAF2dadiSNPinfile',row.names=F,sep='\t',quote=F)

#Create and Write file for DP-SMn

sfs=mafs[,c('Chr','Pos','major','minor','DP','SMn','DPmj','SMnmj')]
sfs$Maj=apply(sfs,1,function(x) paste(c('-',x[3],'-'),collapse=''))
sfs$Min=apply(sfs,1,function(x) paste(c('-',x[4],'-'),collapse=''))
names(sfs)[3]='Allele1'
names(sfs)[4]='Allele2'
sfs=sfs[,c(9:10,3,7:8,4:6,1:2)]#Maj, Min,Allele1,pop1maj,pop2maj,Allele2,pop1min,pop2min,contig,pos
names(sfs)[4:5]=c('DP','SMd')

sfs=sfs[complete.cases(sfs),]

write.table(sfs,file='DP-SMn_MAF2dadiSNPinfile',row.names=F,sep='\t',quote=F)
