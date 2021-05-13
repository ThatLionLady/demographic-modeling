library(stringr)

mafs=read.table('GlobalSNP_HiC_100Pop1670_MinInd100_MinQ20.mafs_fst.vep.txt',header=T,sep='\t')

#multiply the minor allele frequency by number of alleles (individuals*2) to get observed Allele1
mafs$Pop1=round(mafs$Pop1_MAF*mafs$Pop1_nInd*2)
mafs$Pop2=round(mafs$Pop2_MAF*mafs$Pop2_nInd*2)


#1 minus the minor allele frequency multiplied by number of alleles (individuals*2) to get observed Allele2
mafs$Pop1mj=round((1-mafs$Pop1_MAF)*mafs$Pop1_nInd*2)
mafs$Pop2mj=round((1-mafs$Pop2_MAF)*mafs$Pop2_nInd*2)

mafs$Pop2nmj=round((1-mafs$SM_MAF)*mafs$Pop1_nInd*2*(33/88)) #if unequal pop sizes, normalize observations by multiplying by ratio of Pop1/Pop2

#Create and Write file for Pop1-Pop2

sfs=mafs[,c('Chr','Pos','major','minor','Pop1','Pop2','Pop1mj','Pop2mj')]
sfs$Maj=apply(sfs,1,function(x) paste(c('-',x[3],'-'),collapse=''))
sfs$Min=apply(sfs,1,function(x) paste(c('-',x[4],'-'),collapse=''))
names(sfs)[3]='Allele1'
names(sfs)[4]='Allele2'
sfs=sfs[,c(9:10,3,7:8,4:6,1:2)]#Maj, Min,Allele1,pop1maj,pop2maj,Allele2,pop1min,pop2min,contig,pos
names(sfs)[4:5]=c('Pop1','Pop2')

sfs=sfs[complete.cases(sfs),]

write.table(sfs,file='Pop1-Pop2_MAF2dadiSNPinfile',row.names=F,sep='\t',quote=F)