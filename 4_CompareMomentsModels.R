setwd('')

models=c('no_mig','sym_mig','asym_mig','anc_sym_mig','anc_asym_mig','sec_contact_sym_mig','sec_contact_asym_mig','no_mig_size','sym_mig_size','asym_mig_size','anc_sym_mig_size','anc_asym_mig_size','sec_contact_sym_mig_size','sec_contact_asym_mig_size','sym_mig_twoepoch','asym_mig_twoepoch','sec_contact_sym_mig_three_epoch','sec_contact_asym_mig_three_epoch','sec_contact_sym_mig_size_three_epoch','sec_contact_asym_mig_size_three_epoch')

##This loop reads in all the optim files for the models above
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
####if you want to write the ordered AIC output to a file, uncomment and fill in the next line
#write.table(fits,'filename', other arguments)


#this part will estimate the parameters in demographic units for the best model
bestmod=fits[1,]
params=as.numeric(unlist(strsplit(as.character(bestmod$values),",",)))
paramnames=unlist(strsplit(as.character(bestmod$params),"\\.\\."))
paramnames[1]=unlist(strsplit(paramnames[1],"\\."))[2]
paramnames[length(paramnames)]=unlist(strsplit(paramnames[length(paramnames)],"\\."))

type=sapply(strsplit(paramnames,split=''),function(x) x[1])
theta=bestmod$theta

##Fill in the appropriate value in the next three lines
mu=0.3e-8 #mutation rate
s=25000 #sum of SFS (from the line in moments)
seqlgth=2.1e9*(s/19e6) #(total length of genome)*(s/(SNPs in genome))
nref=theta/(4*mu*seqlgth)

##This will 
tparams=params
tparams[which(type=='n')]=params[which(type=='n')]*nref
tparams[which(type=='T')]=params[which(type=='T')]*nref*2
mr=params[4]/(nref*2)#migration rate
paramnames #names of parameters
tparams #values in regular units
mr #migration rate (relative to Ne)
