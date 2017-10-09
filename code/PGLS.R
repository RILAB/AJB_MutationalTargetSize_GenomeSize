library(phytools)
tre<-read.nexus('Soltis_et_al.treorg') #multilocus tree from Soltis et al. 2011 doi: 10.3732/ajb.1000404
vals<-read.csv('alpha_values.csv')
l<-as.character(vals[,1])
l<-c(l,'Stellaria')
#drops tips we don't need, also makes it a species tree with total tree depth = 1
tre<-drop.tip(tre,which(!tre$tip %in% l))
#change stellaria to scheidia
tre$tip.label[which(tre$tip=='Stellaria')]<-'Schiedea' #rename caryophyllaceae representative
#bind major lineages not in the phylogeny
tr2<-bind.tip(tre,'Boechera',edge.length=0.05,position=0.05,which(tre$tip=='Arabidopsis'))
tr3<-bind.tip(tr2,'Capsella',edge.length=0.025,position=0.025,which(tr2$tip=='Boechera'))
tr4<-bind.tip(tr3,'Sorghum',edge.length=0.025,position=0.025,which(tr3$tip=='Zea'))
#now read in the species names
s<-paste(vals[,1],vals[,2],sep="_")
#now find the genera w/ multiple species sampled
addBranches<-sapply(tr4$tip.label,function(x) length(grep(x,s)))

#now make a little for loop to add a bunch of branches
trAdd<-tr4
order<-addBranches[which(addBranches>1)]
for (x in 1:length(order)){
  for (y in 1:(order[x]-1)){
    trAdd<-bind.tip(trAdd,names(order)[x],edge.length=0.01,position=0.01,which(trAdd$tip==names(order)[x])[1])
  }
}

#now assign the tips a species names:
for ( x in unique(trAdd$tip)){
trAdd$tip.label[trAdd$tip.label %in% x]<- s[grep(x,s)]
}
#get rid of 0 edges and ultrametricize
trAdd$edge.length<-trAdd$edge.length+0.001
trAdd<-chronopl(trAdd,5e-5)

#now load in the alpha values:
order<-sapply(trAdd$tip.label,function(x) strsplit(x,'_')[[1]][2])
rownames(vals)<-vals[,2]
valsFin<-vals[order,]
rownames(valsFin)<-trAdd$tip.label

#now do pgls
library(nlme)
library(geiger)
glsFit<-gls(derived_alpha ~ log_genome_1C_mb, correlation = corBrownian(phy = trAdd), data=valsFin,method='ML')
plot(valsFin$log_genome_1C_mb,valsFin$derived_alpha,ylab='Derived Alpha',xlab='Log genome size MB')
abline(coef = coef(glsFit))
ft<-lm(derived_alpha ~ log_genome_1C_mb, data=valsFin)
abline(coef=coef(ft),col='red')




