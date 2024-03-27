## format phenotype information for gemma

## read individual IDs, order as in the genotype file
ids<-read.table("IndIds.txt",header=FALSE)
dim(ids)
#[1] 238   1

## read sample and phenotype data
dat<-read.table("sample_sp_loc_host_morph.dsv",header=TRUE)
dim(dat)
#[1] 1169    5

## note orders not same, so need to sort

nn<-rep(NA,238)
ph<-as.data.frame(cbind(nn,nn,nn,nn,nn))
for(i in 1:238){
	a<-which(dat$sample==ids[i,1])
	ph[i,]<-dat[a,]
}

colnames(ph)<-colnames(dat)
mean(ph[,1]==ids[,1])
#[1] 1
## all good

table(ph$morph)
#
#  M   S   U
# 20 120  98

## melanic vs green
melanic<-as.numeric(ph$morph=="M")
mean(melanic)
#[1] 0.08403361

## get index of green (gs and gus)
gg<-which(ph$morph!="M")
phg<-ph[gg,]
pattern<-as.numeric(phg$morph=="S")
mean(pattern)
#[1] 0.5504587

write.table(melanic,file="ph_color.txt",row.names=FALSE,col.names=FALSE,quote=FALSE) # 1 = melanic
write.table(pattern,file="ph_pattern.txt",row.names=FALSE,col.names=FALSE,quote=FALSE) # 1 = striped
write.table(gg,file="pattern_ids.txt",row.names=FALSE,col.names=FALSE,quote=FALSE) ## ids to keep for pattern geno
