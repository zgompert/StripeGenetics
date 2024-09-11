## format phenotype information for gemma

## read individual IDs, order as in the genotype file
ids<-read.table("IndIds.txt",header=FALSE)
dim(ids)
#[1] 602   1

## read sample and phenotype data
dat<-read.table("2013_FHA_All_raw_phenotypes_plus_with_resids.csv",header=TRUE,sep=",")
dim(dat)
#[1] 592    11

## all inds with trait data in genetic set, but 10 missing from traits
tdat<-which(ids[,1] %in% dat$ID)
## order same, just extras
mean(ids[tdat,1]==dat$ID)
#[1] 1

## get pattern + genetic ids, tdat has all with trait, good for color mapping
## note orders not same, so need to sort
pdat<-tdat[which(dat$colour %in% c("Green","Green-striped"))] ## relative to genetic data
pids<-which(dat$colour %in% c("Green","Green-striped"))## relative to pheno file

color<-as.numeric(dat$colour %in% c("Gray","Red")) # 1 = melanic, 0 = green
pattern_bin<-as.numeric(dat$colour[pids] == "Green-striped")
pattern<-cbind(pattern_bin,dat$percent_stripe[pids],dat$resid_percent_stripe[pids])
pattern[,2]<-(pattern[,2]-mean(pattern[,2]))/sd(pattern[,2])
pattern[,3]<-(pattern[,3]-mean(pattern[,3]))/sd(pattern[,3])


write.table(color,file="ph_color.txt",row.names=FALSE,col.names=FALSE,quote=FALSE) # 1 = melanic
write.table(pattern,file="ph_pattern.txt",row.names=FALSE,col.names=FALSE,quote=FALSE) # 1 = striped
write.table(pdat,file="pattern_ids.txt",row.names=FALSE,col.names=FALSE,quote=FALSE) ## ids to keep for pattern geno
write.table(tdat,file="color_ids.txt",row.names=FALSE,col.names=FALSE,quote=FALSE) ## ids to keep for color geno
