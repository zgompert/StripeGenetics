library(data.table)

## colinearity plots for all homologous chromsomes
## chrom number, gs,gsr1, gsr2 ... for reference, from Refugio analysis
chtab<-matrix(c(1,8483,12,10,
	2,14640,6,8,
	3,42935,2,1,
	4,42912,1,1,
	5,18722,7,2,
	6,9928,8,5,
	7,10660,10,7,
	8,7748,11,9,
	9,16151,5,4,
	10,14160,4,3,
	11,12033,9,6,
	12,12380,13,12,
	13,14101,3,11),nrow=13,ncol=4,byrow=TRUE)


## read synteny dat
#########################################################################
## gsh1 x gsr1
##############
dat_gsh1_gsr1<-fread("out_cactusStripe_TcrGSH1_TcrGSR1.psl",header=FALSE)
dfdat<-as.data.frame(dat_gsh1_gsr1)

## target = GSH1
## query = GSR1

## verify/identify large scaffolds (all should be) 
xx<-table(dfdat[,14])
gsr1Ch<-names(xx)[xx>500]
xx<-table(dfdat[,10])
gshxx<-table(dfdat[,10])
gsh1Ch<-names(xx)[xx>500]
keep<-(dfdat[,14] %in% gsr1Ch) & (dfdat[,10] %in% gsh1Ch)
subDfdat<-dfdat[keep,]## retains all, 13 CH each

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tc_gsh1<-as.numeric(unlist(strsplit(x=rownames(tab),split="_"))[seq(2,103,8)])
tc_gsr1<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(2,103,8)])
sizes<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(8,104,8)])## sizes for refugio


## normalize  
ntab<-tab
for(i in 1:13){
	ntab[i,]<-ntab[i,]/sum(ntab[i,])
}
pdf("SynTcrisStripe_GSR1_GSH1.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. cristinae (H GS1)",ylab="T. cristinae (R GS1)",cex.lab=1.4)
axis(2,at=seq(0,13,length.out=13)/13,tc_gsr1,las=2)
axis(1,at=seq(0,13,length.out=13)/13,tc_gsh1,las=2)
box()
dev.off()



## colinearity plots for all homologous chromsomes
## chrom number, gs,gsr1, gsh1
chtab<-matrix(c(1,8483,12,13,
	2,14640,6,5,
	3,42935,2,3,
	4,42912,1,1,
	5,18722,7,12,
	6,9928,8,4,
	7,10660,10,10,
	8,7748,11,11,
	9,16151,5,8,
	10,14160,4,7,
	11,12033,9,9,
	12,12380,13,6,
	13,14101,3,2),nrow=13,ncol=4,byrow=TRUE)


pdf("AlnPlotTcris_GSR1_GSH1.pdf",width=11,height=11)
par(mfrow=c(4,4))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:13){
	tcr_gsh1<-grep(x=subDfdat[,10],pattern=paste("old_",chtab[i,4],"_",sep="")) 
	tcr_gsr1<-grep(x=subDfdat[,14],pattern=paste("old_",chtab[i,3],"_",sep=""))
	cc<-tcr_gsh1[tcr_gsh1 %in% tcr_gsr1]
	subd<-subDfdat[cc,]
	xub<-max(subd[,13]);yub<-max(subd[,17])	

	plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='l',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="H GS1",ylab="R GS1")
	title(main=paste("Chromosome",i),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 2:N){
		if(subd[j,9]=="++"){
			lines(subd[j,12:13],subd[j,16:17])
		}
		else{
			lines(subd[j,12:13],sizes[which(tc_gsr1==chtab[i,3])]-subd[j,16:17])
		}
	}
}
dev.off()

## read synteny dat
#########################################################################
## gsh1 x gsh2
##############
dat_gsh1_gsh2<-fread("out_cactusStripe_TcrGSH1_TcrGSH2.psl",header=FALSE)
dfdat<-as.data.frame(dat_gsh1_gsh2)

## verify/identify large scaffolds (all should be) 
xx<-table(dfdat[,14])
gsh1Ch<-names(xx)[xx>500]
xx<-table(dfdat[,10])
gshxx<-table(dfdat[,10])
gsh2Ch<-names(xx)[xx>500]
keep<-(dfdat[,14] %in% gsh1Ch) & (dfdat[,10] %in% gsh2Ch)
subDfdat<-dfdat[keep,]## retains all, 13 CH each

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tc_gsh2<-as.numeric(unlist(strsplit(x=rownames(tab),split="_"))[seq(2,103,8)])
tc_gsh1<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(2,103,8)])
sizes<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(8,104,8)])## sizes for h2


## normalize  
ntab<-tab
for(i in 1:13){
	ntab[i,]<-ntab[i,]/sum(ntab[i,])
}
pdf("SynTcrisStripe_GSH2_GSH1.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. cristinae (H GS2)",ylab="T. cristinae (H GS1)",cex.lab=1.4)
axis(2,at=seq(0,13,length.out=13)/13,tc_gsh1,las=2)
axis(1,at=seq(0,13,length.out=13)/13,tc_gsh2,las=2)
box()
dev.off()



## colinearity plots for all homologous chromsomes
## chrom number, gs,gsr1, gsh1, gsh2
chtab<-matrix(c(1,8483,12,13,13,
	2,14640,6,5,6,
	3,42935,2,3,2,
	4,42912,1,1,1,
	5,18722,7,12,12,
	6,9928,8,4,5,
	7,10660,10,10,8,
	8,7748,11,11,4,
	9,16151,5,8,9,
	10,14160,4,7,7,
	11,12033,9,9,10,
	12,12380,13,6,11,
	13,14101,3,2,3),nrow=13,ncol=5,byrow=TRUE)


pdf("AlnPlotTcris_GSH2_GSH1.pdf",width=11,height=11)
par(mfrow=c(4,4))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:13){
	tcr_gsh1<-grep(x=subDfdat[,14],pattern=paste("old_",chtab[i,4],"_",sep="")) 
	tcr_gsh2<-grep(x=subDfdat[,10],pattern=paste("old_",chtab[i,5],"_",sep=""))
	cc<-tcr_gsh1[tcr_gsh1 %in% tcr_gsh2]
	subd<-subDfdat[cc,]
	xub<-max(subd[,13]);yub<-max(subd[,17])	

	plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='l',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="H GS2",ylab="H GS1")
	title(main=paste("Chromosome",i),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 2:N){
		if(subd[j,9]=="++"){
			lines(subd[j,12:13],subd[j,16:17])
		}
		else{
			lines(subd[j,12:13],sizes[which(tc_gsh1==chtab[i,4])]-subd[j,16:17])
		}
	}
}
dev.off()

## read synteny dat
#########################################################################
## gsh1 x gush1
##############
dat_gsh1_gush1<-fread("out_cactusStripe_TcrGSH1_TcrGUSH1.psl",header=FALSE)
dfdat<-as.data.frame(dat_gsh1_gush1)

## target = GSH1
## query = GUSH1

## verify/identify large scaffolds (all should be) 
xx<-table(dfdat[,14])
gsh1Ch<-names(xx)[xx>500]
xx<-table(dfdat[,10])
gush1Ch<-names(xx)[xx>500]
keep<-(dfdat[,14] %in% gsh1Ch) & (dfdat[,10] %in% gush1Ch)
subDfdat<-dfdat[keep,]

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tc_gush1<-as.numeric(unlist(strsplit(x=rownames(tab),split="_"))[seq(2,52,4)])
tc_gsh1<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(2,103,8)])
sizes<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(8,104,8)])## sizes for gsh1


## normalize  
ntab<-tab
for(i in 1:13){
	ntab[i,]<-ntab[i,]/sum(ntab[i,])
}
pdf("SynTcrisStripe_GUSH1_GSH1.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. cristinae (H GUS1)",ylab="T. cristinae (H GS1)",cex.lab=1.4)
axis(2,at=seq(0,13,length.out=13)/13,tc_gsh1,las=2)
axis(1,at=seq(0,13,length.out=13)/13,tc_gush1,las=2)
box()
dev.off()



## colinearity plots for all homologous chromsomes
## chrom number, gs,gsr1, gsh1, gush1
chtab<-matrix(c(1,8483,12,13,22,
	2,14640,6,5,23,
	3,42935,2,3,16,
	4,42912,1,1,64,
	5,18722,7,12,5,
	6,9928,8,4,11,
	7,10660,10,10,54,
	8,7748,11,11,7,
	9,16151,5,8,46,
	10,14160,4,7,15,
	11,12033,9,9,2,
	12,12380,13,6,1,
	13,14101,3,2,36),nrow=13,ncol=5,byrow=TRUE)


pdf("AlnPlotTcris_GUSH1_GSH1.pdf",width=11,height=11)
par(mfrow=c(4,4))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:13){
	tcr_gush1<-grep(x=subDfdat[,10],pattern=paste("6S_",chtab[i,5],"_",sep="")) 
	tcr_gsh1<-grep(x=subDfdat[,14],pattern=paste("old_",chtab[i,4],"_",sep=""))
	cc<-tcr_gsh1[tcr_gsh1 %in% tcr_gush1]
	subd<-subDfdat[cc,]
	xub<-max(subd[,13]);yub<-max(subd[,17])	

	plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='l',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="H GUS1",ylab="H GS1")
	title(main=paste("Chromosome",i),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 2:N){
		if(subd[j,9]=="++"){
			lines(subd[j,12:13],subd[j,16:17])
		}
		else{
			lines(subd[j,12:13],sizes[which(tc_gsh1==chtab[i,4])]-subd[j,16:17])
		}
	}
}
dev.off()

## SV bounds
i<-8 ## ch8
sz<-sizes[which(tc_gsh1==chtab[i,4])]

## external translocation
aa<-subd[which(sz-subd[,16] < 2e7 & subd[,12] > 7e7 & subd[,9]=="+-"),1:17]
## gs coords
max(sz-aa[,16])
#[1] 14801201 ## upper bound
min(sz-aa[,17])
#[1] 11974 ## really 0
## size =  14,789,227
## gus coords
min(aa[,12])
#[1] 74554677
max(aa[,13])
#[1] 89370299
# size = 14,815,622

## total inverted bit
aa<-subd[which(subd[,9]=="+-" & subd[,12] > 1e7 & subd[,13] < 2e7 & (sz-subd[,16]) > 3e7 & (sz-subd[,16] < 4.1e7)),1:17]
## gs coords
max(sz-aa[,16])
#[1] 39030359
min(sz-aa[,17])
#[1] 32490874
# size = 6,539,485
## gus coords
min(aa[,12])
#[1] 10097162
max(aa[,13])
#[1] 18370148
# size = 8,272,986

## read synteny dat
#########################################################################
## gsh1 x gush2
##############
dat_gsh1_gush2<-fread("out_cactusStripe_TcrGSH1_TcrGUSH2.psl",header=FALSE)
dfdat<-as.data.frame(dat_gsh1_gush2)

## target = GSH1
## query = GUSH2

## verify/identify large scaffolds (all should be) 
xx<-table(dfdat[,14])
gsh1Ch<-names(xx)[xx>500]
xx<-table(dfdat[,10])
gush2Ch<-names(xx)[xx>500]
keep<-(dfdat[,14] %in% gsh1Ch) & (dfdat[,10] %in% gush2Ch)
subDfdat<-dfdat[keep,]

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tc_gush2<-as.numeric(unlist(strsplit(x=rownames(tab),split="_"))[seq(2,52,4)])
tc_gsh1<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(2,103,8)])
sizes<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(8,104,8)])## sizes for gsh1


## normalize  
ntab<-tab
for(i in 1:13){
	ntab[i,]<-ntab[i,]/sum(ntab[i,])
}
pdf("SynTcrisStripe_GUSH2_GSH1.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. cristinae (H GUS2)",ylab="T. cristinae (H GS1)",cex.lab=1.4)
axis(2,at=seq(0,13,length.out=13)/13,tc_gsh1,las=2)
axis(1,at=seq(0,13,length.out=13)/13,tc_gush2,las=2)
box()
dev.off()



## colinearity plots for all homologous chromsomes
## chrom number, gs,gsr1, gsh1, gush2
chtab<-matrix(c(1,8483,12,13,15,
	2,14640,6,5,1,
	3,42935,2,3,3,
	4,42912,1,1,35,
	5,18722,7,12,10,
	6,9928,8,4,44,
	7,10660,10,10,7,
	8,7748,11,11,23,
	9,16151,5,8,21,
	10,14160,4,7,16,
	11,12033,9,9,12,
	12,12380,13,6,36,
	13,14101,3,2,8),nrow=13,ncol=5,byrow=TRUE)


pdf("AlnPlotTcris_GUSH2_GSH1.pdf",width=11,height=11)
par(mfrow=c(4,4))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:13){
	tcr_gush2<-grep(x=subDfdat[,10],pattern=paste("5T_",chtab[i,5],"_",sep="")) 
	tcr_gsh1<-grep(x=subDfdat[,14],pattern=paste("old_",chtab[i,4],"_",sep=""))
	cc<-tcr_gsh1[tcr_gsh1 %in% tcr_gush2]
	subd<-subDfdat[cc,]
	xub<-max(subd[,13]);yub<-max(subd[,17])	

	plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='l',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="H GUS2",ylab="H GS1")
	title(main=paste("Chromosome",i),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 2:N){
		if(subd[j,9]=="++"){
			lines(subd[j,12:13],subd[j,16:17])
		}
		else{
			lines(subd[j,12:13],sizes[which(tc_gsh1==chtab[i,4])]-subd[j,16:17])
		}
	}
}
dev.off()

## SV bounds
i<-8 ## ch8
sz<-sizes[which(tc_gsh1==chtab[i,4])]
tcr_gush2<-grep(x=subDfdat[,10],pattern=paste("5T_",chtab[i,5],"_",sep=""))
tcr_gsh1<-grep(x=subDfdat[,14],pattern=paste("old_",chtab[i,4],"_",sep=""))
cc<-tcr_gsh1[tcr_gsh1 %in% tcr_gush2]
subd<-subDfdat[cc,]

pnt_gsh1<-apply(subd[,16:17],1,mean)
pnt_gush2<-apply(subd[,12:13],1,mean)
pnt_gsh1[subd[,9]=="+-"]<-sz-pnt_gsh1[subd[,9]=="+-"]
plot(pnt_gsh1,pnt_gush2,pch=19,col="darkgray")

## inverted translocation
aa<-which(pnt_gsh1 > 2.4e7 & (pnt_gsh1-pnt_gush2) < -.5e7 & pnt_gsh1 < 3.6e7 & (pnt_gsh1-pnt_gush2) > -1.3e7 )
points(pnt_gsh1[aa],pnt_gsh1[aa]-pnt_gush2[aa],col="gray",pch=19)

min(subd[aa,12])
#[1] 33112541
max(subd[aa,13])
#[1] 44121870
b1<-subd[aa,16]
b2<-subd[aa,17]
b1[subd[aa,9]=="+-"]<-sz-b1[subd[aa,9]=="+-"]
b2[subd[aa,9]=="+-"]<-sz-b2[subd[aa,9]=="+-"]
min(c(b1,b2))
#[1] 24457103
max(c(b1,b2))
#[1] 32202650

## rest of translocation
aa<-which(pnt_gsh1 > 3.0e7 & (pnt_gsh1-pnt_gush2) > .1e7 & pnt_gsh1 < 4e7 & (pnt_gsh1-pnt_gush2) < 2e7 )
points(pnt_gsh1[aa],pnt_gsh1[aa]-pnt_gush2[aa],col="gray",pch=19)
min(subd[aa,12])
#[1] 24803527
max(subd[aa,13])
#[1] 32951302
b1<-subd[aa,16]
b2<-subd[aa,17]
b1[subd[aa,9]=="+-"]<-sz-b1[subd[aa,9]=="+-"]
b2[subd[aa,9]=="+-"]<-sz-b2[subd[aa,9]=="+-"]
min(c(b1,b2))
#[1] 32202951
max(c(b1,b2))
#[1] 39030359

########### the rest for completeness ############

## read synteny dat
#########################################################################
## gsh2 x gush1
##############
dat_gsh2_gush1<-fread("out_cactusStripe_TcrGSH2_TcrGUSH1.psl",header=FALSE)
dfdat<-as.data.frame(dat_gsh2_gush1)

## target = GSH2
## query = GUSH1

## verify/identify large scaffolds (all should be) 
xx<-table(dfdat[,14])
gsh2Ch<-names(xx)[xx>500]
xx<-table(dfdat[,10])
gshxx<-table(dfdat[,10])
gush1Ch<-names(xx)[xx>500]
keep<-(dfdat[,14] %in% gsh2Ch) & (dfdat[,10] %in% gush1Ch)
subDfdat<-dfdat[keep,]## retains all, 13 CH each

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tc_gsh2<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(2,103,8)])
tc_gush1<-as.numeric(unlist(strsplit(x=rownames(tab),split="_"))[seq(2,52,4)])
sizes<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(8,104,8)])## sizes for gs h2


## normalize  
ntab<-tab
for(i in 1:13){
	ntab[i,]<-ntab[i,]/sum(ntab[i,])
}
pdf("SynTcrisStripe_GUSH1_GSH2.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. cristinae (H GUS1)",ylab="T. cristinae (H GS2)",cex.lab=1.4)
axis(1,at=seq(0,13,length.out=13)/13,tc_gush1,las=2)
axis(2,at=seq(0,13,length.out=13)/13,tc_gsh2,las=2)
box()
dev.off()

## colinearity plots for all homologous chromsomes
## chrom number, gs, gsr1, gsh2, gush1
chtab<-matrix(c(1,8483,12,13,22,
	2,14640,6,6,23,
	3,42935,2,2,16,
	4,42912,1,1,64,
	5,18722,7,12,5,
	6,9928,8,5,11,
	7,10660,10,8,54,
	8,7748,11,4,7,
	9,16151,5,9,46,
	10,14160,4,7,15,
	11,12033,9,10,2,
	12,12380,13,11,1,
	13,14101,3,3,36),nrow=13,ncol=5,byrow=TRUE)


pdf("AlnPlotTcris_GUSH1_GSH2.pdf",width=11,height=11)
par(mfrow=c(4,4))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:13){
	tcr_gush1<-grep(x=subDfdat[,10],pattern=paste("d6S_",chtab[i,5],"_",sep="")) 
	tcr_gsh2<-grep(x=subDfdat[,14],pattern=paste("old_",chtab[i,4],"_",sep=""))
	cc<-tcr_gsh2[tcr_gsh2 %in% tcr_gush1]
	subd<-subDfdat[cc,]
	xub<-max(subd[,13]);yub<-max(subd[,17])	

	plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='l',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="H GUS1",ylab="H GS2")
	title(main=paste("Chromosome",i),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 2:N){
		if(subd[j,9]=="++"){
			lines(subd[j,12:13],subd[j,16:17])
		}
		else{
			lines(subd[j,12:13],sizes[which(tc_gsh2==chtab[i,4])]-subd[j,16:17])
		}
	}
}
dev.off()

## read synteny dat
#########################################################################
## gsh2 x gush2
##############
dat_gsh2_gush2<-fread("out_cactusStripe_TcrGSH2_TcrGUSH2.psl",header=FALSE)
dfdat<-as.data.frame(dat_gsh2_gush2)

## target = GSH2
## query = GUSH2

## verify/identify large scaffolds (all should be) 
xx<-table(dfdat[,14])
gsh2Ch<-names(xx)[xx>500]
xx<-table(dfdat[,10])
gshxx<-table(dfdat[,10])
gush2Ch<-names(xx)[xx>500]
keep<-(dfdat[,14] %in% gsh2Ch) & (dfdat[,10] %in% gush2Ch)
subDfdat<-dfdat[keep,]## retains all, 13 CH each

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tc_gsh2<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(2,103,8)])
tc_gush2<-as.numeric(unlist(strsplit(x=rownames(tab),split="_"))[seq(2,52,4)])
sizes<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(8,104,8)])## sizes for gs h2


## normalize  
ntab<-tab
for(i in 1:13){
	ntab[i,]<-ntab[i,]/sum(ntab[i,])
}
pdf("SynTcrisStripe_GUSH2_GSH2.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. cristinae (H GUS2)",ylab="T. cristinae (H GS2)",cex.lab=1.4)
axis(1,at=seq(0,13,length.out=13)/13,tc_gush2,las=2)
axis(2,at=seq(0,13,length.out=13)/13,tc_gsh2,las=2)
box()
dev.off()

## colinearity plots for all homologous chromsomes
## chrom number, gs, gsr1, gsh2, gush2
chtab<-matrix(c(1,8483,12,13,15,
	2,14640,6,6,1,
	3,42935,2,2,3,
	4,42912,1,1,35,
	5,18722,7,12,10,
	6,9928,8,5,44,
	7,10660,10,8,7,
	8,7748,11,4,23,
	9,16151,5,9,21,
	10,14160,4,7,16,
	11,12033,9,10,12,
	12,12380,13,11,36,
	13,14101,3,3,8),nrow=13,ncol=5,byrow=TRUE)


pdf("AlnPlotTcris_GUSH2_GSH2.pdf",width=11,height=11)
par(mfrow=c(4,4))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:13){
	tcr_gush2<-grep(x=subDfdat[,10],pattern=paste("45T_",chtab[i,5],"_",sep="")) 
	tcr_gsh2<-grep(x=subDfdat[,14],pattern=paste("old_",chtab[i,4],"_",sep=""))
	cc<-tcr_gsh2[tcr_gsh2 %in% tcr_gush2]
	subd<-subDfdat[cc,]
	xub<-max(subd[,13]);yub<-max(subd[,17])	

	plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='l',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="H GUS2",ylab="H GS2")
	title(main=paste("Chromosome",i),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 2:N){
		if(subd[j,9]=="++"){
			lines(subd[j,12:13],subd[j,16:17])
		}
		else{
			lines(subd[j,12:13],sizes[which(tc_gsh2==chtab[i,4])]-subd[j,16:17])
		}
	}
}
dev.off()
