library(data.table)

## read synteny dat

#################################################################################
## start with refugio GS 1 vs main mountain GS
dat_gsr1_gs<-fread("PslFiles/out_cactusTcrGSR1_TcrGS.psl",header=FALSE)
dfdat<-as.data.frame(dat_gsr1_gs)

## target = GS
## query = GSR1

## identify large scaffolds for GSR!
xx<-table(dfdat[,14])
gsr1Ch<-names(xx)[xx>500]
xx<-table(dfdat[,10])
gsCh<-names(xx)[xx>500]
keep<-(dfdat[,14] %in% gsr1Ch) & (dfdat[,10] %in% gsCh)
subDfdat<-dfdat[keep,]

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tc_gs<-as.numeric(unlist(strsplit(x=rownames(tab),split="_"))[seq(2,50,4)])
tc_gsr1<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(2,103,8)])
sizes<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(8,104,8)])


## normalize with respect to main mountain green stripe 
ntab<-tab
for(i in 1:13){
	ntab[i,]<-ntab[i,]/sum(ntab[i,])
}
ch<-c(7,11,12,13,10,2,9,5,4,3,8,1,6)
## chrom vs scaf
cbind(tc_gs,ch)
#      tc_sc ch
# [1,] 10660  7
# [2,] 12033 11
# [3,] 12380 12
# [4,] 14101 13
# [5,] 14160 10
# [6,] 14640  2
# [7,] 16151  9
# [8,] 18722  5
# [9,] 42912  4
#[10,] 42935  3
#[11,]  7748  8
#[12,]  8483  1
#[13,]  9928  6
pdf("SynTcrisStripe_GSR1_GS.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. cristinae (GS)",ylab="T. cristinae (R GS1)",cex.lab=1.4,breaks=seq(0,1,length.out=13))
axis(2,at=seq(0,13,length.out=13)/13,tc_gsr1,las=2)
axis(1,at=seq(0,13,length.out=13)/13,ch,las=2)
box()
dev.off()



## colinearity plots for all homologous chromsomes
chtab<-matrix(c(8483,12,
	14640,6,
	42935,2,
	42912,1,
	18722,7,
	9928,8,
	10660,10,
	7748,11,
	16151,5,
	14160,4,
	12033,9,
	12380,13,
	14101,3),nrow=13,ncol=2,byrow=TRUE)




pdf("AlnPlotTcris_GSR1_GS.pdf",width=10,height=10)
par(mfrow=c(3,3))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:13){
	tcr_gs<-grep(x=subDfdat[,10],pattern=paste("_",chtab[i,1],"_",sep="")) 
	tcr_gsr1<-grep(x=subDfdat[,14],pattern=paste("old_",chtab[i,2],"_",sep=""))
	cc<-tcr_gs[tcr_gs %in% tcr_gsr1]
	subd<-subDfdat[cc,]
	xub<-max(subd[,13]);yub<-max(subd[,17])	

	plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='l',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="T. cristinae (GS main)",ylab="T. cristinae (R GS1)")
	title(main=paste("Chrom.",i),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 2:N){
		if(subd[j,9]=="++"){
			lines(subd[j,12:13],subd[j,16:17])
		}
		else{
			lines(subd[j,12:13],sizes[which(tc_gsr1==chtab[i,2])]-subd[j,16:17],col="cadetblue")
		}
	}
}
dev.off()

#################################################################################
#################################################################################
## now refugio GS1 vs GS2
dat_gsr1_gsr2<-fread(ff[[2]],header=FALSE)
dfdat<-as.data.frame(dat_gsr1_gsr2)

## target = GSR2
## query = GSR1

## identify large scaffolds for GSR 1 and 2
xx<-table(dfdat[,14])
gsr1Ch<-names(xx)[xx>500]
xx<-table(dfdat[,10])
gsr2Ch<-names(xx)[xx>500]
keep<-(dfdat[,14] %in% gsr1Ch) & (dfdat[,10] %in% gsr2Ch)
subDfdat<-dfdat[keep,]

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tc_gsr2<-as.numeric(unlist(strsplit(x=rownames(tab),split="_"))[seq(2,95,8)])
tc_gsr1<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(2,103,8)])
sizes<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(8,104,8)])


## normalize  
ntab<-tab
for(i in 1:12){
	ntab[i,]<-ntab[i,]/sum(ntab[i,])
}
pdf("SynTcrisStripe_GSR1_GSR2.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. cristinae (R GS2)",ylab="T. cristinae (R GS1)",cex.lab=1.4,breaks=seq(0,1,length.out=13))
axis(2,at=seq(0,13,length.out=13)/13,tc_gsr1,las=2)
axis(1,at=seq(0,12,length.out=12)/12,tc_gsr2,las=2)
box()
dev.off()



## colinearity plots for all homologous chromsomes
## chrom number, gs,gsr1, gsr2
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




pdf("AlnPlotTcris_GSR1_GSR2.pdf",width=10,height=10)
par(mfrow=c(3,3))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:13){
	tcr_gsr2<-grep(x=subDfdat[,10],pattern=paste("old_",chtab[i,4],"_",sep="")) 
	tcr_gsr1<-grep(x=subDfdat[,14],pattern=paste("old_",chtab[i,3],"_",sep=""))
	cc<-tcr_gsr2[tcr_gsr2 %in% tcr_gsr1]
	subd<-subDfdat[cc,]
	xub<-max(subd[,13]);yub<-max(subd[,17])	

	plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='l',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="T. cristinae (R GS2)",ylab="T. cristinae (R GS1)")
	title(main=paste("Chrom.",i),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 2:N){
		if(subd[j,9]=="++"){
			lines(subd[j,12:13],subd[j,16:17])
		}
		else{
			lines(subd[j,12:13],sizes[which(tc_gsr1==chtab[i,3])]-subd[j,16:17],col="cadetblue")
		}
	}
}
dev.off()

#################################################################################
#################################################################################
## now refugio GS1 vs GUS2
dat_gsr1_gusr1<-fread("PslFiles/out_cactusTcrGSR1_TcrGUSR1.psl",header=FALSE)

dfdat<-as.data.frame(dat_gsr1_gusr1)

## target = GUSR1
## query = GSR1

## identify large scaffolds for GSR 1 and GUSR1
xx<-table(dfdat[,14])
gsr1Ch<-names(xx)[xx>500]
xx<-table(dfdat[,10])
gusr1Ch<-names(xx)[xx>500]
keep<-(dfdat[,14] %in% gsr1Ch) & (dfdat[,10] %in% gusr1Ch)
subDfdat<-dfdat[keep,]

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tc_gusr1<-as.numeric(unlist(strsplit(x=rownames(tab),split="_"))[seq(2,95,8)])
tc_gsr1<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(2,103,8)])
sizes<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(8,104,8)])


## normalize  
ntab<-tab
for(i in 1:12){
	ntab[i,]<-ntab[i,]/sum(ntab[i,])
}
pdf("SynTcrisStripe_GSR1_GUSR1.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. cristinae (R GUS1)",ylab="T. cristinae (R GS1)",cex.lab=1.4,breaks=seq(0,1,length.out=13))
axis(2,at=seq(0,13,length.out=13)/13,tc_gsr1,las=2)
axis(1,at=seq(0,12,length.out=12)/12,tc_gusr1,las=2)
box()
dev.off()



## colinearity plots for all homologous chromsomes
## chrom number, gs,gsr1, gusr1
chtab<-matrix(c(1,8483,12,4,
	2,14640,6,11,
	3,42935,2,1,
	4,42912,1,1,
	5,18722,7,6,
	6,9928,8,7,
	7,10660,10,10,
	8,7748,11,3,
	9,16151,5,9,
	10,14160,4,8,
	11,12033,9,5,
	12,12380,13,12,
	13,14101,3,2),nrow=13,ncol=4,byrow=TRUE)




pdf("AlnPlotTcris_GSR1_GUSR1.pdf",width=10,height=10)
par(mfrow=c(3,3))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:13){
	tcr_gusr1<-grep(x=subDfdat[,10],pattern=paste("old_",chtab[i,4],"_",sep="")) 
	tcr_gsr1<-grep(x=subDfdat[,14],pattern=paste("old_",chtab[i,3],"_",sep=""))
	cc<-tcr_gusr1[tcr_gusr1 %in% tcr_gsr1]
	subd<-subDfdat[cc,]
	xub<-max(subd[,13]);yub<-max(subd[,17])	

	plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='l',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="T. cristinae (R GUS1)",ylab="T. cristinae (R GS1)")
	title(main=paste("Chrom.",i),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 2:N){
		if(subd[j,9]=="++"){
			lines(subd[j,12:13],subd[j,16:17])
		}
		else{
			lines(subd[j,12:13],sizes[which(tc_gsr1==chtab[i,3])]-subd[j,16:17],col="cadetblue")
		}
	}
}
dev.off()

#################################################################################
#################################################################################
## now refugio GS1 vs GUS2
dat_gsr1_gusr2<-fread("PslFiles/out_cactusTcrGSR1_TcrGUSR2.psl",,header=FALSE)
dfdat<-as.data.frame(dat_gsr1_gusr2)

## target = GUSR2
## query = GSR1

## identify large scaffolds for GSR 1 and GUSR2
xx<-table(dfdat[,14])
gsr1Ch<-names(xx)[xx>500]
xx<-table(dfdat[,10])
gusr2Ch<-names(xx)[xx>500]
keep<-(dfdat[,14] %in% gsr1Ch) & (dfdat[,10] %in% gusr2Ch)
subDfdat<-dfdat[keep,]

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tc_gusr2<-as.numeric(unlist(strsplit(x=rownames(tab),split="_"))[seq(2,95,8)])
tc_gsr1<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(2,103,8)])
sizes<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(8,104,8)])


## normalize  
ntab<-tab
for(i in 1:12){
	ntab[i,]<-ntab[i,]/sum(ntab[i,])
}
pdf("SynTcrisStripe_GSR1_GUSR2.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. cristinae (R GUS2)",ylab="T. cristinae (R GS1)",cex.lab=1.4,breaks=seq(0,1,length.out=13))
axis(2,at=seq(0,13,length.out=13)/13,tc_gsr1,las=2)
axis(1,at=seq(0,12,length.out=12)/12,tc_gusr1,las=2)
box()
dev.off()



## colinearity plots for all homologous chromsomes
## chrom number, gs,gsr1, gusr2
chtab<-matrix(c(1,8483,12,6,
	2,14640,6,9,
	3,42935,2,1,
	4,42912,1,1,
	5,18722,7,4,
	6,9928,8,5,
	7,10660,10,11,
	8,7748,11,3,
	9,16151,5,8,
	10,14160,4,10,
	11,12033,9,7,
	12,12380,13,12,
	13,14101,3,2),nrow=13,ncol=4,byrow=TRUE)




pdf("AlnPlotTcris_GSR1_GUSR2.pdf",width=10,height=10)
par(mfrow=c(3,3))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:13){
	tcr_gusr2<-grep(x=subDfdat[,10],pattern=paste("old_",chtab[i,4],"_",sep="")) 
	tcr_gsr1<-grep(x=subDfdat[,14],pattern=paste("old_",chtab[i,3],"_",sep=""))
	cc<-tcr_gusr2[tcr_gusr2 %in% tcr_gsr1]
	subd<-subDfdat[cc,]
	xub<-max(subd[,13]);yub<-max(subd[,17])	

	plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='l',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="T. cristinae (R GUS2)",ylab="T. cristinae (R GS1)")
	title(main=paste("Chrom.",i),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 2:N){
		if(subd[j,9]=="++"){
			lines(subd[j,12:13],subd[j,16:17])
		}
		else{
			lines(subd[j,12:13],sizes[which(tc_gsr1==chtab[i,3])]-subd[j,16:17],col="cadetblue")
		}
	}
}
dev.off()

#################################################################################
#################################################################################
## now refugio GSR2 vs GUSR1
dat_gsr2_gusr1<-fread("PslFiles/out_cactusTcrGSR2_TcrGUSR1.psl",header=FALSE)

dfdat<-as.data.frame(dat_gsr2_gusr1)

## target = GUSR1
## query = GSR2

## identify large scaffolds for GSR 1 and 2
xx<-table(dfdat[,14])
gsr2Ch<-names(xx)[xx>500]
xx<-table(dfdat[,10])
gusr1Ch<-names(xx)[xx>500]
keep<-(dfdat[,14] %in% gsr2Ch) & (dfdat[,10] %in% gusr1Ch)
subDfdat<-dfdat[keep,]

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tc_gusr1<-as.numeric(unlist(strsplit(x=rownames(tab),split="_"))[seq(2,95,8)])
tc_gsr2<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(2,95,8)])
sizes<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(8,96,8)])


## normalize  
ntab<-tab
for(i in 1:12){
	ntab[i,]<-ntab[i,]/sum(ntab[i,])
}
pdf("SynTcrisStripe_GSR2_GUSR1.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. cristinae (R GUS1)",ylab="T. cristinae (R GS2)",cex.lab=1.4,breaks=seq(0,1,length.out=13))
axis(2,at=seq(0,12,length.out=12)/12,tc_gsr2,las=2)
axis(1,at=seq(0,12,length.out=12)/12,tc_gusr1,las=2)
box()
dev.off()



## colinearity plots for all homologous chromsomes
## chrom number, gs,gsr2, gusr1
chtab<-matrix(c(1,8483,10,4,
	2,14640,8,11,
	3,42935,1,1,
	5,18722,2,6,
	6,9928,5,7,
	7,10660,7,10,
	8,7748,9,3,
	9,16151,4,9,
	10,14160,3,8,
	11,12033,6,5,
	12,12380,12,12,
	13,14101,11,2),nrow=12,ncol=4,byrow=TRUE)




pdf("AlnPlotTcris_GSR2_GUSR1.pdf",width=10,height=10)
par(mfrow=c(3,3))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:12){
	tcr_gusr1<-grep(x=subDfdat[,10],pattern=paste("old_",chtab[i,4],"_",sep="")) 
	tcr_gsr2<-grep(x=subDfdat[,14],pattern=paste("old_",chtab[i,3],"_",sep=""))
	cc<-tcr_gusr1[tcr_gusr1 %in% tcr_gsr2]
	subd<-subDfdat[cc,]
	xub<-max(subd[,13]);yub<-max(subd[,17])	

	plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='l',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="T. cristinae (R GUS1)",ylab="T. cristinae (R GS2)")
	title(main=paste("Chrom.",chtab[i,1]),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 2:N){
		if(subd[j,9]=="++"){
			lines(subd[j,12:13],subd[j,16:17])
		}
		else{
			lines(subd[j,12:13],sizes[which(tc_gsr2==chtab[i,3])]-subd[j,16:17],col="cadetblue")
		}
	}
}
dev.off()

save(list=ls(),file="StripeAln.rdat")

#################################################################################
#################################################################################
## now refugio GSR2 vs GUSR2
dat_gsr2_gusr2<-fread("PslFiles/out_cactusTcrGSR2_TcrGUSR2.psl",header=FALSE)
dfdat<-as.data.frame(dat_gsr2_gusr2)

## target = GUSR2
## query = GSR2

## identify large scaffolds for GSR 2 and GUSR 2
xx<-table(dfdat[,14])
gsr2Ch<-names(xx)[xx>500]
xx<-table(dfdat[,10])
gusr2Ch<-names(xx)[xx>500]
keep<-(dfdat[,14] %in% gsr2Ch) & (dfdat[,10] %in% gusr2Ch)
subDfdat<-dfdat[keep,]

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tc_gusr2<-as.numeric(unlist(strsplit(x=rownames(tab),split="_"))[seq(2,95,8)])
tc_gsr2<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(2,95,8)])
sizes<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(8,96,8)])


## normalize  
ntab<-tab
for(i in 1:12){
	ntab[i,]<-ntab[i,]/sum(ntab[i,])
}
pdf("SynTcrisStripe_GSR2_GUSR2.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. cristinae (R GUS2)",ylab="T. cristinae (R GS2)",cex.lab=1.4,breaks=seq(0,1,length.out=13))
axis(2,at=seq(0,12,length.out=12)/12,tc_gsr2,las=2)
axis(1,at=seq(0,12,length.out=12)/12,tc_gusr2,las=2)
box()
dev.off()



## colinearity plots for all homologous chromsomes
## chrom number, gs,gsr2, gusr2
chtab<-matrix(c(1,8483,10,6,
	2,14640,8,9,
	3,42935,1,1,
	5,18722,2,4,
	6,9928,5,5,
	7,10660,7,11,
	8,7748,9,3,
	9,16151,4,8,
	10,14160,3,10,
	11,12033,6,7,
	12,12380,12,12,
	13,14101,11,2),nrow=12,ncol=4,byrow=TRUE)




pdf("AlnPlotTcris_GSR2_GUSR2.pdf",width=10,height=10)
par(mfrow=c(3,3))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:12){
	tcr_gusr2<-grep(x=subDfdat[,10],pattern=paste("old_",chtab[i,4],"_",sep="")) 
	tcr_gsr2<-grep(x=subDfdat[,14],pattern=paste("old_",chtab[i,3],"_",sep=""))
	cc<-tcr_gusr2[tcr_gusr2 %in% tcr_gsr2]
	subd<-subDfdat[cc,]
	xub<-max(subd[,13]);yub<-max(subd[,17])	

	plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='l',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="T. cristinae (R GUS2)",ylab="T. cristinae (R GS2)")
	title(main=paste("Chrom.",chtab[i,1]),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 2:N){
		if(subd[j,9]=="++"){
			lines(subd[j,12:13],subd[j,16:17])
		}
		else{
			lines(subd[j,12:13],sizes[which(tc_gsr2==chtab[i,3])]-subd[j,16:17],col="cadetblue")
		}
	}
}
dev.off()

#################################################################################
#################################################################################
## now refugio GUSR1 vs GUS main
dat_gusr1_gus<-fread("PslFiles/out_cactusTcrGUSR1_TcrGUS.psl",header=FALSE)

dfdat<-as.data.frame(dat_gusr1_gus)

## target = GUS main
## query = GUSR1

## identify large scaffolds for GUSR 1 and GUS main
xx<-table(dfdat[,14])
gusrCh<-names(xx)[xx>500]
xx<-table(dfdat[,10])
gusCh<-names(xx)[xx>500]
keep<-(dfdat[,14] %in% gusr1Ch) & (dfdat[,10] %in% gusCh)
subDfdat<-dfdat[keep,]

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tc_gus<-as.numeric(unlist(strsplit(x=rownames(tab),split="_"))[seq(2,40,3)])
tc_gusr1<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(2,95,8)])

#ss<-read.table("../t_crist_gus/scafLen.txt",header=FALSE)
#sizes<-rep(NA,13)
#for(i in 1:13){
#	sizes[i]<-ss[grep(pattern=paste("_",tc_gus[i],"_",sep=""),x=ss[,1]),2]
#}
sizes<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(8,96,8)])




## normalize  
ntab<-tab
for(i in 1:12){
	ntab[,i]<-ntab[,i]/sum(ntab[,i])
}
pdf("SynTcrisStripe_GUSR1_GUS.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. cristinae (GUS main)",ylab="T. cristinae (R GUS1)",cex.lab=1.4,breaks=seq(0,1,length.out=13))
axis(2,at=seq(0,12,length.out=12)/12,tc_gusr1,las=2)
axis(1,at=seq(0,13,length.out=13)/13,tc_gus,las=2)
box()
dev.off()



## colinearity plots for all homologous chromsomes
## chrom number, gs, gus, gusr1
chtab<-matrix(c(1,8483,11166,4,
	2,14640,2255,11,
	3,42935,11174,1,
	4,42912,12,1,
	5,18722,1,6,
	6,9928,935,7,
	7,10660,100,10,
	8,7748,3,3,
	9,16151,1001,9,
	10,14160,24,8,
	11,12033,3508,5,
	12,12380,11137,12,
	13,14101,9884,2),nrow=13,ncol=4,byrow=TRUE)




pdf("AlnPlotTcris_GUSR1_GUS.pdf",width=10,height=10)
par(mfrow=c(3,3))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:13){
	tcr_gus<-grep(x=subDfdat[,10],pattern=paste("_",chtab[i,3],"_",sep="")) 
	tcr_gusr1<-grep(x=subDfdat[,14],pattern=paste("old_",chtab[i,4],"_",sep=""))
	cc<-tcr_gus[tcr_gus %in% tcr_gusr1]
	subd<-subDfdat[cc,]
	xub<-max(subd[,13]);yub<-max(subd[,17])	

	plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='l',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="T. cristinae (GUS main)",ylab="T. cristinae (R GUS1)")
	title(main=paste("Chrom.",chtab[i,1]),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 2:N){
		if(subd[j,9]=="++"){
			lines(subd[j,12:13],subd[j,16:17])
		}
		else{
			lines(subd[j,12:13],sizes[which(tc_gusr1==chtab[i,4])]-subd[j,16:17],col="cadetblue")
		}
	}
}
dev.off()

#################################################################################
#################################################################################
## now refugio GUSR1 vs GUSR2 main
dat_gusr1_gusr2<-fread("PslFiles/out_cactusTcrGUSR1_TcrGUSR2.psl",header=FALSE)
dfdat<-as.data.frame(dat_gusr1_gusr2)

## target = GUSR2
## query = GUSR1

## identify large scaffolds for GUSR 1 and GUS main
xx<-table(dfdat[,14])
gusr1Ch<-names(xx)[xx>500]
xx<-table(dfdat[,10])
gusr2Ch<-names(xx)[xx>500]
keep<-(dfdat[,14] %in% gusr1Ch) & (dfdat[,10] %in% gusr2Ch)
subDfdat<-dfdat[keep,]

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tc_gusr1<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(2,95,8)])
tc_gusr2<-as.numeric(unlist(strsplit(x=rownames(tab),split="_"))[seq(2,95,8)])

sizes<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(8,96,8)])




## normalize  
ntab<-tab
for(i in 1:12){
	ntab[i,]<-ntab[i,]/sum(ntab[i,])
}
pdf("SynTcrisStripe_GUSR1_GUSR2.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. cristinae (R GUS2)",ylab="T. cristinae (R GUS1)",cex.lab=1.4,breaks=seq(0,1,length.out=13))
axis(2,at=seq(0,12,length.out=12)/12,tc_gusr1,las=2)
axis(1,at=seq(0,12,length.out=12)/12,tc_gusr2,las=2)
box()
dev.off()

pdf("SynTcrisStripe_GUSR1_GUSR2_legend.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image.plot(ntab,axes=FALSE,col = hcl.colors(12, "YlOrRd", rev = TRUE),xlab="T. cristinae (R GUS2)",ylab="T. cristinae (R GUS1)",cex.lab=1.4,breaks=seq(0,1,length.out=13))
axis(2,at=seq(0,12,length.out=12)/12,tc_gusr1,las=2)
axis(1,at=seq(0,12,length.out=12)/12,tc_gusr2,las=2)
box()
dev.off()


## colinearity plots for all homologous chromsomes
## chrom number, gs, gusr2, gusr1
chtab<-matrix(c(1,8483,6,4,
	2,14640,9,11,
	3,42935,1,1,
	5,18722,4,6,
	6,9928,5,7,
	7,10660,11,10,
	8,7748,3,3,
	9,16151,8,9,
	10,14160,10,8,
	11,12033,7,5,
	12,12380,12,12,
	13,14101,2,2),nrow=12,ncol=4,byrow=TRUE)


pdf("AlnPlotTcris_GUSR1_GUSR2.pdf",width=10,height=10)
par(mfrow=c(3,3))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:12){
	tcr_gusr2<-grep(x=subDfdat[,10],pattern=paste("old_",chtab[i,3],"_",sep="")) 
	tcr_gusr1<-grep(x=subDfdat[,14],pattern=paste("old_",chtab[i,4],"_",sep=""))
	cc<-tcr_gusr2[tcr_gusr2 %in% tcr_gusr1]
	subd<-subDfdat[cc,]
	xub<-max(subd[,13]);yub<-max(subd[,17])	

	plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='l',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="T. cristinae (R GUS2)",ylab="T. cristinae (R GUS1)")
	title(main=paste("Chrom.",chtab[i,1]),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 2:N){
		if(subd[j,9]=="++"){
			lines(subd[j,12:13],subd[j,16:17])
		}
		else{
			lines(subd[j,12:13],sizes[which(tc_gusr1==chtab[i,4])]-subd[j,16:17],col="cadetblue")
		}
	}
}
dev.off()

save(list=ls(),file="StripeAln.rdat")

