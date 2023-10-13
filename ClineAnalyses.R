## estimating clines, selection on SV
library(rstan)
library(data.table)
library(scales)
library(codep)


## definition of locus
## gs, translocated = 56898118 to 65729704
## rest = 22442098 to 56898118
sv_gs<-c(22442098,65729704)

## gus, translocated = 22220178 to 31457192
## rest = 31523304 to 65829835

## read location data
ldat<-read.table("/uufs/chpc.utah.edu/common/home/gompert-group2/projects/tcr_nfds/Tcris_master_32.csv",header=TRUE,sep=",")

ref<-grep(pattern="^R",x=ldat$location)

lat<-tapply(X=ldat$latitude[ref],INDEX=ldat$location[ref],mean,na.rm=TRUE)
long<-tapply(X=ldat$longitude[ref],INDEX=ldat$location[ref],mean,na.rm=TRUE)

##################### striped ###############################

Ggs<-as.matrix(fread("g_tcr_refugio_gs.txt",sep=",",header=FALSE))
dim(Ggs)
#[1]   238 85558

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

snps<-read.table("../gwa_refugio/snps_gs.txt",header=FALSE)
sv_snps_gs<-which(snps[,1]==11 & snps[,2] >= sv_gs[1] & snps[,2] <= sv_gs[2])



aa<-which(names(lat) %in% unique(ph$loc))
plot(long[aa],lat[aa],pch=19,col=alpha("gray",.8))
text(long[aa],lat[aa],names(lat)[aa])

tapply(X=ph$host=="A",INDEX=ph$loc,mean)
#      R10       R11       R12       R13       R14       R17       R18       R19
#0.0000000 0.0000000 0.2580645 0.0000000 0.0000000 0.9824561 0.9000000 1.0000000
#      R23        R4        R5        R6        R7        R8        R9
#1.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
## basically all A until R 17, then mostly C, R12 slight outlier

## ordinate geo to get 1 D transect
locs<-cbind(lat[aa],long[aa])
## change to great circle geo distance
D<-as.matrix(gcd.vife(locs)) ## distances in km
pc<-prcomp(D,center=TRUE,scale=FALSE)
summary(pc)
#Importance of components:
#                            PC1      PC2      PC3      PC4       PC5       PC6
#Standard deviation     0.7259 0.4723 0.21899 0.09997 0.05794 0.03631 0.03581
#Proportion of Variance 0.6453 0.2733 0.05874 0.01224 0.00411 0.00161 0.00157
#Cumulative Proportion  0.6453 0.9186 0.97730 0.98954 0.99365 0.99526 0.99683
Dpc<-as.matrix(dist(pc$x[,1],method="euclidean",upper=TRUE,diag=TRUE))
Dpcv<-as.vector(Dpc[upper.tri(Dpc)])
Dv<-as.vector(D[upper.tri(D)])

## get relationship, will need this to conver scales between PC and original
o<-lm(Dv ~ Dpcv)
summary(o)
#            Estimate Std. Error t value Pr(>|t|)
#(Intercept)  0.10316    0.01360   7.587 1.51e-11 ***

#Dpcv         0.36624    0.01324  27.651  < 2e-16 ***
#Residual standard error: 0.0808 on 103 degrees of freedom
#Multiple R-squared:  0.8813,	Adjusted R-squared:  0.8801
#F-statistic: 764.6 on 1 and 103 DF,  p-value: < 2.2e-16
o$coefficients
#(Intercept)        Dpcv
#  0.1031575   0.3662368

## now genetics get the cline
pcg<-prcomp(Ggs[,sv_snps_gs],center=TRUE,scale=FALSE)
ko<-kmeans(x=pcg$x[,1:2],centers=6,iter.max=100,nstart=25)
ko$cluster
#  [1] 5 5 4 2 5 4 2 2 4 2 5 6 4 4 2 4 1 2 4 2 5 6 4 4 2 2 5 1 4 4 2 2 2 5 6 4 6
# [38] 1 4 4 4 2 5 5 6 6 6 6 3 2 2 2 6 5 5 5 2 6 6 1 6 4 4 2 2 2 1 4 4 2 5 2 2 2
# [75] 6 6 6 2 1 2 4 1 4 2 2 1 3 1 4 6 4 5 2 1 5 3 3 1 6 2 4 5 2 2 2 2 2 4 4 5 2
#[112] 4 5 2 2 4 2 4 2 4 2 2 1 5 5 3 3 3 3 3 3 3 6 6 3 3 5 6 3 6 3 6 3 3 1 6 3 6
#[149] 5 6 1 3 3 3 3 3 6 3 3 3 3 1 6 3 3 6 3 3 3 6 3 3 3 6 6 3 6 6 6 6 3 6 3 6 6
#[186] 6 3 6 6 3 3 6 6 3 6 5 6 2 2 5 4 3 6 1 4 1 2 2 6 2 3 6 3 3 3 6 6 6 3 6 6 3
#[223] 6 6 5 5 3 1 1 4 2 4 1 1 2 1 5 1
cs<-rep("cadetblue",238)
cs[ph$morph=="U"]<-"forestgreen"
cs[ph$morph=="M"]<-"gray10"
plot(pcg$x[,1],pcg$x[,2],pch=19,col=alpha(cs,.5))

## 5 = MM
## 4 = UU
## 3 = SS
## 6 = MS
## 2 = MU
## 1 = US

## use NA for melanic, not counting for now
svGen<-matrix(NA,nrow=238,ncol=2)
svGen[ko$cluster==4,]<-c(0,0)
svGen[ko$cluster==3,]<-c(1,1)
svGen[ko$cluster==6,]<-c(1,NA)
svGen[ko$cluster==2,]<-c(0,NA)
svGen[ko$cluster==1,]<-c(0,1)

## location allele freqs
y<-tapply(INDEX=ph$loc,X=apply(svGen,1,sum,na.rm=TRUE),sum)
n<-tapply(INDEX=ph$loc,X=apply(is.na(svGen)==FALSE,1,sum),sum)
p<-y/n

plot(pc$x[,1],p,pch=19)

## see https://onlinelibrary.wiley.com/doi/pdfdirect/10.1111/j.0014-3820.2006.tb01139.x
## also good for LD based estimates of dispersal
# y = center, w = wdith, x = position
cf<-function(x,y,w){
		
	p<-(1 + tanh(2* (x - y)/w))/2

	return(p)
}

cdat<-list(N=length(y),x=pc$x[,1],y=y,n=n)
cfit<-stan("cline.stan",data=cdat)

#Inference for Stan model: cline.
#4 chains, each with iter=2000; warmup=1000; thin=1; 
#post-warmup draws per chain=1000, total post-warmup draws=4000.

#        mean se_mean   sd   2.5%    25%    50%    75%  97.5% n_eff Rhat
#c       0.02    0.00 0.05  -0.08  -0.02   0.02   0.05   0.11  2983    1
#w       1.29    0.00 0.15   1.03   1.19   1.28   1.38   1.60  2510    1
#p[1]    0.39    0.00 0.04   0.32   0.37   0.39   0.42   0.47  3225    1
#p[2]    0.32    0.00 0.04   0.25   0.30   0.32   0.35   0.40  3179    1
#p[3]    0.34    0.00 0.04   0.27   0.32   0.34   0.37   0.42  3203    1
#p[4]    0.49    0.00 0.04   0.42   0.47   0.49   0.52   0.57  3128    1
#p[5]    0.63    0.00 0.04   0.55   0.60   0.63   0.65   0.70  2911    1
#p[6]    0.93    0.00 0.02   0.88   0.91   0.93   0.94   0.96  2564    1
#p[7]    0.95    0.00 0.02   0.91   0.94   0.95   0.96   0.98  2540    1
#p[8]    0.95    0.00 0.02   0.92   0.94   0.95   0.97   0.98  2531    1
#p[9]    0.97    0.00 0.01   0.94   0.96   0.97   0.98   0.99  2494    1
#p[10]   0.03    0.00 0.01   0.01   0.02   0.02   0.03   0.05  2248    1
#p[11]   0.03    0.00 0.01   0.01   0.03   0.03   0.04   0.07  2303    1
#p[12]   0.07    0.00 0.02   0.03   0.05   0.06   0.08   0.11  2453    1
#p[13]   0.25    0.00 0.04   0.19   0.23   0.25   0.28   0.33  3038    1
#p[14]   0.38    0.00 0.04   0.31   0.35   0.38   0.40   0.45  3223    1
#p[15]   0.46    0.00 0.04   0.39   0.44   0.46   0.49   0.54  3174    1
#lp__  -40.43    0.02 1.03 -43.31 -40.77 -40.12 -39.71 -39.45  1720    1


pdf("RefugioCline.pdf",width=7,height=5)
par(mar=c(5,5,1,1))
xx<-seq(-1.2,1.2,0.01)
plot(xx,cf(xx,0.02,1.29),type='l',xlab="Transect location",ylab="Stripe SV frequency")
points(pc$x[,1],p,pch=19,col=alpha("black",.7))
dev.off()

## see https://books.google.com/books?hl=en&lr=&id=aFJFkVKskYIC&oi=fnd&pg=PA13&ots=MFk-dnK9NI&sig=r7KfAnHvJyLgyUJsMeLs7jpp33c#v=onepage&q&f=false
## barton and gale, genetic analysis of hybrid zones
## width = 1.732 sigma/ sqrt(s)
## sd of total distance moved = sqrt(2) * sigma
disp<-scan("dispersal.txt")
#sig<-sd(disp)/sqrt(2)
#[1] 25.60094
sig<-mean(disp)/sqrt(2) ## this makes more sense, double check, but doesn't matter
#[1] 22.67093

## put width on original scale
w<-0.1031575 + 1.29 * 0.3662368

# w = 1.732 sig / sqrt(s)
# sqrt(s) = 1.732 * sig / w
# s = (1.732 * sig / w)^2
## divide by 1000 as sig in meters not km
(1.732 * (sig/1000)/w)^2
#[1] 0.004653593

# w = sqrt(8 or 3 sig^2)/s

save(list=ls(),file="cline.rdat")


## other clines
## define parents
popA<-c("R4","R5","R6")
popB<-c("R17","R18","R19","R23")
idA<-which(ph$loc %in% popA)
idB<-which(ph$loc %in% popB)
no8_snps_gs<-which(snps[,1]!=11)
## parent af, no 8
p_A<-apply(Ggs[idA,no8_snps_gs],2,mean)/2
p_B<-apply(Ggs[idB,no8_snps_gs],2,mean)/2
anc<-no8_snps_gs[which(abs(p_A-p_B) > .2)]
#table(snps[anc,1])
#
# 1  2  3  4  5  6  7  8  9 10 12 13
# 7  7 15  3  2 11  5  5  1  3  4  3

w_raw<-numeric(length(anc))
for(i in 1:length(anc)){
	a<-tapply(INDEX=ph$loc,X=round(Ggs[,anc[i]],0),sum)
	b<-tapply(INDEX=ph$loc,X=(Ggs[,anc[i]] >= 0) * 2,sum)
	ps<-a/b
	if(cor(pc$x[,1],ps) < 0){
		a<-b-a
		ps<-1-ps
	}
	sdat<-list(N=length(a),x=pc$x[,1],y=a,n=b)
	sfit<-stan("llcline.stan",data=sdat,warmup=3000,iter=5000)
	w_raw[i]<-abs(mean(extract(sfit,"w")[[1]]))
}

w_snps<-0.1031575 + w_raw * 0.3662368

hist(w_snps,xlim=c(0,7))
abline(v=w)

## ld
## location allele freqs
Gsv<-apply(svGen,1,sum,na.rm=TRUE)
Nsv<-apply(is.na(svGen)==FALSE,1,sum)
idH<-c(1:238)[-c(idA,idB)]
ldK<-which(Nsv == 2 & (1:238) %in% idH)

ld<-numeric(length(anc))
for(i in 1:length(anc)){
	ld[i]<-abs(cor(Gsv[ldK],round(Ggs[ldK,anc[i]],0)))
}

p_hyb_snps<-apply(round(Ggs[ldK,anc],0),2,mean)/2
p_hyb_sv<-mean(Gsv[ldK])/2

mnw<-mean(w_snps+w)/2
mnR<-mean(ld)
con<-1/sqrt(mean(p_hyb_sv * p_hyb_snps * (1-p_hyb_snps) * (1-p_hyb_sv)))
## 4.7
sig_gen<-mnw * sqrt(mnR*.5/con)

(1.732 * sig_gen/w)^2
## gives
##[1] 0.3780105


