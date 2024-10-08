## summarize TE information for the inverted translocation 

## Chr8
## h154
## gs1 = Scaffold_11__2_contigs__length_91821751
## gus2 = ScrX45T_23_HRSCAF_264
## refugio
## gs1 = Scaffold_11__3_contigs__length_97865747
## gus1 = Scaffold_3__1_contigs__length_98311894


## bounds of region
## h154
## gs1 (inverted bit than rest) = 24457103,32202650 and 32202951,39030359
## gus2 (inverted bit than rest) = 33112541,44121870 and 24803527,32951302
## refugio
## gs1 (inverted bit than rest) = 56898118,65729704 and 22442098,56898118
## gus1 (inverted bit than rest) = 22220178,31457192 and 31523304,65829835

gff_h154_gs_h1<-read.table("t_crist_gs_hap_cen4119/HiRise/Hap1/final_assembly.fasta.out.gff", header=FALSE)
gff_h154_gs_h2<-read.table("t_crist_gs_hap_cen4119/HiRise/Hap2/final_assembly.fasta.out.gff", header=FALSE)
gff_h154_gus_h1<-read.table("t_crist_gus_hap_cen4280/HiRise/Hap1/ojincantatabio-cen4280-hap1-mb-hirise-ig5ps__01-30-2024__hic_output.fasta.out.gff",header=FALSE)
gff_h154_gus_h2<-read.table("t_crist_gus_hap_cen4280/HiRise/Hap2/ojincantatabio-cen4280-hap2-mb-hirise-i2xb7__01-30-2024__hic_output.fasta.out.gff",header=FALSE)

gff_refug_gus_h1<-read.table("t_crist_refug_green/HiRise/hap1/ojincantatabio-cen4120-hap1-mb-hirise-wlbll__08-15-2023__final_assembly.fasta.out.gff", header=FALSE)
gff_refug_gus_h2<-read.table("t_crist_refug_green/HiRise/hap2/ojincantatabio-cen4120-hap2-mb-hirise-bn0ko__08-15-2023__final_assembly.fasta.out.gff", header=FALSE)
gff_refug_gs_h1<-read.table("t_crist_refug_stripe/HiRise/hap1/ojincantatabio-cen4122-hap1-mb-hirise-g4hzf__08-10-2023__final_assembly.fasta.out.gff", header=FALSE)
gff_refug_gs_h2<-read.table("t_crist_refug_stripe/HiRise/hap2/ojincantatabio-cen4122-hap2-mb-hirise-14fv0__08-10-2023__final_assembly.fasta.out.gff", header=FALSE)


## indexes for ch8, focus on four "core" genomes for now
ch8i_h154_gus_h2<-grep(pattern="ScrX45T_23;HRSCAF=264",x=gff_h154_gus_h2$V1,fixed=TRUE)
ch8i_h154_gs_h1<-grep(pattern="Scaffold_11__2_contigs__length_91821751",x=gff_h154_gs_h1$V1,fixed=TRUE)
ch8i_refug_gus_h1<-grep(pattern="Scaffold_3__1_contigs__length_98311894",x=gff_refug_gus_h1$V1,fixed=TRUE)
ch8i_refug_gs_h1<-grep(pattern="Scaffold_11__3_contigs__length_97865747",x=gff_refug_gs_h1$V1,fixed=TRUE)


## subset to ch8
ch8_h154_gus<-gff_h154_gus_h2[ch8i_h154_gus_h2,]
ch8_h154_gs<-gff_h154_gs_h1[ch8i_h154_gs_h1,]
ch8_refug_gus<-gff_refug_gus_h1[ch8i_refug_gus_h1,]
ch8_refug_gs<-gff_refug_gs_h1[ch8i_refug_gs_h1,]

## bounds whole region
## buffer 
buffer<-50000
h154_gs_lb<-24457103-buffer;h154_gs_ub<-39030359+buffer
h154_gus_lb<-24803527-buffer;h154_gus_ub<-44121870+buffer

refug_gs_lb<-22442098-buffer;refug_gs_ub<-65729704+buffer
refug_gus_lb<-22220178-buffer;refug_gus_ub<-65829835+buffer

## inverted trans only
h154_gs_it_lb<-24457103;h154_gs_it_ub<-32202650
h154_gus_it_lb<-33112541;h154_gus_it_ub<-44121870

refug_gs_it_lb<-56898118;refug_gs_it_ub<-65729704
refug_gus_it_lb<-22220178;refug_gus_it_ub<-31457192

## indexes for region
li_h154_gus<-which(ch8_h154_gus$V4 >= h154_gus_lb & ch8_h154_gus$V5 <= h154_gus_ub)
li_h154_gs<-which(ch8_h154_gs$V4 >= h154_gs_lb & ch8_h154_gs$V5 <= h154_gs_ub)
li_refug_gus<-which(ch8_refug_gus$V4 >= refug_gus_lb & ch8_refug_gus$V5 <= refug_gus_ub)
li_refug_gs<-which(ch8_refug_gs$V4 >= refug_gs_lb & ch8_refug_gs$V5 <= refug_gs_ub)


loc_h154_gus<-ch8_h154_gus[li_h154_gus,]
loc_h154_gs<-ch8_h154_gs[li_h154_gs,]
loc_refug_gus<-ch8_refug_gus[li_refug_gus,]
loc_refug_gs<-ch8_refug_gs[li_refug_gs,]

########### H154 #####################

## totals
sum(loc_h154_gus$V5-loc_h154_gus$V4)
sum(loc_h154_gs$V5-loc_h154_gs$V4)
h154_gus_ub-h154_gus_lb
h154_gs_ub-h154_gs_lb
## locus is about 1.3 x larger for green than stripe, and about 1.4 x for repeat content for green than stripe, so similar proportions
sum(loc_h154_gus$V5-loc_h154_gus$V4)/(h154_gus_ub-h154_gus_lb)
#[1] 0.6227906= 62% repeat
sum(loc_h154_gs$V5-loc_h154_gs$V4)/(h154_gs_ub-h154_gs_lb)
#[1] 0.5932954 = 59% repeat


## but what about types
## indexes for simple repeats
gs_h154_sim_rep<-grep(pattern=")n",x=loc_h154_gs$V10)
gus_h154_sim_rep<-grep(pattern=")n",x=loc_h154_gus$V10)
sum(loc_h154_gus$V5[gus_h154_sim_rep]-loc_h154_gus$V4[gus_h154_sim_rep])/(h154_gus_ub-h154_gus_lb)
sum(loc_h154_gs$V5[gs_h154_sim_rep]-loc_h154_gs$V4[gs_h154_sim_rep])/(h154_gs_ub-h154_gs_lb)
loc_h154_gs$V10[gs_h154_sim_rep]<-"SimpleRep"
loc_h154_gus$V10[gus_h154_sim_rep]<-"SimpleRep"


## answer is not simple repeats as per above, only about 1% simple repeats period, but other forms of simple repeats still included below

## try rest
gs_types<-gsub(pattern="^Motif:([a-zA-Z-/]+)[0-9]+_[0-9a-zA-Z-]+",x=loc_h154_gs$V10,perl=TRUE,replacement="\\1",fixed=FALSE)
gs_types<-gsub(pattern="^Motif:([a-zA-Z0-9/]+)[0-9a-zA-Z-]+_Tcr",x=gs_types,perl=TRUE,replacement="\\1",fixed=FALSE)
gs_types<-gsub(pattern="-$",x=gs_types,perl=TRUE,replacement="",fixed=FALSE)
gs_types<-gsub(pattern="^Motif:",x=gs_types,perl=TRUE,replacement="",fixed=FALSE)
gs_h154_sran_rep<-which(gs_types %in% c("SimpleRep","A-rich","G-rich","GA-rich","rnd"))
gs_type_cnts<-tapply(X=loc_h154_gs$V5[-gs_h154_sran_rep]-loc_h154_gs$V4[-gs_h154_sran_rep],INDEX=gs_types[-gs_h154_sran_rep],sum)


gus_types<-gsub(pattern="^Motif:([a-zA-Z-/]+)[0-9]+_[0-9a-zA-Z-]+",x=loc_h154_gus$V10,perl=TRUE,replacement="\\1",fixed=FALSE)
gus_types<-gsub(pattern="^Motif:([a-zA-Z0-9/]+)[0-9a-zA-Z-]+_Tcr",x=gus_types,perl=TRUE,replacement="\\1",fixed=FALSE)
gus_types<-gsub(pattern="-$",x=gus_types,perl=TRUE,replacement="",fixed=FALSE)
gus_types<-gsub(pattern="^Motif:",x=gus_types,perl=TRUE,replacement="",fixed=FALSE)
gus_h154_sran_rep<-which(gus_types %in% c("SimpleRep","A-rich","G-rich","GA-rich","rnd"))
gus_type_cnts<-tapply(X=loc_h154_gus$V5[-gus_h154_sran_rep]-loc_h154_gus$V4[-gus_h154_sran_rep],INDEX=gus_types[-gus_h154_sran_rep],sum)


dotchart(gs_type_cnts, cex = 0.5) 
dotchart(gus_type_cnts,cex=.5)
################################################## for main figure ##########
pdf("TE.pdf",width=9,height=6)
par(mfrow=c(2,3))
par(mar=c(4.5,5.5,2.5,1.5))
cl<-1.45
ca<-1.0
cm<-1.3
dotchart(gs_type_cnts, cex = 0.5,pch=19,cex.lab=cl,cex.axis=ca,xlab="Abundance (bps)")
title(main="(E) Stripe repetitive elements",cex.main=cm)
dotchart(gus_type_cnts,cex=.5,pch=19,cex.lab=cl,cex.axis=ca,xlab="Abundance (bps)")
title(main="(F) Green repetitive elements",cex.main=cm)
pbuffer<-50000
lb<-32202650-pbuffer;ub<-32202951+pbuffer
plot(c(lb,ub),c(0,1),type='n',xlab="Position (kbps)",ylab="",axes=FALSE,cex.lab=cl)
axis(1,at=seq(lb,ub,50000),round(seq(lb,ub,50000)/1000,0))
xxi<-which(loc_h154_gs$V4 < ub & loc_h154_gs$V5 > lb)
xx<-gs_types[which(loc_h154_gs$V4 < ub & loc_h154_gs$V5 > lb)]
a<-which(xx %in%  c("SimpleRep","A-rich","G-rich","GA-rich","rnd"))
b<-which(!xx %in%  c("SimpleRep","A-rich","G-rich","GA-rich","rnd"))
A<-length(a)
B<-length(b)
#for(i in 1:A){
#	lines(c(loc_h154_gs$V4[xxi[a][i]],loc_h154_gs$V5[xxi[a][i]]),c(.2,.2))
#}
for(i in 1:B){
	lines(c(loc_h154_gs$V4[xxi[b][i]],loc_h154_gs$V5[xxi[b][i]]),c(.5,.5),col="firebrick")
}
genes<-which(genes_h154_gs_h1$V3=="gene" & genes_h154_gs_h1$V1 == "Scaffold_11__2_contigs__length_91821751" & genes_h154_gs_h1$V4 < ub & genes_h154_gs_h1$V5 > lb)
if(length(genes) > 0){
	for(i in 1:length(genes)){
		lines(c(genes_h154_gs_h1$V4[genes[i]],genes_h154_gs_h1$V5[genes[i]]),c(.7,.7),col="cadetblue")
	}
}
abline(v=c(32202650,32202951),lty=2)

title(main="(G) Breakpoint boundary",cex.main=cm)
dev.off()
###############################################################
L_rdat_gs_h1<-vector("list",3)

## h154 gs h1
genes_h154_gs_h1<-read.table("Annotation/t_crist_hyw154_stripe_h1/braker/braker.gff3",header=FALSE,fill=TRUE,comment.char="#")

pbuffer<-10000
lb<-24457103-pbuffer;ub<-24457103+pbuffer
plot(c(lb,ub),c(0,1),type='n',xlab="Position (bps)",ylab="")
xxi<-which(loc_h154_gs$V4 < ub & loc_h154_gs$V5 > lb)
xx<-gs_types[which(loc_h154_gs$V4 < ub & loc_h154_gs$V5 > lb)]
a<-which(xx %in%  c("SimpleRep","A-rich","G-rich","GA-rich","rnd"))
b<-which(!xx %in%  c("SimpleRep","A-rich","G-rich","GA-rich","rnd"))
A<-length(a)
B<-length(b)
for(i in 1:A){
	lines(c(loc_h154_gs$V4[xxi[a][i]],loc_h154_gs$V5[xxi[a][i]]),c(.2,.2))
}
for(i in 1:B){
	lines(c(loc_h154_gs$V4[xxi[b][i]],loc_h154_gs$V5[xxi[b][i]]),c(.4,.4),col="firebrick")
}
genes<-which(genes_h154_gs_h1$V3=="gene" & genes_h154_gs_h1$V1 == "Scaffold_11__2_contigs__length_91821751" & genes_h154_gs_h1$V4 < ub & genes_h154_gs_h1$V5 > lb)
if(length(genes) > 0){
	for(i in 1:length(genes)){
		lines(c(genes_h154_gs_h1$V4[genes[i]],genes_h154_gs_h1$V5[genes[i]]),c(.7,.7),col="cadetblue")
	}
}
abline(v=24457103,lty=2)

dat<-cbind(loc_h154_gs$V4[xxi[b]],loc_h154_gs$V5[xxi[b]])
rdat<-dat
lb<-24457103;ub<-24457103
rdat[dat<lb] <- dat[dat<lb]-lb
rdat[dat>ub] <- dat[dat>ub]-ub

## need to flip so this is "right" boundary not left and needs negative
L_rdat_gs_h1[[3]]<--1 * rdat




lb<-32202650-pbuffer;ub<-32202951+pbuffer
plot(c(lb,ub),c(0,1),type='n',xlab="Position (bps)",ylab="")
xxi<-which(loc_h154_gs$V4 < ub & loc_h154_gs$V5 > lb)
xx<-gs_types[which(loc_h154_gs$V4 < ub & loc_h154_gs$V5 > lb)]
a<-which(xx %in%  c("SimpleRep","A-rich","G-rich","GA-rich","rnd"))
b<-which(!xx %in%  c("SimpleRep","A-rich","G-rich","GA-rich","rnd"))
A<-length(a)
B<-length(b)
for(i in 1:A){
	lines(c(loc_h154_gs$V4[xxi[a][i]],loc_h154_gs$V5[xxi[a][i]]),c(.2,.2))
}
for(i in 1:B){
	lines(c(loc_h154_gs$V4[xxi[b][i]],loc_h154_gs$V5[xxi[b][i]]),c(.4,.4),col="firebrick")
}
genes<-which(genes_h154_gs_h1$V3=="gene" & genes_h154_gs_h1$V1 == "Scaffold_11__2_contigs__length_91821751" & genes_h154_gs_h1$V4 < ub & genes_h154_gs_h1$V5 > lb)
if(length(genes) > 0){
	for(i in 1:length(genes)){
		lines(c(genes_h154_gs_h1$V4[genes[i]],genes_h154_gs_h1$V5[genes[i]]),c(.7,.7),col="cadetblue")
	}
}
abline(v=c(32202650,32202951),lty=2)
## HAT then Gypsy, HAT is directly adjacent and all are touching

dat<-cbind(loc_h154_gs$V4[xxi[b]],loc_h154_gs$V5[xxi[b]])
rdat<-dat
lb<-32202650;ub<-32202951
rdat[dat<lb] <- dat[dat<lb]-lb
rdat[dat>ub] <- dat[dat>ub]-ub
rdat[dat >= lb & dat <=ub]<-0

## need to flip, still internal but needs negative
L_rdat_gs_h1[[2]]<--1 * rdat

lb<-39030359-pbuffer;ub<-39030359+pbuffer
plot(c(lb,ub),c(0,1),type='n',xlab="Position (bps)",ylab="")
xxi<-which(loc_h154_gs$V4 < ub & loc_h154_gs$V5 > lb)
xx<-gs_types[which(loc_h154_gs$V4 < ub & loc_h154_gs$V5 > lb)]
a<-which(xx %in%  c("SimpleRep","A-rich","G-rich","GA-rich","rnd"))
b<-which(!xx %in%  c("SimpleRep","A-rich","G-rich","GA-rich","rnd"))
A<-length(a)
B<-length(b)
for(i in 1:A){
	lines(c(loc_h154_gs$V4[xxi[a][i]],loc_h154_gs$V5[xxi[a][i]]),c(.2,.2))
}
for(i in 1:B){
	lines(c(loc_h154_gs$V4[xxi[b][i]],loc_h154_gs$V5[xxi[b][i]]),c(.4,.4),col="firebrick")
}
if(length(genes) > 0){
	for(i in 1:length(genes)){
		lines(c(genes_h154_gs_h1$V4[genes[i]],genes_h154_gs_h1$V5[genes[i]]),c(.7,.7),col="cadetblue")
	}
}
abline(v=c(39030359),lty=2)
## RTE and Mariner, but not touching bounds
dat<-cbind(loc_h154_gs$V4[xxi[b]],loc_h154_gs$V5[xxi[b]])
rdat<-dat
lb<-39030359;ub<-39030359
rdat[dat<lb] <- dat[dat<lb]-lb
rdat[dat>ub] <- dat[dat>ub]-ub
rdat[dat >= lb & dat <=ub]<-0

## need to flip, left not right and needs negative
L_rdat_gs_h1[[1]]<--1 * rdat


## h154 gus h2
L_rdat_gus_h2<-vector("list",3)

genes_h154_gus_h2<-read.table("Annotation/t_crist_hyw154_green_h2/braker/braker.gff3",header=FALSE,fill=TRUE,comment.char="#")

pbuffer<-10000
lb<-24803527-pbuffer;ub<-24803527+pbuffer
plot(c(lb,ub),c(0,1),type='n',xlab="Position (bps)",ylab="")
xxi<-which(loc_h154_gus$V4 < ub & loc_h154_gus$V5 > lb)
xx<-gus_types[which(loc_h154_gus$V4 < ub & loc_h154_gus$V5 > lb)]
a<-which(xx %in%  c("SimpleRep","A-rich","G-rich","GA-rich","rnd"))
b<-which(!xx %in%  c("SimpleRep","A-rich","G-rich","GA-rich","rnd"))
A<-length(a)
B<-length(b)
for(i in 1:A){
	lines(c(loc_h154_gus$V4[xxi[a][i]],loc_h154_gus$V5[xxi[a][i]]),c(.2,.2))
}
for(i in 1:B){
	lines(c(loc_h154_gus$V4[xxi[b][i]],loc_h154_gus$V5[xxi[b][i]]),c(.4,.4),col="firebrick")
}
genes<-which(genes_h154_gus_h2$V3=="gene" & genes_h154_gus_h2$V1 == "ScrX45T_23_HRSCAF_264" & genes_h154_gus_h2$V4 < ub & genes_h154_gus_h2$V5 > lb)
if(length(genes) > 0){
	for(i in 1:length(genes)){
		lines(c(genes_h154_gus_h2$V4[genes[i]],genes_h154_gus_h2$V5[genes[i]]),c(.7,.7),col="cadetblue")
	}
}
abline(v=24803527,lty=2)
## gene not known, far left TE is LINE/SINEHRC rest HAT

dat<-cbind(loc_h154_gus$V4[xxi[b]],loc_h154_gus$V5[xxi[b]])
rdat<-dat
lb<-24803527;ub<-24803527
rdat[dat<lb] <- dat[dat<lb]-lb
rdat[dat>ub] <- dat[dat>ub]-ub
rdat[dat >= lb & dat <=ub]<-0

## need to flip, right not left and needs negative
L_rdat_gus_h2[[3]]<--1 * rdat


lb<-32951302-pbuffer;ub<-32951302+pbuffer
plot(c(lb,ub),c(0,1),type='n',xlab="Position (bps)",ylab="")
xxi<-which(loc_h154_gus$V4 < ub & loc_h154_gus$V5 > lb)
xx<-gus_types[which(loc_h154_gus$V4 < ub & loc_h154_gus$V5 > lb)]
a<-which(xx %in%  c("SimpleRep","A-rich","G-rich","GA-rich","rnd"))
b<-which(!xx %in%  c("SimpleRep","A-rich","G-rich","GA-rich","rnd"))
A<-length(a)
B<-length(b)
for(i in 1:A){
	lines(c(loc_h154_gus$V4[xxi[a][i]],loc_h154_gus$V5[xxi[a][i]]),c(.2,.2))
}
for(i in 1:B){
	lines(c(loc_h154_gus$V4[xxi[b][i]],loc_h154_gus$V5[xxi[b][i]]),c(.4,.4),col="firebrick")
}
genes<-which(genes_h154_gus_h2$V3=="gene" & genes_h154_gus_h2$V1 == "ScrX45T_23_HRSCAF_264" & genes_h154_gus_h2$V4 < ub & genes_h154_gus_h2$V5 > lb)
if(length(genes) > 0){
	for(i in 1:length(genes)){
		lines(c(genes_h154_gus_h2$V4[genes[i]],genes_h154_gus_h2$V5[genes[i]]),c(.7,.7),col="cadetblue")
	}
}
abline(v=c(32951302,32951302),lty=2)
#[1] "Academ"   "Gypsy"    "Gypsy"    "GypsyHRC"
# but not touchi-ng
dat<-cbind(loc_h154_gus$V4[xxi[b]],loc_h154_gus$V5[xxi[b]])
rdat<-dat
lb<-32951302;ub<-32951302
rdat[dat<lb] <- dat[dat<lb]-lb
rdat[dat>ub] <- dat[dat>ub]-ub
rdat[dat >= lb & dat <=ub]<-0

## need to flip, still internal but needs negative
L_rdat_gus_h2[[2]]<--1 * rdat


lb<-44121870-pbuffer;ub<-44121870+pbuffer
plot(c(lb,ub),c(0,1),type='n',xlab="Position (bps)",ylab="")
xxi<-which(loc_h154_gus$V4 < ub & loc_h154_gus$V5 > lb)
xx<-gus_types[which(loc_h154_gus$V4 < ub & loc_h154_gus$V5 > lb)]
a<-which(xx %in%  c("SimpleRep","A-rich","G-rich","GA-rich","rnd"))
b<-which(!xx %in%  c("SimpleRep","A-rich","G-rich","GA-rich","rnd"))
A<-length(a)
B<-length(b)
for(i in 1:A){
	lines(c(loc_h154_gus$V4[xxi[a][i]],loc_h154_gus$V5[xxi[a][i]]),c(.2,.2))
}
for(i in 1:B){
	lines(c(loc_h154_gus$V4[xxi[b][i]],loc_h154_gus$V5[xxi[b][i]]),c(.4,.4),col="firebrick")
}
genes<-which(genes_h154_gus_h2$V3=="gene" & genes_h154_gus_h2$V1 == "ScrX45T_23_HRSCAF_264" & genes_h154_gus_h2$V4 < ub & genes_h154_gus_h2$V5 > lb)
if(length(genes) > 0){
	for(i in 1:length(genes)){
		lines(c(genes_h154_gus_h2$V4[genes[i]],genes_h154_gus_h2$V5[genes[i]]),c(.7,.7),col="cadetblue")
	}
}
abline(v=c(44121870),lty=2)
## [1] "HAT"      "HAT"      "HAT"      "GypsyHRC"
## but not touching
## gene = lipoyl synthase like

dat<-cbind(loc_h154_gus$V4[xxi[b]],loc_h154_gus$V5[xxi[b]])
rdat<-dat
lb<-44121870;ub<-44121870
rdat[dat<lb] <- dat[dat<lb]-lb
rdat[dat>ub] <- dat[dat>ub]-ub
rdat[dat >= lb & dat <=ub]<-0

## need to flip, left not right and needs negative
L_rdat_gus_h2[[1]]<--1 * rdat



########### refugio #####################

## totals
sum(loc_refug_gus$V5-loc_refug_gus$V4)
sum(loc_refug_gs$V5-loc_refug_gs$V4)
refug_gus_ub-refug_gus_lb
refug_gs_ub-refug_gs_lb
## locus is about 1.3 x larger for green than stripe, and about 1.4 x for repeat content for green than stripe, so similar proportions
sum(loc_refug_gus$V5-loc_refug_gus$V4)/(refug_gus_ub-refug_gus_lb)
#[1] 0.5984133 = 60% repeat
sum(loc_refug_gs$V5-loc_refug_gs$V4)/(refug_gs_ub-refug_gs_lb)
#[1] 0.5982367 = 60% repeat


## but what about types
## indexes for simple repeats
gs_refug_sim_rep<-grep(pattern=")n",x=loc_refug_gs$V10)
gus_refug_sim_rep<-grep(pattern=")n",x=loc_refug_gus$V10)
sum(loc_refug_gus$V5[gus_refug_sim_rep]-loc_refug_gus$V4[gus_refug_sim_rep])/(refug_gus_ub-refug_gus_lb)
#[1] 0.01400432
sum(loc_refug_gs$V5[gs_refug_sim_rep]-loc_refug_gs$V4[gs_refug_sim_rep])/(refug_gs_ub-refug_gs_lb)
#[1] 0.01292157
loc_refug_gs$V10[gs_refug_sim_rep]<-"SimpleRep"
loc_refug_gus$V10[gus_refug_sim_rep]<-"SimpleRep"

## answer is not simple repeats as per above, only about 1% simple repeats period, but other forms of simple repeats still included below

## try rest
gs_types<-gsub(pattern="^Motif:([a-zA-Z-/]+)[0-9]+_[0-9a-zA-Z-]+",x=loc_refug_gs$V10,perl=TRUE,replacement="\\1",fixed=FALSE)
gs_types<-gsub(pattern="^Motif:([a-zA-Z0-9/]+)[0-9a-zA-Z-]+_Tcr",x=gs_types,perl=TRUE,replacement="\\1",fixed=FALSE)
gs_types<-gsub(pattern="-$",x=gs_types,perl=TRUE,replacement="",fixed=FALSE)
gs_types<-gsub(pattern="^Motif:",x=gs_types,perl=TRUE,replacement="",fixed=FALSE)
gs_refug_sran_rep<-which(gs_types %in% c("SimpleRep","A-rich","G-rich","GA-rich","rnd"))
gs_type_cnts<-tapply(X=loc_refug_gs$V5[-gs_refug_sran_rep]-loc_refug_gs$V4[-gs_refug_sran_rep],INDEX=gs_types[-gs_refug_sran_rep],sum)


gus_types<-gsub(pattern="^Motif:([a-zA-Z-/]+)[0-9]+_[0-9a-zA-Z-]+",x=loc_refug_gus$V10,perl=TRUE,replacement="\\1",fixed=FALSE)
gus_types<-gsub(pattern="^Motif:([a-zA-Z0-9/]+)[0-9a-zA-Z-]+_Tcr",x=gus_types,perl=TRUE,replacement="\\1",fixed=FALSE)
gus_types<-gsub(pattern="-$",x=gus_types,perl=TRUE,replacement="",fixed=FALSE)
gus_types<-gsub(pattern="^Motif:",x=gus_types,perl=TRUE,replacement="",fixed=FALSE)
gus_refug_sran_rep<-which(gus_types %in% c("SimpleRep","A-rich","G-rich","GA-rich","rnd"))
gus_type_cnts<-tapply(X=loc_refug_gus$V5[-gus_refug_sran_rep]-loc_refug_gus$V4[-gus_refug_sran_rep],INDEX=gus_types[-gus_refug_sran_rep],sum)

dotchart(gs_type_cnts, cex = 0.5) 
dotchart(gus_type_cnts,cex=.5)



## refug gs h1
genes_refug_gs_h1<-read.table("Annotation/t_crist_refug_stripe_h1/braker/braker.gff3",header=FALSE,fill=TRUE,comment.char="#")

L_rdat_gs_r1<-vector("list",3)

pbuffer<-10000
lb<-22442098-pbuffer;ub<-22442098+pbuffer
plot(c(lb,ub),c(0,1),type='n',xlab="Position (bps)",ylab="")
xxi<-which(loc_refug_gs$V4 < ub & loc_refug_gs$V5 > lb)
xx<-gs_types[which(loc_refug_gs$V4 < ub & loc_refug_gs$V5 > lb)]
a<-which(xx %in%  c("SimpleRep","A-rich","G-rich","GA-rich","rnd"))
b<-which(!xx %in%  c("SimpleRep","A-rich","G-rich","GA-rich","rnd"))
A<-length(a)
B<-length(b)
for(i in 1:A){
        lines(c(loc_refug_gs$V4[xxi[a][i]],loc_refug_gs$V5[xxi[a][i]]),c(.2,.2))
}
for(i in 1:B){
        lines(c(loc_refug_gs$V4[xxi[b][i]],loc_refug_gs$V5[xxi[b][i]]),c(.4,.4),col="firebrick")
}
genes<-which(genes_refug_gs_h1$V3=="gene" & genes_refug_gs_h1$V1 == "Scaffold_11__3_contigs__length_97865747" & genes_refug_gs_h1$V4 < ub & genes_refug_gs_h1$V5 > lb)
if(length(genes) > 0){
        for(i in 1:length(genes)){
                lines(c(genes_refug_gs_h1$V4[genes[i]],genes_refug_gs_h1$V5[genes[i]]),c(.7,.7),col="cadetblue")
        }
}
abline(v=22442098,lty=2)
#[1] "RTE"         "Copia"       "Copia"       "Copia"       "Penelope"   
#[6] "HelitronHRC" "HelitronHRC"
## copia on break

dat<-cbind(loc_refug_gs$V4[xxi[b]],loc_refug_gs$V5[xxi[b]])
rdat<-dat
lb<-22442098;ub<-22442098
rdat[dat<lb] <- dat[dat<lb]-lb
rdat[dat>ub] <- dat[dat>ub]-ub
rdat[dat >= lb & dat <=ub]<-0

## good as is
L_rdat_gs_r1[[1]]<-rdat


lb<-56898118-pbuffer;ub<-56898118+pbuffer
plot(c(lb,ub),c(0,1),type='n',xlab="Position (bps)",ylab="")
xxi<-which(loc_refug_gs$V4 < ub & loc_refug_gs$V5 > lb)
xx<-gs_types[which(loc_refug_gs$V4 < ub & loc_refug_gs$V5 > lb)]
a<-which(xx %in%  c("SimpleRep","A-rich","G-rich","GA-rich","rnd"))
b<-which(!xx %in%  c("SimpleRep","A-rich","G-rich","GA-rich","rnd"))
A<-length(a)
B<-length(b)
for(i in 1:A){
        lines(c(loc_refug_gs$V4[xxi[a][i]],loc_refug_gs$V5[xxi[a][i]]),c(.2,.2))
}
for(i in 1:B){
        lines(c(loc_refug_gs$V4[xxi[b][i]],loc_refug_gs$V5[xxi[b][i]]),c(.4,.4),col="firebrick")
}
genes<-which(genes_refug_gs_h1$V3=="gene" & genes_refug_gs_h1$V1 == "Scaffold_11__3_contigs__length_97865747" & genes_refug_gs_h1$V4 < ub & genes_refug_gs_h1$V5 > lb)
if(length(genes) > 0){
        for(i in 1:length(genes)){
                lines(c(genes_refug_gs_h1$V4[genes[i]],genes_refug_gs_h1$V5[genes[i]]),c(.7,.7),col="cadetblue")
        }
}
abline(v=56898118,lty=2)
## nothing
dat<-cbind(loc_refug_gs$V4[xxi[b]],loc_refug_gs$V5[xxi[b]])
rdat<-dat
lb<-56898118;ub<-56898118
rdat[dat<lb] <- dat[dat<lb]-lb
rdat[dat>ub] <- dat[dat>ub]-ub
rdat[dat >= lb & dat <=ub]<-0

## good as is
L_rdat_gs_r1[[2]]<-rdat

lb<-65729704-pbuffer;ub<-65729704+pbuffer
plot(c(lb,ub),c(0,1),type='n',xlab="Position (bps)",ylab="")
xxi<-which(loc_refug_gs$V4 < ub & loc_refug_gs$V5 > lb)
xx<-gs_types[which(loc_refug_gs$V4 < ub & loc_refug_gs$V5 > lb)]
a<-which(xx %in%  c("SimpleRep","A-rich","G-rich","GA-rich","rnd"))
b<-which(!xx %in%  c("SimpleRep","A-rich","G-rich","GA-rich","rnd"))
A<-length(a)
B<-length(b)
for(i in 1:A){
        lines(c(loc_refug_gs$V4[xxi[a][i]],loc_refug_gs$V5[xxi[a][i]]),c(.2,.2))
}
for(i in 1:B){
        lines(c(loc_refug_gs$V4[xxi[b][i]],loc_refug_gs$V5[xxi[b][i]]),c(.4,.4),col="firebrick")
}
genes<-which(genes_refug_gs_h1$V3=="gene" & genes_refug_gs_h1$V1 == "Scaffold_11__3_contigs__length_97865747" & genes_refug_gs_h1$V4 < ub & genes_refug_gs_h1$V5 > lb)
if(length(genes) > 0){
        for(i in 1:length(genes)){
                lines(c(genes_refug_gs_h1$V4[genes[i]],genes_refug_gs_h1$V5[genes[i]]),c(.7,.7),col="cadetblue")
        }
}
abline(v=c(65729704),lty=2)
## [1] "RTE"     "Mariner" "HAT"     "HAT"     "HAT" 
## HAT Is near boundary, a few kb away
dat<-cbind(loc_refug_gs$V4[xxi[b]],loc_refug_gs$V5[xxi[b]])
rdat<-dat
lb<-65729704;ub<-65729704
rdat[dat<lb] <- dat[dat<lb]-lb
rdat[dat>ub] <- dat[dat>ub]-ub
rdat[dat >= lb & dat <=ub]<-0

## good as is
L_rdat_gs_r1[[3]]<-rdat

## refug gus h1
genes_refug_gus_h1<-read.table("Annotation/t_crist_refug_green_h1/braker/braker.gff3",header=FALSE,fill=TRUE,comment.char="#")
L_rdat_gus_r1<-vector("list",3)

pbuffer<-10000
lb<-22220178-pbuffer;ub<-22220178+pbuffer
plot(c(lb,ub),c(0,1),type='n',xlab="Position (bps)",ylab="")
xxi<-which(loc_refug_gus$V4 < ub & loc_refug_gus$V5 > lb)
xx<-gus_types[which(loc_refug_gus$V4 < ub & loc_refug_gus$V5 > lb)]
a<-which(xx %in%  c("SimpleRep","A-rich","G-rich","GA-rich","rnd"))
b<-which(!xx %in%  c("SimpleRep","A-rich","G-rich","GA-rich","rnd"))
A<-length(a)
B<-length(b)
for(i in 1:A){
        lines(c(loc_refug_gus$V4[xxi[a][i]],loc_refug_gus$V5[xxi[a][i]]),c(.2,.2))
}
for(i in 1:B){
        lines(c(loc_refug_gus$V4[xxi[b][i]],loc_refug_gus$V5[xxi[b][i]]),c(.4,.4),col="firebrick")
}
genes<-which(genes_refug_gus_h1$V3=="gene" & genes_refug_gus_h1$V1 == "Scaffold_11__3_contigs__length_97865747" & genes_refug_gus_h1$V4 < ub & genes_refug_gus_h1$V5 > lb)
if(length(genes) > 0){
        for(i in 1:length(genes)){
                lines(c(genes_refug_gus_h1$V4[genes[i]],genes_refug_gus_h1$V5[genes[i]]),c(.7,.7),col="cadetblue")
        }
}
abline(v=22220178,lty=2)
# gus_types[xxi[b]]
#[1] "Academ"   "RTE"      "Gypsy"    "Polinton" "Polinton" "HAT"      "Polinton"
## spanned by polinton
dat<-cbind(loc_refug_gus$V4[xxi[b]],loc_refug_gus$V5[xxi[b]])
rdat<-dat
lb<-22220178;ub<-22220178
rdat[dat<lb] <- dat[dat<lb]-lb
rdat[dat>ub] <- dat[dat>ub]-ub
rdat[dat >= lb & dat <=ub]<-0

## good as is
L_rdat_gus_r1[[1]]<-rdat



lb<-31457192-pbuffer;ub<-31523304+pbuffer
plot(c(lb,ub),c(0,1),type='n',xlab="Position (bps)",ylab="")
xxi<-which(loc_refug_gus$V4 < ub & loc_refug_gus$V5 > lb)
xx<-gus_types[which(loc_refug_gus$V4 < ub & loc_refug_gus$V5 > lb)]
a<-which(xx %in%  c("SimpleRep","A-rich","G-rich","GA-rich","rnd"))
b<-which(!xx %in%  c("SimpleRep","A-rich","G-rich","GA-rich","rnd"))
A<-length(a)
B<-length(b)
for(i in 1:A){
        lines(c(loc_refug_gus$V4[xxi[a][i]],loc_refug_gus$V5[xxi[a][i]]),c(.2,.2))
}
for(i in 1:B){
        lines(c(loc_refug_gus$V4[xxi[b][i]],loc_refug_gus$V5[xxi[b][i]]),c(.4,.4),col="firebrick")
}
genes<-which(genes_refug_gus_h1$V3=="gene" & genes_refug_gus_h1$V1 == "Scaffold_11__3_contigs__length_97865747" & genes_refug_gus_h1$V4 < ub & genes_refug_gus_h1$V5 > lb)
if(length(genes) > 0){
        for(i in 1:length(genes)){
                lines(c(genes_refug_gus_h1$V4[genes[i]],genes_refug_gus_h1$V5[genes[i]]),c(.7,.7),col="cadetblue")
        }
}
abline(v=c(31457192,31523304),lty=2)
unique(gus_types[xxi[b]])
#[1] "HelitronHRC" "HAT"         "RTEHRC"      "Gypsy"       "Harbinger"  
#[6] "Mariner"     "Ginger"      "Copia"       "Penelope"
## copia touching left bounds
dat<-cbind(loc_refug_gus$V4[xxi[b]],loc_refug_gus$V5[xxi[b]])
rdat<-dat
lb<-31457192;ub<-31523304
rdat[dat<lb] <- dat[dat<lb]-lb
rdat[dat>ub] <- dat[dat>ub]-ub
rdat[dat >= lb & dat <=ub]<-0

## good as is
L_rdat_gus_r1[[2]]<-rdat


lb<-65829835-pbuffer;ub<-65829835+pbuffer
plot(c(lb,ub),c(0,1),type='n',xlab="Position (bps)",ylab="")
xxi<-which(loc_refug_gus$V4 < ub & loc_refug_gus$V5 > lb)
xx<-gus_types[which(loc_refug_gus$V4 < ub & loc_refug_gus$V5 > lb)]
a<-which(xx %in%  c("SimpleRep","A-rich","G-rich","GA-rich","rnd"))
b<-which(!xx %in%  c("SimpleRep","A-rich","G-rich","GA-rich","rnd"))
A<-length(a)
B<-length(b)
for(i in 1:A){
        lines(c(loc_refug_gus$V4[xxi[a][i]],loc_refug_gus$V5[xxi[a][i]]),c(.2,.2))
}
for(i in 1:B){
        lines(c(loc_refug_gus$V4[xxi[b][i]],loc_refug_gus$V5[xxi[b][i]]),c(.4,.4),col="firebrick")
}
genes<-which(genes_refug_gus_h1$V3=="gene" & genes_refug_gus_h1$V1 == "Scaffold_11__3_contigs__length_97865747" & genes_refug_gus_h1$V4 < ub & genes_refug_gus_h1$V5 > lb)
if(length(genes) > 0){
        for(i in 1:length(genes)){
                lines(c(genes_refug_gus_h1$V4[genes[i]],genes_refug_gus_h1$V5[genes[i]]),c(.7,.7),col="cadetblue")
        }
}
abline(v=c(65829835),lty=2)
## unique(gus_types[xxi[b]])
#[1] "HAT"     "Mariner" "Copia"   "Kolobok"

## nothing quite touching
dat<-cbind(loc_refug_gus$V4[xxi[b]],loc_refug_gus$V5[xxi[b]])
rdat<-dat
lb<-65829835;ub<-65829835
rdat[dat<lb] <- dat[dat<lb]-lb
rdat[dat>ub] <- dat[dat>ub]-ub
rdat[dat >= lb & dat <=ub]<-0

## good as is
L_rdat_gus_r1[[3]]<-rdat

save(list=ls(),file="TE.rdat")
################################################## for main figure, version 2 ##########
pdf("TE.pdf",width=9,height=3)
par(mfrow=c(1,3))
par(mar=c(4.5,5.5,2.5,1.5))
cl<-1.4
ca<-1.0
cm<-1.25
lw<-4
bnd<-5000

plot(c(-1*bnd,bnd),c(.8,4.2),type='n',xlab="Position (bps)",ylab="Reference genome",axes=FALSE,cex.lab=cl)
axis(1,at=seq(-1*bnd,bnd,bnd/2))
k<-1
N<-dim(L_rdat_gs_r1[[k]])[1]
for(i in 1:N){
	lines(c(L_rdat_gs_r1[[k]][i,1],L_rdat_gs_r1[[k]][i,2]),c(4,4),col="firebrick",lwd=lw)
}
k<-1
N<-dim(L_rdat_gus_r1[[k]])[1]
for(i in 1:N){
	lines(c(L_rdat_gus_r1[[k]][i,1],L_rdat_gus_r1[[k]][i,2]),c(3,3),col="firebrick",lwd=lw)
}
k<-1
N<-dim(L_rdat_gs_h1[[k]])[1]
for(i in 1:N){
	lines(c(L_rdat_gs_h1[[k]][i,1],L_rdat_gs_h1[[k]][i,2]),c(2,2),col="firebrick",lwd=lw)
}
k<-1
N<-dim(L_rdat_gus_h2[[k]])[1]
for(i in 1:N){
	lines(c(L_rdat_gus_h2[[k]][i,1],L_rdat_gus_h2[[k]][i,2]),c(1,1),col="firebrick",lwd=lw)
}
abline(v=0,lty=2,lwd=2)

axis(2,at=c(1:4),c("154 G","154 S","Ref G","Ref S"),cex.axis=ca)

title(main="(E) Breakpoint A",cex.main=cm)

plot(c(-1*bnd,bnd),c(.8,4.2),type='n',xlab="Position (bps)",ylab="Reference genome",axes=FALSE,cex.lab=cl)
axis(1,at=seq(-1*bnd,bnd,bnd/2))
k<-2
N<-dim(L_rdat_gs_r1[[k]])[1]
for(i in 1:N){
	lines(c(L_rdat_gs_r1[[k]][i,1],L_rdat_gs_r1[[k]][i,2]),c(4,4),col="firebrick",lwd=lw)
}
k<-2
N<-dim(L_rdat_gus_r1[[k]])[1]
for(i in 1:N){
	lines(c(L_rdat_gus_r1[[k]][i,1],L_rdat_gus_r1[[k]][i,2]),c(3,3),col="firebrick",lwd=lw)
}
k<-2
N<-dim(L_rdat_gs_h1[[k]])[1]
for(i in 1:N){
	lines(c(L_rdat_gs_h1[[k]][i,1],L_rdat_gs_h1[[k]][i,2]),c(2,2),col="firebrick",lwd=lw)
}
k<-2
N<-dim(L_rdat_gus_h2[[k]])[1]
for(i in 1:N){
	lines(c(L_rdat_gus_h2[[k]][i,1],L_rdat_gus_h2[[k]][i,2]),c(1,1),col="firebrick",lwd=lw)
}
abline(v=0,lty=2,lwd=2)
axis(2,at=c(1:4),c("154 G","154 S","Ref G","Ref S"),cex.axis=ca)

title(main="(F) Breakpoint B",cex.main=cm)

plot(c(-1*bnd,bnd),c(.8,4.2),type='n',xlab="Position (bps)",ylab="Reference genome",axes=FALSE,cex.lab=cl)
axis(1,at=seq(-1*bnd,bnd,bnd/2))
k<-3
N<-dim(L_rdat_gs_r1[[k]])[1]
for(i in 1:N){
	lines(c(L_rdat_gs_r1[[k]][i,1],L_rdat_gs_r1[[k]][i,2]),c(4,4),col="firebrick",lwd=lw)
}
k<-3
N<-dim(L_rdat_gus_r1[[k]])[1]
for(i in 1:N){
	lines(c(L_rdat_gus_r1[[k]][i,1],L_rdat_gus_r1[[k]][i,2]),c(3,3),col="firebrick",lwd=lw)
}
k<-3
N<-dim(L_rdat_gs_h1[[k]])[1]
for(i in 1:N){
	lines(c(L_rdat_gs_h1[[k]][i,1],L_rdat_gs_h1[[k]][i,2]),c(2,2),col="firebrick",lwd=lw)
}
k<-3
N<-dim(L_rdat_gus_h2[[k]])[1]
for(i in 1:N){
	lines(c(L_rdat_gus_h2[[k]][i,1],L_rdat_gus_h2[[k]][i,2]),c(1,1),col="firebrick",lwd=lw)
}
abline(v=0,lty=2,lwd=2)

title(main="(F) Breakpoint C",cex.main=cm)
axis(2,at=c(1:4),c("154 G","154 S","Ref G","Ref S"),cex.axis=ca)
dev.off()
