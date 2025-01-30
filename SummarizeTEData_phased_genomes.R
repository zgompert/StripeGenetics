## summarize TE information for entire genome, comparison with SV 

gff_h154_gs_h1<-read.table("t_crist_gs_hap_cen4119/HiRise/Hap1/final_assembly.fasta.out.gff", header=FALSE)
gff_h154_gs_h2<-read.table("t_crist_gs_hap_cen4119/HiRise/Hap2/final_assembly.fasta.out.gff", header=FALSE)
gff_h154_gus_h1<-read.table("t_crist_gus_hap_cen4280/HiRise/Hap1/ojincantatabio-cen4280-hap1-mb-hirise-ig5ps__01-30-2024__hic_output.fasta.out.gff",header=FALSE)
gff_h154_gus_h2<-read.table("t_crist_gus_hap_cen4280/HiRise/Hap2/ojincantatabio-cen4280-hap2-mb-hirise-i2xb7__01-30-2024__hic_output.fasta.out.gff",header=FALSE)

gff_refug_gus_h1<-read.table("t_crist_refug_green/HiRise/hap1/ojincantatabio-cen4120-hap1-mb-hirise-wlbll__08-15-2023__final_assembly.fasta.out.gff", header=FALSE)
gff_refug_gus_h2<-read.table("t_crist_refug_green/HiRise/hap2/ojincantatabio-cen4120-hap2-mb-hirise-bn0ko__08-15-2023__final_assembly.fasta.out.gff", header=FALSE)
gff_refug_gs_h1<-read.table("t_crist_refug_stripe/HiRise/hap1/ojincantatabio-cen4122-hap1-mb-hirise-g4hzf__08-10-2023__final_assembly.fasta.out.gff", header=FALSE)
gff_refug_gs_h2<-read.table("t_crist_refug_stripe/HiRise/hap2/ojincantatabio-cen4122-hap2-mb-hirise-14fv0__08-10-2023__final_assembly.fasta.out.gff", header=FALSE)

## rename but no actual subsetting
loc_h154_gus<-gff_h154_gus_h2
loc_h154_gs<-gff_h154_gs_h1
loc_refug_gus<-gff_refug_gs_h1
loc_refug_gs<-gff_refug_gus_h1
########### H154 #####################


## totals
sum(loc_h154_gus$V5-loc_h154_gus$V4)
## 713,704,933 bps 
sum(loc_h154_gus$V5-loc_h154_gus$V4)/1244673105
#[1] 0.5734075
sum(loc_h154_gs$V5-loc_h154_gs$V4)
## 711,246,201 bps  
sum(loc_h154_gs$V5-loc_h154_gs$V4)/1248345166
##  0.5697512

## but what about types
## indexes for simple repeats
gs_h154_sim_rep<-grep(pattern=")n",x=loc_h154_gs$V10)
gus_h154_sim_rep<-grep(pattern=")n",x=loc_h154_gus$V10)
sum(loc_h154_gus$V5[gus_h154_sim_rep]-loc_h154_gus$V4[gus_h154_sim_rep])
sum(loc_h154_gs$V5[gs_h154_sim_rep]-loc_h154_gs$V4[gs_h154_sim_rep]
loc_h154_gs$V10[gs_h154_sim_rep]<-"SimpleRep"
loc_h154_gus$V10[gus_h154_sim_rep]<-"SimpleRep"

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

pdf("svSFX_tesums_genome.pdf",width=9,height=9)
par(mfrow=c(2,2))
par(mar=c(4.5,8,3,1))
cl<-1.4;ca<-1;cm<-1.4

dotchart(gus_type_cnts,cex=.5,,pch=19,xlab="Base pairs",cex.lab=cl,cex.axis=ca)
title(main="(C) Hwy 154 green",cex.main=cm)
dotchart(gs_type_cnts, cex = 0.5,pch=19,xlab="Base pairs",cex.lab=cl,cex.axis=ca) 
title(main="(D) Hwy 154 striped",cex.main=cm)
dev.off()
########### refugio #####################

## totals
sum(loc_refug_gus$V5-loc_refug_gus$V4)
# 725,898,450 bps 
sum(loc_refug_gus$V5-loc_refug_gus$V4)/1363418462
#[1] 0.5324106
sum(loc_refug_gs$V5-loc_refug_gs$V4)
# 733,214,197 bps
sum(loc_refug_gs$V5-loc_refug_gs$V4)/1248345166
#[1] 0.5873489


## but what about types
## indexes for simple repeats
gs_refug_sim_rep<-grep(pattern=")n",x=loc_refug_gs$V10)
gus_refug_sim_rep<-grep(pattern=")n",x=loc_refug_gus$V10)
sum(loc_refug_gus$V5[gus_refug_sim_rep]-loc_refug_gus$V4[gus_refug_sim_rep])
sum(loc_refug_gs$V5[gs_refug_sim_rep]-loc_refug_gs$V4[gs_refug_sim_rep])
loc_refug_gs$V10[gs_refug_sim_rep]<-"SimpleRep"
loc_refug_gus$V10[gus_refug_sim_rep]<-"SimpleRep"

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

pdf("svSFX_tesumsR_genome.pdf",width=9,height=9)
par(mfrow=c(2,2))
par(mar=c(4.5,8,3,1))
cl<-1.4;ca<-1;cm<-1.4

dotchart(gus_type_cnts,cex=.5,,pch=19,xlab="Base pairs",cex.lab=cl,cex.axis=ca)
title(main="(A) Refugio green",cex.main=cm)
dotchart(gs_type_cnts, cex = 0.5,pch=19,xlab="Base pairs",cex.lab=cl,cex.axis=ca)
title(main="(B) Refugio striped",cex.main=cm)
dev.off()

