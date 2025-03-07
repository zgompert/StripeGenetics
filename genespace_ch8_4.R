library(GENESPACE)

#############################################
## see format_genespace_ch8.sh

###############################################
wd<-"/scratch/general/nfs1/u6000989/TimemaCh8_4/"
path2mcscanx<-"~/bin/"

# -- initalize the run and QC the inputs
gpar<-init_genespace(wd=wd,path2mcscanx=path2mcscanx)

## need to set this
gpar$shellCalls$orthofinder<-"orthofinder"

# -- accomplish the run
out <- run_genespace(gpar, overwrite = T)
#`load('/scratch/general/nfs1/u6000989/TimemaCh8_5//results/gsParams.rda',


# define region of interest, chr8, t_crist_refug_green_h1, t_crist_refug_stripe_h1
# green = Scaffold_3__1_contigs__length_98311894 
# stripe = Scaffold_11__3_contigs__length_97865747
roi<-data.frame(
		genome=c("t_crist_h154_green_h2",
			 "t_crist_h154_stripe_h1",
		"t_crist_refug_green_h1","t_crist_refug_stripe_h1"),
		chr=c("ScrX45T_23_HRSCAF_264",
			"Scaffold_11__2_contigs__length_91821751",
			"Scaffold_3__1_contigs__length_98311894",
			"Scaffold_11__3_contigs__length_97865747"),
		start=c(0,0,0,0),end=c(Inf,Inf,Inf,Inf))

invchr <- data.frame(
  genome = c("t_crist_h154_green_h2","t_crist_h154_stripe_h1"), 
  chr = c("ScrX45T_23_HRSCAF_264","Scaffold_11__2_contigs__length_91821751"))


ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white"))
customPal <- colorRampPalette(
  #c("salmon"))
  c("darkorange", "skyblue", "darkblue", "purple", "darkred", "salmon"))

pdf("syn4way.pdf",width=9,height=5.6)
ripd <- plot_riparian(
	gsParam = out,
	palette = customPal,
        highlightBed = roi,
	braidAlpha = .3,
	useOrder=TRUE,
	chrFill = "lightgrey",
        addThemes = ggthemes,
	useRegions = FALSE,
	invertTheseChrs = invchr,
  	refGenome = "t_crist_h154_green_h2",
	genomeIDs = c("t_crist_h154_stripe_h1","t_crist_h154_green_h2",
                "t_crist_refug_green_h1","t_crist_refug_stripe_h1"),
	#genomeIDs = c("t_crist_h154_green_h2","t_crist_h154_stripe_h1",
        #        "t_crist_refug_green_h1","t_crist_refug_stripe_h1"),
	backgroundColor = NULL)
dev.off()

#hits<-read_allBlast(filepath="/scratch/general/nfs1/u6000989/TimemaCh8/syntenicHits/t_crist_refug_green_h1_vs_t_crist_refug_stripe_h1.allBlast.txt.gz")
#ggdotplot(hits = hits, type = "syntenic", verbose = FALSE)
#hits<-read_allBlast(filepath="/scratch/general/nfs1/u6000989/TimemaCh8/syntenicHits/t_crist_h154_green_h2_vs_t_crist_h154_stripe_h1.allBlast.txt.gz")
#ggdotplot(hits = hits, type = "syntenic", verbose = FALSE)

qreturn <- query_hits(gsParam = out, bed = roi, synOnly = TRUE)
xvx <- subset(qreturn[["t_crist_h154_green_h2, ScrX45T_23_HRSCAF_264: 0-Inf"]], genome2 == "t_crist_h154_stripe_h1")
gghits(xvx, useOrder = F)

xvx_h154 <- subset(qreturn[["t_crist_h154_green_h2, ScrX45T_23_HRSCAF_264: 0-Inf"]], genome2 == "t_crist_h154_stripe_h1")
xvx_refugio <- subset(qreturn[["t_crist_refug_green_h1, Scaffold_3__1_contigs__length_98311894: 0-Inf"]], genome2 == "t_crist_refug_stripe_h1")
xvx_stripe <- subset(qreturn[["t_crist_refug_stripe_h1, Scaffold_11__3_contigs__length_97865747: 0-Inf"]], genome2 == "t_crist_h154_stripe_h1")
xvx_green<- subset(qreturn[["t_crist_refug_green_h1, Scaffold_3__1_contigs__length_98311894: 0-Inf"]], genome2 == "t_crist_h154_green_h2")

## notes are from manual inspection with genomeview browser 

## pro-corazonin
## g900 = stripe refugio = 2 exons, 12,005 bp intron... but only 4207 in stripe hap 2 
## so not a consistent difference with green... all others are in fact ~4000 bp, always appears single copy
## g4303 = green refugio = 2 exons, 4201 bp intron
## g933 = stripe h154 = 2 exons, 3469 bp intron
## g3288 = green h154 = 2 exons, 4382 bp intron


corz_refugio<-which(xvx_refugio$id1=="g4303")
corz_h154<-which(xvx_h154$id1=="g3288")
corz_stripe<-which(xvx_stripe$id1=="g900")
corz_green<-which(xvx_green$id1=="g4303")

## Ecdysteroid kinase-like family
## g890 and 691 = stripe refugio
## both 1 exon
## g4118 and 4081 = green refugio ... Ecdysone receptor 4201 but only in this genome
## both 1 exon, receptor has 8 exons
## g919 and 893 = stripe h154
## both 1 exon
## g3302 = green h154
## 1 exon
ecdys_refugio<-which((xvx_refugio$id1=="g4118") | (xvx_refugio$id1=="g4081"))
ecdys_h154<-which((xvx_h154$id1=="g3302") | (xvx_h154$id1=="g3337"))
ecdys_stripe<-which((xvx_stripe$id1=="g890") | (xvx_stripe$id1=="g691"))
ecdys_green<-which((xvx_green$id1=="g4118") | (xvx_green$id1=="g4081"))

pdf("geneDotPlots.pdf",width=9,height=9)
par(mfrow=c(2,2))
par(mar=c(4.5,5.5,2.5,1.5))
plot(xvx_refugio$ord1,xvx_refugio$ord2,pch=19,col=alpha("darkgray",.7),xlab="Gene order GUS R",ylab="Gene order GS R")
points(xvx_refugio$ord1[corz_refugio],xvx_refugio$ord2[corz_refugio],pch=19,col="firebrick",cex=1.5)
title(main="(A) GUS x GS Refugio")

## flip for consistency
max1<-max(xvx_h154$ord1);max2<-max(xvx_h154$ord2)
plot(max1-xvx_h154$ord1,max2-xvx_h154$ord2,pch=19,col=alpha("darkgray",.7),xlab="Gene order GUS H",ylab="Gene order GS H")
points(max1-xvx_h154$ord1[corz_h154],max2-xvx_h154$ord2[corz_h154],pch=19,col="firebrick",cex=1.5)
title(main="(B) GUS x GS H154")

max3<-max(xvx_stripe$ord2)
plot(xvx_stripe$ord1,max3-xvx_stripe$ord2,pch=19,col=alpha("darkgray",.7),xlab="Gene order GS R",ylab="Gene order GS H")
points(xvx_stripe$ord1[corz_stripe],max3-xvx_stripe$ord2[corz_stripe],pch=19,col="firebrick",cex=1.5)
title(main="(C) GS Refugio x H154")

max4<-max(xvx_green$ord2)
plot(xvx_green$ord1,max4-xvx_green$ord2,pch=19,col=alpha("darkgray",.7),xlab="Gene order GUS R",ylab="Gene order GUS H")
points(xvx_green$ord1[corz_green],max4-xvx_green$ord2[corz_green],pch=19,col="firebrick",cex=1.5)
title(main="(D) GUS Refugio x H154")
dev.off()

cl<-1.45
ca<-1.0
cm<-1.3
pdf("svF4_genes.pdf",width=9,height=6)
par(mfrow=c(2,3))
par(mar=c(4.5,5.5,2.5,1.5))
plot(xvx_refugio$ord2,xvx_refugio$ord1,pch=19,col=alpha("darkgray",.7),xlab="Gene order stripe",ylab="Gene order green",cex.axis=ca,cex.lab=cl)
points(xvx_refugio$ord2[corz_refugio],xvx_refugio$ord1[corz_refugio],pch=19,col="firebrick",cex=1.5)
points(xvx_refugio$ord2[ecdys_refugio],xvx_refugio$ord1[ecdys_refugio],pch=19,col="cadetblue",cex=1.5)
title(main="(B) Refugio, stripe vs green",cex.main=cm)
legend(3900,790,c("Pro-corazonin","Ecdysteroid kinase"),pch=19,bty='n',col=c("firebrick","cadetblue"))

## flip for consistency
max1<-max(xvx_h154$ord1);max2<-max(xvx_h154$ord2)
plot(max2-xvx_h154$ord2,max1-xvx_h154$ord1,pch=19,col=alpha("darkgray",.7),xlab="Gene order stripe",ylab="Gene order green",cex.axis=ca,cex.lab=cl)
points(max2-xvx_h154$ord2[ecdys_h154],max1-xvx_h154$ord1[ecdys_h154],pch=19,col="cadetblue",cex=1.5)
points(max2-xvx_h154$ord2[corz_h154],max1-xvx_h154$ord1[corz_h154],pch=19,col="firebrick",cex=1.5)
title(main="(C) Hwy 154, stripe vs green",cex.main=cm)

max4<-max(xvx_green$ord2)
plot(xvx_green$ord1,max4-xvx_green$ord2,pch=19,col=alpha("darkgray",.7),xlab="Gene order Refugio",ylab="Gene order Hwy 154",cex.axis=ca,cex.lab=cl)
points(xvx_green$ord1[corz_green],max4-xvx_green$ord2[corz_green],pch=19,col="firebrick",cex=1.5)
points(xvx_green$ord1[ecdys_green],max4-xvx_green$ord2[ecdys_green],pch=19,col="cadetblue",cex=1.5)
title(main="(D) Green, Refugio vs Hwy 154",cex.main=cm)
dev.off()

save(list=ls(),file="syn4way.rdat")
