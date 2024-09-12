library(GENESPACE)

##############################################
## see format_genespace_ch8.sh

###############################################
wd<-"/scratch/general/nfs1/u6000989/TimemaCh8/"
path2mcscanx<-"~/bin/"

# -- initalize the run and QC the inputs
gpar<-init_genespace(wd=wd,path2mcscanx=path2mcscanx)

## need to set this
gpar$shellCalls$orthofinder<-"orthofinder"

# -- accomplish the run
out <- run_genespace(gpar, overwrite = T)


# define region of interest, chr8, t_crist_refug_green_h1, t_crist_refug_stripe_h1
# green = Scaffold_3__1_contigs__length_98311894 
# stripe = Scaffold_11__3_contigs__length_97865747
roi<-data.frame(
		genome=c("t_crist_h154_green_h2","t_crist_h154_stripe_h1",
		"t_crist_refug_green_h1","t_crist_refug_stripe_h1"),
		chr=c("ScrX45T_23_HRSCAF_264",
			"Scaffold_11__2_contigs__length_91821751",
			"Scaffold_3__1_contigs__length_98311894",
			"Scaffold_11__3_contigs__length_97865747"),
		start=c(0,0),end=c(Inf,Inf))

invchr <- data.frame(
  genome = c("t_crist_h154_green_h2","t_crist_h154_stripe_h1"), 
  chr = c("ScrX45T_23_HRSCAF_264","Scaffold_11__2_contigs__length_91821751"))

ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white"))
customPal <- colorRampPalette(
  c("darkorange", "skyblue", "darkblue", "purple", "darkred", "salmon"))

pdf("syn4way_orderA.pdf",width=9,height=7)
ripd <- plot_riparian(
	gsParam = out,
	useRegions = FALSE,
	invertTheseChrs = invchr,
  	refGenome = "t_crist_h154_green_h2",
  	highlightBed = roi,
	genomeIDs = c("t_crist_h154_green_h2","t_crist_h154_stripe_h1",
                "t_crist_refug_green_h1","t_crist_refug_stripe_h1"),
	backgroundColor = NULL)
dev.off()

pdf("syn4way_orderB.pdf",width=9,height=7)
ripd <- plot_riparian(
	gsParam = out,
	useRegions = FALSE,
	invertTheseChrs = invchr,
  	refGenome = "t_crist_h154_green_h2",
  	highlightBed = roi,
	genomeIDs = c("t_crist_h154_green_h2","t_crist_refug_green_h1","t_crist_refug_stripe_h1","t_crist_h154_stripe_h1"),
	backgroundColor = NULL)
dev.off()

save(list=ls(),file="gs_4way.rdat")

hits<-read_allBlast(filepath="/scratch/general/nfs1/u6000989/TimemaCh8/syntenicHits/t_crist_refug_green_h1_vs_t_crist_refug_stripe_h1.allBlast.txt.gz")
ggdotplot(hits = hits, type = "syntenic", verbose = FALSE)
hits<-read_allBlast(filepath="/scratch/general/nfs1/u6000989/TimemaCh8/syntenicHits/t_crist_h154_green_h2_vs_t_crist_h154_stripe_h1.allBlast.txt.gz")
ggdotplot(hits = hits, type = "syntenic", verbose = FALSE)

qreturn <- query_hits(gsParam = out, bed = roi, synOnly = TRUE)
xvx <- subset(qreturn[["t_crist_h154_green_h2, ScrX45T_23_HRSCAF_264: 0-Inf"]], genome2 == "t_crist_h154_stripe_h1")
gghits(xvx, useOrder = F)

xvx_h154 <- subset(qreturn[["t_crist_h154_green_h2, ScrX45T_23_HRSCAF_264: 0-Inf"]], genome2 == "t_crist_h154_stripe_h1")
xvx_refugio <- subset(qreturn[["t_crist_refug_green_h1, Scaffold_3__1_contigs__length_98311894: 0-Inf"]], genome2 == "t_crist_refug_stripe_h1")
xvx_stripe <- subset(qreturn[["t_crist_refug_stripe_h1, Scaffold_11__3_contigs__length_97865747: 0-Inf"]], genome2 == "t_crist_h154_stripe_h1")
xvx_green<- subset(qreturn[["t_crist_refug_green_h1, Scaffold_3__1_contigs__length_98311894: 0-Inf"]], genome2 == "t_crist_h154_green_h2")

## pro-corazonin
## g900 = stripe refugio
## g4303 = green refugio
## g933 = stripe h154
## g3288 = green h154

corz_refugio<-which(xvx_refugio$id1=="g4303")
corz_h154<-which(xvx_h154$id1=="g3288")
corz_stripe<-which(xvx_stripe$id1=="g900")
corz_green<-which(xvx_green$id1=="g4303")

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
