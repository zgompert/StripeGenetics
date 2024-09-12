library(GENESPACE)

##############################################
## see format_genespace_ch8.sh

###############################################
wd<-"/scratch/general/nfs1/u6000989/TimemaCh8_8/"
path2mcscanx<-"~/bin/"

# -- initalize the run and QC the inputs
gpar<-init_genespace(wd=wd,path2mcscanx=path2mcscanx)

## need to set this
gpar$shellCalls$orthofinder<-"orthofinder"

# -- accomplish the run
out <- run_genespace(gpar, overwrite = T)
#`load('/scratch/general/nfs1/u6000989/TimemaCh8_8//results/gsParams.rda',


# define region of interest, chr8, t_crist_refug_green_h1, t_crist_refug_stripe_h1
# green = Scaffold_3__1_contigs__length_98311894 
# stripe = Scaffold_11__3_contigs__length_97865747
roi<-data.frame(
		genome=c("t_crist_h154_green_h1","t_crist_h154_green_h2",
			 "t_crist_h154_stripe_h1","t_crist_h154_stripe_h2",
		"t_crist_refug_green_h1","t_crist_refug_green_h2",
		"t_crist_refug_stripe_h1","t_crist_refug_stripe_h2"),
		chr=c("ScIzd6S_7_HRSCAF_90","ScrX45T_23_HRSCAF_264",
			"Scaffold_11__2_contigs__length_91821751",
			"Scaffold_4__1_contigs__length_97222829",
			#"Scaffold_6__1_contigs__length_78844258",## is this right?? NO
			"Scaffold_3__1_contigs__length_98311894",
			"Scaffold_3__1_contigs__length_97218928",
			"Scaffold_11__3_contigs__length_97865747",
			"Scaffold_9__2_contigs__length_97179199"),
		start=c(0,0,0,0,0,0,0,0),end=c(Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf))

invchr <- data.frame(
  genome = c("t_crist_h154_green_h1","t_crist_h154_green_h2","t_crist_h154_stripe_h1",
	     "t_crist_refug_stripe_h2","t_crist_h154_stripe_h2"), 
  chr = c("ScIzd6S_7_HRSCAF_90","ScrX45T_23_HRSCAF_264","Scaffold_11__2_contigs__length_91821751",
	  "Scaffold_9__2_contigs__length_97179199","Scaffold_4__1_contigs__length_97222829"))


ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white"))
customPal <- colorRampPalette(
  c("skyblue"))
  #c("darkorange", "skyblue", "darkblue", "purple", "darkred", "salmon"))

pdf("syn8way.pdf",width=9,height=7)
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
	genomeIDs = c("t_crist_h154_green_h1","t_crist_h154_green_h2","t_crist_h154_stripe_h1","t_crist_h154_stripe_h2",
                "t_crist_refug_green_h1","t_crist_refug_green_h2","t_crist_refug_stripe_h1","t_crist_refug_stripe_h2"),
	backgroundColor = NULL)
dev.off()

save(list=ls(),file="syn8way.rdat")
