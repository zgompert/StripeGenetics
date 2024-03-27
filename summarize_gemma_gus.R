## summarize gemma results, gus genome
library(data.table)
library(scales)

## read in data
gwa_p<-fread("output/o_ref_pattern_gus.assoc.txt",header=TRUE)
gwa_p_no8<-fread("output/o_ref_pattern_no8_gus.assoc.txt",header=TRUE)
gwa_c<-fread("output/o_ref_color_gus.assoc.txt",header=TRUE)
gwa_c_no8<-fread("output/o_ref_color_no8_gus.assoc.txt",header=TRUE)

## read snps
snp<-read.table("snps_gus.txt",header=FALSE)

## signal only on chr8 = scaffold 3 for gs
ch8<-which(snp[,1]==3)


cl<-1.45
pdf("gwa_refugio_gus.pdf",height=10,width=8)
par(mfrow=c(4,1))
par(mar=c(4.5,4.5,2.5,1))
plot(snp[ch8,2],-1 * log10(gwa_p_no8$p_lrt[ch8]),pch=19,col=alpha("black",.6),xlab="Position (bp) chromosome 8",ylab="-log10(P)",cex.lab=cl)
title("(a) Pattern, GUS Refugio genome, k without 8")
plot(snp[ch8,2],-1 * log10(gwa_p$p_lrt[ch8]),pch=19,col=alpha("black",.6),xlab="Position (bp) chromosome 8",ylab="-log10(P)",cex.lab=cl)
title("(b) Pattern, GUS Refugio genome, k with 8")
plot(snp[ch8,2],-1 * log10(gwa_c_no8$p_lrt[ch8]),pch=19,col=alpha("black",.6),xlab="Position (bp) chromosome 8",ylab="-log10(P)",cex.lab=cl)
title("(c) Color, GUS Refugio genome, k without 8")
plot(snp[ch8,2],-1 * log10(gwa_c$p_lrt[ch8]),pch=19,col=alpha("black",.6),xlab="Position (bp) chromosome 8",ylab="-log10(P)",cex.lab=cl)
title("(d) Color, GUS Refugio genome, k with 8")
dev.off()
