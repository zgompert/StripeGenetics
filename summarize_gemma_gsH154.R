## summarize gemma results, gs genome
library(data.table)
library(scales)

## read in data
gwa_p<-fread("output/o_h154_pattern_gs.assoc.txt",header=TRUE)
gwa_p_no8<-fread("output/o_h154_pattern_no8_gs.assoc.txt",header=TRUE)
gwa_c<-fread("output/o_h154_color_gs.assoc.txt",header=TRUE)
gwa_c_no8<-fread("output/o_h154_color_no8_gs.assoc.txt",header=TRUE)

## read snps
snp<-read.table("snps_gs.txt",header=FALSE)

## signal only on chr8 = scaffold 11 for gs
ch8<-which(snp[,1]==11)
inv_trans_gs<-c(24457103,32202650)
rest_trans_gs<-c(32202951,39030359)
## g933 = stripe h154
coraz_gs<-c(34833964,34837819)

cl<-1.45
pdf("gwa_h154_gs.pdf",height=10,width=8)
par(mfrow=c(4,1))
par(mar=c(4.5,4.5,2.5,1))
plot(snp[ch8,2],-1 * log10(gwa_p_no8$p_lrt[ch8]),pch=19,col=alpha("black",.6),xlab="Position (bp) chromosome 8",ylab="-log10(P)",cex.lab=cl)
polygon(c(inv_trans_gs,rev(inv_trans_gs)),c(-10,-10,100,100),col=alpha("tan",.2),border=NA)
polygon(c(rest_trans_gs,rev(rest_trans_gs)),c(-10,-10,100,100),col=alpha("darkgray",.2),border=NA)
abline(v=coraz_gs,lty=3,col="red")
title("(a) Pattern, GS H154 genome, k without 8")

plot(snp[ch8,2],-1 * log10(gwa_p$p_lrt[ch8]),pch=19,col=alpha("black",.6),xlab="Position (bp) chromosome 8",ylab="-log10(P)",cex.lab=cl)
polygon(c(inv_trans_gs,rev(inv_trans_gs)),c(-10,-10,100,100),col=alpha("tan",.2),border=NA)
polygon(c(rest_trans_gs,rev(rest_trans_gs)),c(-10,-10,100,100),col=alpha("darkgray",.2),border=NA)
abline(v=coraz_gs,lty=3,col="red")
title("(b) Pattern, GS H154 genome, k with 8")

plot(snp[ch8,2],-1 * log10(gwa_c_no8$p_lrt[ch8]),pch=19,col=alpha("black",.6),xlab="Position (bp) chromosome 8",ylab="-log10(P)",cex.lab=cl)
polygon(c(inv_trans_gs,rev(inv_trans_gs)),c(-10,-10,100,100),col=alpha("tan",.2),border=NA)
polygon(c(rest_trans_gs,rev(rest_trans_gs)),c(-10,-10,100,100),col=alpha("darkgray",.2),border=NA)
abline(v=coraz_gs,lty=3,col="red")
title("(c) Color, GS H154 genome, k without 8")

plot(snp[ch8,2],-1 * log10(gwa_c$p_lrt[ch8]),pch=19,col=alpha("black",.6),xlab="Position (bp) chromosome 8",ylab="-log10(P)",cex.lab=cl)
polygon(c(inv_trans_gs,rev(inv_trans_gs)),c(-10,-10,100,100),col=alpha("tan",.2),border=NA)
polygon(c(rest_trans_gs,rev(rest_trans_gs)),c(-10,-10,100,100),col=alpha("darkgray",.2),border=NA)
abline(v=coraz_gs,lty=3,col="red")
title("(d) Color, GS H154 genome, k with 8")
dev.off()
