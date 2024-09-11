## summarize gemma results, gus genome
library(data.table)
library(scales)

## read in data
gwa_p<-fread("output/o_h154_pattern_gus.assoc.txt",header=TRUE)
gwa_p_no8<-fread("output/o_h154_pattern_no8_gus.assoc.txt",header=TRUE)
gwa_c<-fread("output/o_h154_color_gus.assoc.txt",header=TRUE)
gwa_c_no8<-fread("output/o_h154_color_no8_gus.assoc.txt",header=TRUE)

## read snps
snp<-read.table("snps_gus.txt",header=FALSE)

## signal only on chr8 = scaffold 23 for gus
ch8<-which(snp[,1]==23)
inv_trans_gus<-c(33112541,44121870)
rest_trans_gus<-c(24803527,32951302)
coraz_gus<-c(32607929,32612697)

cl<-1.45
pdf("gwa_h154_gus.pdf",height=10,width=8)
par(mfrow=c(4,1))
par(mar=c(4.5,4.5,2.5,1))
plot(snp[ch8,2],-1 * log10(gwa_p_no8$p_lrt[ch8]),pch=19,col=alpha("black",.6),xlab="Position (bp) chromosome 8",ylab="-log10(P)",cex.lab=cl)
polygon(c(inv_trans_gus,rev(inv_trans_gus)),c(-10,-10,200,200),col=alpha("tan",.2),border=NA)
polygon(c(rest_trans_gus,rev(rest_trans_gus)),c(-10,-10,200,200),col=alpha("darkgray",.2),border=NA)
abline(v=coraz_gus,lty=3,col="red")
title("(a) Pattern, GUS H154 genome, k without 8")
plot(snp[ch8,2],-1 * log10(gwa_p$p_lrt[ch8]),pch=19,col=alpha("black",.6),xlab="Position (bp) chromosome 8",ylab="-log10(P)",cex.lab=cl)
polygon(c(inv_trans_gus,rev(inv_trans_gus)),c(-10,-10,200,200),col=alpha("tan",.2),border=NA)
polygon(c(rest_trans_gus,rev(rest_trans_gus)),c(-10,-10,200,200),col=alpha("darkgray",.2),border=NA)
abline(v=coraz_gus,lty=3,col="red")
title("(b) Pattern, GUS H154 genome, k with 8")
plot(snp[ch8,2],-1 * log10(gwa_c_no8$p_lrt[ch8]),pch=19,col=alpha("black",.6),xlab="Position (bp) chromosome 8",ylab="-log10(P)",cex.lab=cl)
polygon(c(inv_trans_gus,rev(inv_trans_gus)),c(-10,-10,100,100),col=alpha("tan",.2),border=NA)
polygon(c(rest_trans_gus,rev(rest_trans_gus)),c(-10,-10,100,100),col=alpha("darkgray",.2),border=NA)
abline(v=coraz_gus,lty=3,col="red")
title("(c) Color, GUS H154 genome, k without 8")
plot(snp[ch8,2],-1 * log10(gwa_c$p_lrt[ch8]),pch=19,col=alpha("black",.6),xlab="Position (bp) chromosome 8",ylab="-log10(P)",cex.lab=cl)
polygon(c(inv_trans_gus,rev(inv_trans_gus)),c(-10,-10,100,100),col=alpha("tan",.2),border=NA)
polygon(c(rest_trans_gus,rev(rest_trans_gus)),c(-10,-10,100,100),col=alpha("darkgray",.2),border=NA)
abline(v=coraz_gus,lty=3,col="red")
title("(d) Color, GUS H154 genome, k with 8")
dev.off()
