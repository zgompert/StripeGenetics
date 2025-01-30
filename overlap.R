## get unique gene names
ff<-list.files(pattern="uniq_SV_genes")
funL<-vector("list",4)
for(k in 1:4){
	funL[[k]]<-read.table(ff[k],header=FALSE,sep="\t")
}


sum(funL[[1]][,1] %in% funL[[2]][,1] & funL[[1]][,1] %in% funL[[3]][,1] & funL[[1]][,1] %in% funL[[4]][,1])
## 23
aa<-which(funL[[1]][,1] %in% funL[[2]][,1] & funL[[1]][,1] %in% funL[[3]][,1] & funL[[1]][,1] %in% funL[[4]][,1])
funL[[1]][aa,1]
