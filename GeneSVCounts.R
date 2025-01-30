library(GENESPACE)
load("/scratch/general/nfs1/u6000989/TimemaCh8_4/gspace.rdat")

## from 4way comparison

gid2_h154<-as.numeric(gsub(pattern="g",replacement="",x=xvx_h154$id2))
gid1_h154<-as.numeric(gsub(pattern="g",replacement="",x=xvx_h154$id1))
length(which(gid2_h154 >=887 & gid2_h154 <=936 & gid1_h154 >= 3262 & gid1_h154 <= 3345))
h154_match<-which(gid2_h154 >=887 & gid2_h154 <=936 & gid1_h154 >= 3262 & gid1_h154 <= 3345)
## 42 in both
(936-887+1-length(unique(gid2_h154[h154_match]))) + (3345-3262+1-length(unique(gid1_h154[h154_match]))) + max(length(unique(gid2_h154[h154_match])),length(unique(gid1_h154[h154_match])))
## 97

gid2_ref<-as.numeric(gsub(pattern="g",replacement="",x=xvx_refugio$id2))
gid1_ref<-as.numeric(gsub(pattern="g",replacement="",x=xvx_refugio$id1))
length(which(gid2_ref >= 681 & gid2_ref <= 900 & gid1_ref >= 4059 & gid1_ref <= 4303))
ref_match<-which(gid2_ref >= 681 & gid2_ref <= 900 & gid1_ref >= 4059 & gid1_ref <= 4303)
## 212 in both
## so total is number in each plus number combined:
(900-681+1-length(unique(gid2_ref[ref_match]))) + (4303-4059+1-length(unique(gid2_ref[ref_match]))) + max(length(unique(gid2_ref[ref_match])),length(unique(gid1_ref[ref_match])))
## 299

##### USED ABOVE, PROTEIN MATCHES FOR OVERLAP #####


gid2_str<-as.numeric(gsub(pattern="g",replacement="",x=xvx_stripe$id2)) ## h154
gid1_str<-as.numeric(gsub(pattern="g",replacement="",x=xvx_stripe$id1)) ## refugio
length(which(gid2_str >= 887 & gid2_str <= 936 & gid1_str >= 681 & gid1_str <= 900))
str_match<-which(gid2_str >= 887 & gid2_str <= 936 & gid1_str >= 681 & gid1_str <= 900)
## 26 in both
sort(gid2_str[str_match])
## range 2 = 888-936 = 49
sort(gid1_str[str_match])
## range 1 = 683:698 + 870-900 = 38
49+38-26
## 61 total unique


gid2_grn<-as.numeric(gsub(pattern="g",replacement="",x=xvx_green$id2)) ## h154
gid1_grn<-as.numeric(gsub(pattern="g",replacement="",x=xvx_green$id1)) ## refugio
length(which(gid2_grn >= 3262 & gid2_grn <= 3345 & gid1_grn >= 4059 & gid1_grn <= 4303))
grn_match<-which(gid2_grn >= 3262 & gid2_grn <= 3345 & gid1_grn >= 4059 & gid1_grn <= 4303)
## 38 in both, 4068-4156 + 4303 for grn, 3289-3338 + 3288
sort(gid2_grn[grn_match])
## 3288:3345 = 58
sort(gid1_grn[grn_match])
## 4068:4124, 4156, 4303 = 59
58+59-38
## 79 total unique


#### another way, only counts actual matches

## start with 38, add any of the 26 that are not already accounted for
## stripe matching that show up in h154 match
mm<-gid2_str[str_match][(gid2_str[str_match] %in% gid2_h154[h154_match])]
## green number for these 
grn_mm<-gid1_h154[which(gid2_h154 %in% mm)]
## stripe matching that show up in ref match
mm<-gid1_str[str_match][(gid1_str[str_match] %in% gid2_ref[ref_match])]
## green number for these 
grn_mm2<-gid1_ref[which(gid2_ref %in% mm)]

sum((gid1_grn[grn_match] %in% grn_mm2)  | (gid2_grn[grn_match] %in% grn_mm))
#[1] 25
## suggest 25 of 26 are redundant, so total is 38+1 = 39.

## trying something based on the pangenome

roi<-data.frame(
                genome=c("t_crist_h154_green_h2",
                         "t_crist_h154_stripe_h1",
                "t_crist_refug_green_h1","t_crist_refug_stripe_h1"),
                chr=c("ScrX45T_23_HRSCAF_264",
                        "Scaffold_11__2_contigs__length_91821751",
                        "Scaffold_3__1_contigs__length_98311894",
                        "Scaffold_11__3_contigs__length_97865747"),
                start=c(24803527,24457103,22220178,22442098),end=c(44121870,39030359,65829835,65729704))

pan4<-query_pangenes(gsParam=out,bed=roi,showArrayMem = TRUE)

sort(table(c(pan4[[1]]$pgID,pan4[[2]]$pgID,pan4[[3]]$pgID,pan4[[4]]$pgID)))

