#################
# WisAsp data

library(rNeighborGWAS)
library(gaston)
library(tidyverse)
library(patchwork)
library(RAINBOWR)
source("covnei_lmm.R")

setwd("./aspen")

###############
# reshape input data

# PLINK genotypes from Barker et al. 2019 Mol Ecol
g = gaston::read.bed.matrix(basename="WisAsp")

table(substr(g@ped[,"id"],12,14))
table(substr(g@ped[,"id"],17,17))
table(substr(g@ped[,"id"],18,19))

# anntation from Riehl et al. 2022 ECE
ann = read.csv("Riehl2023_ece310541-sup-0002-files2_sheet2.csv")
ann = ann[is.na(ann$Exclusion.Code),]

# trait data from Barker et al. 2018 PLoS ONE
y2014 = read.csv("WisAsp_Data2014.csv",header=TRUE)
y2014$ID = gsub(" ", "_",y2014$ID)
genetID = y2014[(y2014$ID %in% ann$Correct.Sequenced.Ramet),"GENO"]

gIDs = rep(NA,nrow(y2014))
for(i in 1:nrow(y2014)){
  focal = which(ann$Correct.Sequenced.Ramet==y2014$ID[i])
  print(focal)
  if(length(focal)==1){
    gIDs[which(y2014$GENO==y2014$GENO[i])] = ann$ID.RapidGenomics[focal]
  }
}
y2014 = data.frame(gIDs,y2014)
y2014 = y2014[(y2014$gIDs %in% gsub(".bam","",g@ped$id)),]

gwhich = c()
for(i in y2014$gIDs) {
  focal = which(gsub(".bam","",g@ped$id)==i)
  gwhich = c(gwhich,focal)
}

g = g[gwhich,]

gsub(".bam","",g@ped$id) == y2014$gIDs
y2014$ROW = as.numeric(factor(y2014$ROW))

map = data.frame(x=y2014$ROW,y=y2014$POS)

g = rNeighborGWAS::gaston2neiGWAS(g)
AF = apply(g$geno,2,mean)
g$geno = g$geno[,which(AF < 0.95 & AF > 0.05)]
g$gmap = g$gmap[which(AF < 0.95 & AF > 0.05),]

######
# PVE including covariance
y2014[,8:37] = log(y2014[,8:37]+0.5)
traits = colnames(y2014[38:ncol(y2014)])
EM3y2014res12 = c()
for(trait in traits){ print(trait)

  pos = g$gmap
  geno = g$geno
  pheno = y2014
  smap = map
  
  g_nei = nei_coval(geno,smap,scale=1+0.01,grouping=pheno$BLK,n_core=1L)
  # calc. K
  selfCross = tcrossprod(geno)
  neiCross = tcrossprod(g_nei)
  
  q = ncol(geno)
  selfK = (q+selfCross)/(2*(q-1))
  neiK = neiCross/(q-1)
  
  y = as.matrix(pheno[,trait])
  ind_names = paste0("ind",1:length(y))
  rownames(y) = ind_names
  Z = diag(length(y))
  rownames(Z) = ind_names
  colnames(Z) = ind_names
  rownames(selfK) = ind_names
  colnames(selfK) = ind_names
  rownames(neiK) = ind_names
  colnames(neiK) = ind_names
  ZETA1 <- list(A = list(Z = Z, K = selfK),
                Int = list(Z = Z, K = neiK))
  
  K12 = tcrossprod(geno,g_nei) / (q-1)
  K21 = t(K12)
  rownames(K12) = ind_names
  colnames(K12) = ind_names
  rownames(K21) = ind_names
  colnames(K21) = ind_names
  
  covList1 <- rep(list(rep(list(NULL), 2)), 2)
  covList1[[1]][[2]] <- K12
  covList1[[2]][[1]] <- K21
  
  X0 = model.matrix(~factor(pheno$BLK))
  rownames(X0) = ind_names
  
  EM3.res11 <- EM3.cpp(y = y,
                       X0 = X0,
                       ZETA = ZETA1)
  w = EM3.res11$weights
  if(sum(w==TRUE)==1) {
    w[w==1] = 0.99
    w[w==0] = 0.01
  }
  EM3.res12 <- EM3.cov(y = y,
                       X0 = X0,
                       ZETA = ZETA1,
                       covList = covList1,
                       optimizeWeights = TRUE,
                       parInitForWeights = w,
                       forceApproxK = TRUE,
                       nIterOptimization=10^5, traceInside=1)
  
  EM3y2014res11 = c(EM3y2014res11,list(EM3.res11))
  EM3y2014res12 = c(EM3y2014res12,list(EM3.res12))
}
names(EM3y2014res11) = traits
names(EM3y2014res12) = traits

saveRDS(EM3y2014res11,file="../output/aspen/aspen_2014_EM3res11_s1.rds")
saveRDS(EM3y2014res12,file="../output/aspen/aspen_2014_EM3res12_s1.rds")

# show results for Table S1
for(i in 1:length(traits)) print(c(traits[i],EM3y2014res11[[i]]$Vu,EM3y2014res11[[i]]$Ve,EM3y2014res11[[i]]$weights,EM3y2014res11[[i]]$LL))
for(i in 1:length(traits)) print(c(traits[i],EM3y2014res12[[i]]$Vu,EM3y2014res12[[i]]$Ve,EM3y2014res12[[i]]$weights,EM3y2014res12[[i]]$rhosMat[1,2],EM3y2014res12[[i]]$LL))

EM3y2014res11 = readRDS(file="../output/aspen/aspen_2014_EM3res11_s1.rds")
EM3y2014res12 = readRDS(file="../output/aspen/aspen_2014_EM3res12_s1.rds")

# plot PVEs
sum2014 = c()
for(i in 1:length(EM3y2014res12)) sum2014 = rbind(sum2014,c(EM3y2014res12[[i]]$Vu,EM3y2014res12[[i]]$Ve,EM3y2014res12[[i]]$weights,EM3y2014res12[[i]]$rhosMat[1,2],EM3y2014res12[[i]]$LL))
sum2014 = as.data.frame(sum2014)
rownames(sum2014) = names(EM3y2014res12)
colnames(sum2014) = c("Vu","Ve","w1","w2","rho12","LL")

h2 = sum2014$Vu / (sum2014$Vu + sum2014$Ve)
rhow = h2*sum2014$rho12*sqrt(sum2014$w1*sum2014$w2)

pve2014 = data.frame(trait=factor(rownames(sum2014),levels=rownames(sum2014)),DGE=h2*sum2014$w1,IGE=h2*sum2014$w2,cov=h2*sum2014$rho12*sqrt(sum2014$w1*sum2014$w2))
p2014 = pve2014 %>%
  gather(c("DGE","IGE","cov"),key="factor",value="PVE") %>%
  ggplot(aes(x=trait, y=PVE, fill=factor)) + ylim(-0.25,1) + labs(subtitle="2014") +
  geom_bar(stat="identity", position=position_dodge()) + coord_flip() + theme_minimal()


##################
# for GWAS with cov
covneilmm2014 = c()
for(trait in traits){
  print(trait)

  pos = g$gmap
  geno = g$geno
  pheno = y2014
  smap = map
  g_nei = nei_coval(geno,smap,scale=1+0.01,grouping=pheno$BLK,n_core=1L)
  
  X0 = model.matrix(~factor(pheno$BLK))
  res = covnei_lmm(geno,g_nei,pheno=pheno[,trait],X0=X0,EM3.res12=EM3y2014res12[[trait]],K1=selfK,K2=neiK,K12=K12,K21=K21,n_core=4L)
  
  covneilmm2014 = c(covneilmm2014,list(res))
}
names(covneilmm2014) = traits
saveRDS(covneilmm2014,file="../output/aspen/aspen_covneilmm2014.rds")

covneilmm2014 = readRDS(file="../output/aspen/aspen_covneilmm2014.rds")
for(i in 1:21) {
  print(names(covneilmm2014)[i])
  print(sum(p.adjust(covneilmm2014[[i]]$p2,method="fdr")<0.05,na.rm=TRUE))
}


####
# cross-validation
trait = "BAI2014_2015"

pos = g$gmap
geno = g$geno
pheno = y2014
smap = map

g_nei = nei_coval(geno,smap,scale=1+0.01,grouping=pheno$BLK,n_core=1L)

# calc. K
selfCross = tcrossprod(geno)
neiCross = tcrossprod(g_nei)

q = ncol(geno)
selfK = (q+selfCross)/(2*(q-1))
neiK = neiCross/(q-1)

y = as.matrix(pheno[,trait])
ind_names = paste0("ind",1:length(y))
rownames(y) = ind_names
Z = diag(length(y))
rownames(Z) = ind_names
colnames(Z) = ind_names
rownames(selfK) = ind_names
colnames(selfK) = ind_names
rownames(neiK) = ind_names
colnames(neiK) = ind_names
ZETA1 <- list(A = list(Z = Z, K = selfK),
              Int = list(Z = Z, K = neiK))

K12 = tcrossprod(geno,g_nei) / (q-1)
K21 = t(K12)
rownames(K12) = ind_names
colnames(K12) = ind_names
rownames(K21) = ind_names
colnames(K21) = ind_names

covList1 <- rep(list(rep(list(NULL), 2)), 2)
covList1[[1]][[2]] <- K12
covList1[[2]][[1]] <- K21

X0 = model.matrix(~factor(pheno$BLK))
rownames(X0) = ind_names

nfold = 10
set.seed(1234)
foldid = sample(1:nfold,length(ind_names),replace=TRUE)
org_y = y

cv_i = function(i){
  y[foldid==i] = NA
  
  EM3.res11 <- EM3.cpp(y = y,
                       X0 = X0,
                       ZETA = ZETA1)
  w = EM3.res11$weights
  if(sum(w==TRUE)==1) {
    w[w==1] = 0.99
    w[w==0] = 0.01
  }
  EM3.res12 <- EM3.cov(y = y,
                       X0 = X0,
                       ZETA = ZETA1,
                       covList = covList1,
                       optimizeWeights = TRUE,
                       parInitForWeights = w,
                       forceApproxK = TRUE,
                       nIterOptimization=10^5, traceInside=1)
  
  R11 = cor.test(EM3.res11$y.pred[foldid==i],org_y[foldid==i])$estimate^2
  R12 = cor.test(EM3.res12$y.pred[foldid==i],org_y[foldid==i])$estimate^2
  
  mae11 = mean(abs(EM3.res11$y.pred[foldid==i]-org_y[foldid==i]),na.rm=TRUE)
  mae12 = mean(abs(EM3.res12$y.pred[foldid==i]-org_y[foldid==i]),na.rm=TRUE)
  
  return(c(R11,R12,mae11,mae12))
}

cv2014 = parallel::mcmapply(cv_i,1:10,mc.cores=10L)
cv2014 = t(cv2014)
colnames(cv2014) = c("R2_noncov","R2_cov","MAE_noncov","MAE_cov")
cv2014 = as.data.frame(cv2014)

r2p14 = cv2014 %>%
  gather(c("R2_noncov","R2_cov"),key="model",value="R2") %>%
  ggplot(aes(x=model,y=R2)) + geom_boxplot() + geom_jitter() + labs(title="2014") + ylim(0,0.5) + theme_bw()

maep14 = cv2014 %>%
  gather(c("MAE_noncov","MAE_cov"),key="model",value="MAE") %>%
  ggplot(aes(x=model,y=MAE)) + geom_boxplot() + geom_jitter() + theme_bw()

ggsave(r2p14 + maep14, filename=paste0("../../cv_y2014_",trait,".pdf"),width = 6,height = 3)

setwd("../")


