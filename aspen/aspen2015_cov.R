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
y2015 = read.csv("WisAsp_Data2015.csv",header=TRUE)
y2015$ID = gsub(" ", "_",y2015$ID)
genetID = y2015[(y2015$ID %in% ann$Correct.Sequenced.Ramet),"GENO"]

gIDs = rep(NA,nrow(y2015))
for(i in 1:nrow(y2015)){
  focal = which(ann$Correct.Sequenced.Ramet==y2015$ID[i])
  print(focal)
  if(length(focal)==1){
    gIDs[which(y2015$GENO==y2015$GENO[i])] = ann$ID.RapidGenomics[focal]
  }
}
y2015 = data.frame(gIDs,y2015)
y2015 = y2015[(y2015$gIDs %in% gsub(".bam","",g@ped$id)),]

gwhich = c()
for(i in y2015$gIDs) {
  focal = which(gsub(".bam","",g@ped$id)==i)
  gwhich = c(gwhich,focal)
}

g = g[gwhich,]

gsub(".bam","",g@ped$id) == y2015$gIDs
y2015$ROW = as.numeric(factor(y2015$ROW))

map = data.frame(x=y2015$ROW,y=y2015$POS)

g = rNeighborGWAS::gaston2neiGWAS(g)
AF = apply(g$geno,2,mean)
g$geno = g$geno[,which(AF<0.95 & AF>0.05)]
g$gmap = g$gmap[which(AF<0.95 & AF>0.05),]

######
# PVE including covariance
y2015[,8:44] = log(y2015[,8:44]+0.5)
traits = colnames(y2015[,c(8:10,45:ncol(y2015))])
EM3y2015res11 = c()
EM3y2015res12 = c()
EM3y2015res12se = c()
for(trait in traits){ print(trait)

  pos = g$gmap
  geno = g$geno
  pheno = y2015
  smap = map
  
  # refer to the nearest neighbors only
  g_nei = nei_coval(geno,smap,scale=1+0.01,grouping=pheno$BLK,n_core=1L)
  
  # calc. K
  selfCross = tcrossprod(geno)
  neiCross = tcrossprod(g_nei)
  
  q = ncol(geno)
  selfK = (q+selfCross)/(2*(q-1))
  neiK = neiCross/(q-1)
  
  # 20250619 added
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
  
  X0 = model.matrix(~factor(pheno$BLK)+poly(pheno$ROW,2)+poly(pheno$POS,2))
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
  
  EM3cov.SE.res <- try(calcParamSE.cov(
    y = y,
    EM3.cov.res = EM3.res12,
    ZETA = ZETA1,
    covList = covList1,
    X0 = X0,
    REML = TRUE
  ))
  
  EM3y2015res11 = c(EM3y2015res11,list(EM3.res11))
  EM3y2015res12 = c(EM3y2015res12,list(EM3.res12))
  EM3y2015res12se = c(EM3y2015res12se,list(EM3cov.SE.res))
  
}
names(EM3y2015res11) = traits
names(EM3y2015res12) = traits
names(EM3y2015res12se) = traits

saveRDS(EM3y2015res11,file="../output/aspen/aspen_2015_EM3res11_rc.rds")
saveRDS(EM3y2015res12,file="../output/aspen/aspen_2015_EM3res12_rc.rds")
saveRDS(EM3y2015res12se,file="../output/aspen/aspen_2015_EM3res12se_rc.rds")


# reload
EM3y2015res11 = readRDS(file="../output/aspen/aspen_2015_EM3res11_rc.rds")
EM3y2015res12 = readRDS(file="../output/aspen/aspen_2015_EM3res12_rc.rds")
EM3y2015res12se = readRDS(file="../output/aspen/aspen_2015_EM3res12se_rc.rds")


####
# for conventional IGE model
traits = colnames(y2015[,c(8:10,45:ncol(y2015))])
EM3y2015res21 = c()
EM3y2015res22 = c()
for(trait in traits){ print(trait)
  pos = g$gmap
  geno = g$geno
  pheno = y2015
  smap = map
  
  # calc. K
  selfCross = tcrossprod(geno)
  q = ncol(geno)
  selfK = (q+selfCross)/(2*(q-1))
  
  ZAdj = matrix(rep(0,nrow(pheno)^2),nrow(pheno),nrow(pheno))
  grouping = factor(pheno$BLK)
  for(i in levels(grouping)){
    sub = pheno[grouping==i,]
    m = as.matrix(dist(cbind(sub$ROW,sub$POS)))
    Zsub = mapply(FUN=function(x){as.numeric(m[,x]<=1.01&m[,x]>0)},1:nrow(m))
    ZAdj[grouping==i,grouping==i] = Zsub
  }
  
  y = as.matrix(pheno[,trait])
  ind_names = paste0("ind",1:length(y))
  rownames(y) = ind_names
  Z = diag(length(y))
  rownames(Z) = ind_names
  colnames(Z) = ind_names
  rownames(selfK) = ind_names
  colnames(selfK) = ind_names
  
  # for IGE model
  rownames(ZAdj) = ind_names
  colnames(ZAdj) = ind_names
  ZETA2 <- list(A = list(Z = Z, K = selfK),
                Adj = list(Z = ZAdj, K = selfK))
  
  # for IGE model
  covList2 <- rep(list(rep(list(NULL), 2)), 2)
  covList2[[1]][[2]] <- selfK
  covList2[[2]][[1]] <- selfK
  
  X0 = model.matrix(~factor(pheno$BLK)+poly(smap[,1],2)+poly(smap[,2],2))
  rownames(X0) = ind_names
  
  # for IGE model
  EM3.res21 <- EM3.cpp(y = y,
                       X0 = X0,
                       ZETA = ZETA2)
  w = EM3.res21$weights
  if(sum(w==TRUE)==1) {
    w[w==1] = 0.99
    w[w==0] = 0.01
  }
  EM3.res22 <- try(EM3.cov(y = y,
                           X0 = X0,
                           ZETA = ZETA2,
                           covList = covList2,
                           optimizeWeights = TRUE,
                           parInitForWeights = w,
                           forceApproxK = TRUE,
                           nIterOptimization = 10^5, traceInside = 1))
  
  EM3y2015res21 = c(EM3y2015res21,list(EM3.res21))
  EM3y2015res22 = c(EM3y2015res22,list(EM3.res22))
}
names(EM3y2015res21) = traits
names(EM3y2015res22) = traits

saveRDS(EM3y2015res21,file="../../output/aspen_2015_EM3res21_rc1.rds")
saveRDS(EM3y2015res22,file="../../output/aspen_2015_EM3res22_rc1.rds")


####
# LL comparison
EM3y2015res21 = readRDS(file="../output/aspen/aspen_2015_EM3res21_rc1.rds")
EM3y2015res22 = readRDS(file="../output/aspen/aspen_2015_EM3res22_rc1.rds")
ll12 = c(); ll21 = c(); ll22 = c()
for(i in names(EM3y2015res22)){
  ll12 = c(ll12,EM3y2015res12[[i]]$LL)
  ll21 = c(ll21,EM3y2015res21[[i]]$LL)
  ll22 = c(ll22,EM3y2015res22[[i]]$LL)
  print(c(EM3y2015res22[[i]]$weights[2],EM3y2015res22[[i]]$rhosMat[1,2],EM3y2015res22[[i]]$LL))
}
boxplot(ll12-ll22)
wilcox.test(ll12,ll22,paired=TRUE)

ll1 = c(ll12_2014,ll12)
ll2 = c(ll22_2014,ll22)
boxplot(ll1-ll2)
wilcox.test(ll1,ll2,paired=TRUE)
aspenLLp = ggplot(data=data.frame(diffLL=ll1-ll2),aes(x=diffLL)) + geom_histogram(fill=grey(0.5,0.5)) + geom_boxplot(colour="black",lwd=1.5) + theme_bw() + xlab("delta_LL")
ggsave(aspenLLp,filename="../../figs/aspenLLp.pdf",height=6,width=6)

# show results for Table S1
se2015 = c()
for(i in traits){
  chi = -2*(EM3y2015res11[[i]]$LL-EM3y2015res12[[i]]$LL)
  pchi = pchisq(chi,df=1,lower.tail=FALSE)
  if(class(EM3y2015res12se[[i]])!="try-error"){
    se2015 = rbind(se2015, c(EM3y2015res12se[[i]]$parameterTable$Estimate[1],EM3y2015res12se[[i]]$parameterTable$SE[1],
                             EM3y2015res12se[[i]]$parameterTable$Estimate[3],EM3y2015res12se[[i]]$parameterTable$SE[3],
                             EM3y2015res12se[[i]]$parameterTable$Estimate[2],EM3y2015res12se[[i]]$parameterTable$SE[2],
                             EM3y2015res12se[[i]]$parameterTable$Estimate[4],EM3y2015res12se[[i]]$parameterTable$SE[4],
                             chi,pchi)
    )
  } else {
    se2015 = rbind(se2015, c(EM3y2015res12[[i]]$Vu*EM3y2015res12[[i]]$weights[1],NA,
                             EM3y2015res12[[i]]$Vu*EM3y2015res12[[i]]$weights[2],NA,
                             EM3y2015res12[[i]]$rhosMat[1,2],NA,
                             EM3y2015res12[[i]]$Ve,NA,
                             chi,pchi)
    )
  }
}
colnames(se2015) = c("sigma1","sigma1_se","sigma2","sigma2_se","sigma12","sigma12_se","sigmae","sigmae_se","Chi","p")
rownames(se2015) = traits
N = apply(is.na(y2015[,c(8:10,45:ncol(y2015))])==FALSE,2,sum)
se2015 = data.frame(N,se2015)
write.csv(se2015,"../../aspen2015_TableS1b.csv")


# plot PVEs
sum2015 = c()
for(i in 1:length(EM3y2015res12)) sum2015 = rbind(sum2015,c(EM3y2015res12[[i]]$Vu,EM3y2015res12[[i]]$Ve,EM3y2015res12[[i]]$weights,EM3y2015res12[[i]]$rhosMat[1,2],EM3y2015res12[[i]]$LL))
sum2015 = as.data.frame(sum2015)
rownames(sum2015) = names(EM3y2015res12)
colnames(sum2015) = c("Vu","Ve","w1","w2","rho12","LL")

h2 = sum2015$Vu / (sum2015$Vu + sum2015$Ve)
rhow = h2*sum2015$rho12*sqrt(sum2015$w1*sum2015$w2)

pve2015 = data.frame(trait=factor(rownames(sum2015),levels=rownames(sum2015)),DGE=h2*sum2015$w1,IGE=h2*sum2015$w2,cov=h2*sum2015$rho12*sqrt(sum2015$w1*sum2015$w2))
p2015 = pve2015 %>%
  gather(c("DGE","IGE","cov"),key="factor",value="PVE") %>%
  ggplot(aes(x=trait, y=PVE, fill=factor)) + ylim(-0.25,1) + labs(subtitle="2015") +
  geom_bar(stat="identity", position=position_dodge()) + coord_flip() + theme_minimal()

# export Figure 3
ggsave(p2014+p2015,filename="../../figs/aspen_pve.pdf",width=8,height=4)

##################
# for GWAS with cov
covneilmm2015 = c()
for(trait in traits){
  print(trait)

  pos = g$gmap
  geno = g$geno
  pheno = y2015
  smap = map
  g_nei = nei_coval(geno,smap,scale=1+0.01,grouping=pheno$BLK,n_core=1L)
  
  X0 = model.matrix(~factor(pheno$BLK)+poly(pheno$ROW,2)+poly(pheno$POS,2))
  res = covnei_lmm(geno,g_nei,pheno=pheno[,trait],X0=X0,EM3.res12=EM3y2015res12[[trait]],K1=selfK,K2=neiK,K12=K12,K21=K21,n_core=8L)
  
  covneilmm2015 = c(covneilmm2015,list(res))
}

names(covneilmm2015) = traits
saveRDS(covneilmm2015,file="../output/aspen/aspen_covneilmm2015rc.rds")

# no. of significant SNPs
covneilmm2015 = readRDS(file="../output/aspen/aspen_covneilmm2015rc.rds")
nsnps = c()
gifs = c()
for(i in 1:length(covneilmm2015)) {
  print(names(covneilmm2015)[i])
  gifs = c(gifs,median(log10(covneilmm2015[[i]]$p2),na.rm=TRUE)/median(log10(ppoints(length(covneilmm2015[[i]]$p2)))))
  nsnps = c(nsnps,sum(p.adjust(covneilmm2015[[i]]$p2,method="fdr")<0.05,na.rm=TRUE))
}
write.csv(cbind(nsnps,gifs),"../../tmp.csv") # export for Table S1


####
# cross-validation
trait = "BAI2014_2015"

pos = g$gmap
geno = g$geno
pheno = y2015
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

set.seed(1234)
cv2015 = parallel::mcmapply(cv_i,1:10,mc.cores=10L)
cv2015 = t(cv2015)
colnames(cv2015) = c("R2_noncov","R2_cov","MAE_noncov","MAE_cov")
cv2015 = as.data.frame(cv2015)

r2p15 = cv2015 %>%
  gather(c("R2_noncov","R2_cov"),key="model",value="R2") %>%
  ggplot(aes(x=model,y=R2)) + geom_boxplot() + geom_jitter() + labs(title="2015") + ylim(0,0.5) + theme_bw()

maep15 = cv2015 %>%
  gather(c("MAE_noncov","MAE_cov"),key="model",value="MAE") %>%
  ggplot(aes(x=model,y=MAE)) + geom_boxplot() + geom_jitter() + theme_bw()

ggsave(r2p15 + maep15, filename=paste0("../../cv_y2015_",trait,".pdf"),width = 6,height = 3)

# expoert a supplementary figure
cvp = (r2p14 + r2p15) / (maep14 + maep15) + plot_annotation(tag_levels="a")
ggsave(cvp, filename=paste0("../../cv_",trait,"rc.pdf"),width=5,height=5)

setwd("../")
