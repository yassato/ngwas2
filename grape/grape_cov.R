############
# grape INNOVINE data (Flutre et al. 2022 G3)
library(rNeighborGWAS)
library(gaston)
source("covnei_lmm.R")
library(RAINBOWR)
library(tidyverse)
library(patchwork)

# load genotypes
geno = read.vcf("./grape/Flutre2022data/snp-chip-gbs-genotypes_concat-fill_beagle_tag-filter-ld9_rm-dup.vcf.gz",convert.chr = FALSE)
geno[,which(geno@snps$maf>0.05)] # MAF cut-off 5%
geno@snps$chr[geno@snps$chr == "chrUkn"] = "chr100"
geno@snps$chr = as.numeric(gsub("chr","",geno@snps$chr))
geno = gaston2neiGWAS(geno)
rownames(geno$geno) = paste0("X",rownames(geno$geno))


# load phenotypes and annotations
pheno = read.delim("./grape/Flutre2022data/DLVitis_phenotypes.tsv")
pheno$citr.norm[which(pheno$citr.norm>80)] = NA
pheno$block = as.numeric(factor(pheno$block))
ann = read.delim("./grape/Flutre2022data/DLVitis_genotype-identifiers.tsv")

cltvr_no = c()
for(i in pheno$code.intro){
  ID = ann[which(ann$code.intro == i),]$cultivar.number
  if(length(ID)==1) {
    cltvr_no = c(cltvr_no, ID)
  } else {
    cltvr_no = c(cltvr_no, NA)
  }
}
cltvr_no = paste0("X",cltvr_no)

pheno = data.frame(cltvr_no,pheno)
pheno = pheno[is.na(pheno$cltvr_no)!=TRUE,]
pheno = pheno[(pheno$cltvr_no %in% row.names(geno$geno))==TRUE,]
smap = data.frame(pheno$rank,pheno$location)
colnames(smap) = c("x","y")

g_self = geno$geno[pheno$cltvr_no,]
group = paste0(pheno$harvest,"_",pheno$block,"_",pheno$rank) #location 1m < rank 2.5m
g_nei = nei_coval(geno=g_self,smap=smap,grouping=group,scale=1+0.01) #location 1m < rank 2.5m

# convert date to days
mm = as.numeric(substr(gsub("-","",pheno$datever),6,6))*31 - 7*31
dd = as.numeric(substr(gsub("-","",pheno$datever),7,8))
pheno$datever = mm + dd


########
# PVE including cov
for(trait in colnames(pheno)[8:26]) { print(trait)
  EM3res11 = c()
  EM3res12 = c()
  for(i in c(2010,2011,2012)){ print(i)
    #subset into years
    pheno2 = pheno[pheno$harvest==i,]
    g_self2 = g_self[pheno$harvest==i,]
    g_nei2 = g_nei[pheno$harvest==i,]
    
    y = as.matrix(pheno2[,trait])
    
    # calc. K
    selfCross = tcrossprod(g_self2)
    neiCross = tcrossprod(g_nei2)
    
    q = ncol(g_self2)
    selfK = (q+selfCross)/(2*(q-1))
    neiK = neiCross/(q-1)
    
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
    
    K12 = tcrossprod(g_self2,g_nei2) / (q-1)
    K21 = t(K12)
    rownames(K12) = ind_names
    colnames(K12) = ind_names
    rownames(K21) = ind_names
    colnames(K21) = ind_names
    
    covList1 <- rep(list(rep(list(NULL), 2)), 2)
    covList1[[1]][[2]] <- K12
    covList1[[2]][[1]] <- K21
    
    X0 = model.matrix(~factor(pheno2$block)+poly(pheno2$rank,2)+poly(pheno2$location,2))
    rownames(X0) = ind_names
    
    EM3.res11 <- EM3.cpp(y = y,
                         X0 = X0,
                         ZETA = ZETA1)
    w = EM3.res11$weights
    if(sum(w==TRUE)==1) {
      w[w==1] = 0.99
      w[w==0] = 0.01
    }
    EM3.res12 <- try(EM3.cov(y = y,
                             X0 = X0,
                             ZETA = ZETA1,
                             covList = covList1,
                             optimizeWeights = TRUE,
                             parInitForWeights = w,
                             forceApproxK = TRUE,
                             nIterOptimization = 10^5, traceInside = 1))
    
    EM3res11 = c(EM3res11,list(EM3.res11))
    EM3res12 = c(EM3res12,list(EM3.res12))
  }
  
  names(EM3res11) = c("harvest2010","harvest2011","harvest2012")
  names(EM3res12) = c("harvest2010","harvest2011","harvest2012")
  
  saveRDS(EM3res11,file=paste0("./output/grape/EM3res11_",trait,".rds"))
  saveRDS(EM3res12,file=paste0("./output/grape/EM3res12_",trait,".rds"))
}

#########
#SE calc
for(trait in colnames(pheno)[8:26]) { print(trait)
  EM3res11 = try(readRDS(paste0("./output/grape/EM3res11_",trait,".rds")))
  EM3res12 = try(readRDS(paste0("./output/grape/EM3res12_",trait,".rds")))
  
  EM3res11se = c()
  EM3res12se = c()
  for(i in c(2010,2011,2012)){ print(i)
    #subset into years
    pheno2 = pheno[pheno$harvest==i,]
    g_self2 = g_self[pheno$harvest==i,]
    g_nei2 = g_nei[pheno$harvest==i,]
    
    y = as.matrix(pheno2[,trait])
    
    # calc. K
    selfCross = tcrossprod(g_self2)
    neiCross = tcrossprod(g_nei2)
    
    q = ncol(g_self2)
    selfK = (q+selfCross)/(2*(q-1))
    neiK = neiCross/(q-1)
    
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
    
    K12 = tcrossprod(g_self2,g_nei2) / (q-1)
    K21 = t(K12)
    rownames(K12) = ind_names
    colnames(K12) = ind_names
    rownames(K21) = ind_names
    colnames(K21) = ind_names
    
    covList1 <- rep(list(rep(list(NULL), 2)), 2)
    covList1[[1]][[2]] <- K12
    covList1[[2]][[1]] <- K21
    
    X0 = model.matrix(~factor(pheno2$block)+poly(pheno2$rank,2)+poly(pheno2$location,2))
    rownames(X0) = ind_names
    
    EM3.res = EM3res11[[paste0("harvest",i)]]
    EM3.SE.res <- try(calcParamSE(
      y = y,
      MM.res = EM3.res,
      ZETA = ZETA1,
      X0 = X0
    ))
    
    EM3cov.res = EM3res12[[paste0("harvest",i)]]
    EM3cov.SE.res <- try(calcParamSE.cov(
      y = y,
      EM3.cov.res = EM3cov.res,
      ZETA = ZETA1,
      covList = covList1,
      X0 = X0,
      REML = TRUE
    ))
    
    EM3res11se = c(EM3res11se,list(EM3.SE.res))
    EM3res12se = c(EM3res12se,list(EM3cov.SE.res))
  }
  
  names(EM3res11se) = c("harvest2010","harvest2011","harvest2012")
  names(EM3res12se) = c("harvest2010","harvest2011","harvest2012")
  
  saveRDS(EM3res11se,file=paste0("./output/grape/EM3res11se_",trait,".rds"))
  saveRDS(EM3res12se,file=paste0("./output/grape/EM3res12se_",trait,".rds"))
}

########
# for IGE model
library(RAINBOWR)
for(trait in colnames(pheno)[8:26]) { print(trait)
  
  ZAdj = matrix(rep(0,nrow(pheno)^2),nrow(pheno),nrow(pheno))
  grouping = factor(paste0(pheno$harvest,"_",pheno$block,"_",pheno$rank))
  for(i in levels(grouping)){
    sub = pheno[grouping==i,]
    m = as.matrix(dist(sub$location))
    Zsub = mapply(FUN=function(x){as.numeric(m[,x]<=1.01&m[,x]>0)},1:nrow(m))
    ZAdj[grouping==i,grouping==i] = Zsub
  }
  
  EM3res11 = c()
  EM3res12 = c()
  for(i in c(2010,2011,2012)){ print(i)
    #subset into years
    pheno2 = pheno[pheno$harvest==i,]
    g_self2 = g_self[pheno$harvest==i,]
    ZAdj2 = ZAdj[pheno$harvest==i,pheno$harvest==i]
    
    y = as.matrix(pheno2[,trait])
    
    # calc. K
    selfCross = tcrossprod(g_self2)
    
    q = ncol(g_self2)
    selfK = (q+selfCross)/(2*(q-1))
    
    ind_names = paste0("ind",1:length(y))
    rownames(y) = ind_names
    Z = diag(length(y))
    rownames(Z) = ind_names
    colnames(Z) = ind_names
    rownames(selfK) = ind_names
    colnames(selfK) = ind_names
    
    # for IGE model
    rownames(ZAdj2) = ind_names
    colnames(ZAdj2) = ind_names
    ZETA2 <- list(A = list(Z = Z, K = selfK),
                  Adj = list(Z = ZAdj2, K = selfK))
    
    # for IGE model
    covList2 <- rep(list(rep(list(NULL), 2)), 2)
    covList2[[1]][[2]] <- selfK
    covList2[[2]][[1]] <- selfK
    
    X0 = model.matrix(~factor(pheno2$block)+poly(pheno2$rank,2)+poly(pheno2$location,2))
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
    
    EM3res21 = c(EM3res21,list(EM3.res21))
    EM3res22 = c(EM3res22,list(EM3.res22))
  }
  
  names(EM3res21) = c("harvest2010","harvest2011","harvest2012")
  names(EM3res22) = c("harvest2010","harvest2011","harvest2012")
  
  saveRDS(EM3res21,file=paste0("./output/grape/EM3res21_",trait,".rds"))
  saveRDS(EM3res22,file=paste0("./output/grape/EM3res22_",trait,".rds"))
}

#####################################
# compare LL between res12 and res22
ll12 = c(); ll21 = c(); ll22 = c()
rho22 = c(); w22 = c()
for(i in colnames(pheno)[8:26]){ print(i)
  EM3.res12 = readRDS(paste0("../output/grape/EM3res12_",i,".rds"))
  EM3.res21 = readRDS(paste0("../output/grape/EM3res21_",i,".rds"))
  EM3.res22 = readRDS(paste0("../output/grape/EM3res22_",i,".rds"))
  for(j in 1:3){
    if(is.na(EM3.res22[[j]])[1]!=TRUE){
      ll12 = c(ll12,EM3.res12[[j]]$LL)
      ll21 = c(ll21,EM3.res21[[j]]$LL)
      ll22 = c(ll22,EM3.res22[[j]]$LL)
      w22 = c(w22,EM3.res22[[j]]$weights[2])
      rho22 = c(rho22,EM3.res22[[j]]$rhosMat[1,2])
    } else {
      ll12 = c(ll12,NA)
      ll21 = c(ll21,NA)
      ll22 = c(ll22,NA)
      w22 = c(w22,NA)
      rho22 = c(rho22,NA)
    }
  }
}
which(pchisq(-2*(ll21-ll22),1,lower.tail=FALSE)<0.05) # no significant case

wilcox.test(ll12,ll22,paired=TRUE)
grapeLLp = ggplot(data=data.frame(diffLL=ll12-ll22),aes(x=diffLL)) + geom_histogram(fill=grey(0.5,0.5)) + geom_boxplot(colour="black",lwd=1.5) + theme_bw() + xlab("delta_LL")
ggsave(grapeLLp,filename="../figs/grapeLLp.pdf",height=6,width=6)

# export for Table S3
chi=-2*(ll21-ll22)
pval = pchisq(-2*(ll21-ll22),1,lower.tail=FALSE)
write.csv(data.frame(ll12,w22,rho22,ll22,chi,pval),"../tmp.csv")

# check SE
pve_res = c()
for(trait in colnames(pheno)[8:26]){
  EM3.res12 = try(readRDS(paste0("../output/grape/EM3res12_",trait,".rds")))
  EM3.res11 = try(readRDS(paste0("../output/grape/EM3res11_",trait,".rds")))
  EM3.res12se = readRDS(file=paste0("./output/grape/EM3res12se_",trait,".rds"))
  
  for(j in 1:3){
    if(is.na(EM3.res12[[j]])[1]!=TRUE){
      LL0 = EM3.res11[[j]]$LL
      chi = 2*(EM3.res12[[j]]$LL - LL0)
      pchi = pchisq(chi,1,lower.tail=FALSE)
      if(class(EM3.res12se[[j]])!="try-error"){
        pve_res = rbind(pve_res,c(EM3.res12se[[j]]$parameterTable[1,3],EM3.res12se[[j]]$parameterTable[1,4],
                                  EM3.res12se[[j]]$parameterTable[3,3],EM3.res12se[[j]]$parameterTable[3,4],
                                  EM3.res12se[[j]]$parameterTable[2,3],EM3.res12se[[j]]$parameterTable[2,4],
                                  EM3.res12se[[j]]$parameterTable[4,3],EM3.res12se[[j]]$parameterTable[4,4],
                                  chi,pchi))
      } else {
        pve_res = rbind(pve_res,c(EM3.res12[[j]]$Vu*EM3.res12[[j]]$weights[1],NA,
                                  EM3.res12[[j]]$Vu*EM3.res12[[j]]$weights[2],NA,
                                  EM3.res12[[j]]$rhosMat[1,2],NA,
                                  EM3.res12[[j]]$Ve,NA,
                                  chi,pchi))
      }
      
    } else {
      pve_res = rbind(pve_res,c(rep(NA,10)))
    }
  }
}

# show results for Table S3
colnames(pve_res) = c("sigma1","sigma1_se","sigma2","sigma2_se","sigma12","sigma12_se","sigmae","sigmae_se","chisq","pval")
rownames(pve_res) = paste0(rep(colnames(pheno)[8:26],each=3),"_",rep(names(EM3.res12),length(colnames(pheno)[8:26])))
pve_res = as.data.frame(pve_res)
write.csv(pve_res,"../grape_tableS3.csv")

# plot Figure 5
pve_res = c()
traits = colnames(pheno)[8:26]
for(trait in traits){
  EM3.res12 = try(readRDS(paste0("./output/grape/EM3res12_",trait,".rds")))
  
  for(j in 1:3){
    if(is.na(EM3.res12[[j]])[1]!=TRUE){
      tmp = c(EM3.res12[[j]]$Vu,EM3.res12[[j]]$Ve,EM3.res12[[j]]$weights,EM3.res12[[j]]$rhosMat[1,2])
      pve_res = rbind(pve_res,tmp)
    } else {
      pve_res = rbind(pve_res,c(rep(NA,5)))
    }
  }
}
pve_res = as.data.frame(pve_res)
colnames(pve_res) = c("Vu","Ve","w1","w2","rho12")
rownames(pve_res) = c(1:nrow(pve_res))
trait = rep(traits,each=3)
year = rep(names(EM3.res12),length(traits))
h2 = pve_res$Vu / (pve_res$Vu+pve_res$Ve)
pve_res = data.frame(trait,year,pve_res,DGE=h2*pve_res$w1,IGE=h2*pve_res$w2,cov=h2*pve_res$rho12*sqrt(pve_res$w1*pve_res$w2))
pve_res$trait=factor(pve_res$trait, levels=colnames(pheno)[8:26])

pve_grape2010 = pve_res %>%
  filter(year=="harvest2010") %>%
  gather(c("DGE","IGE","cov"),key="factor",value="PVE") %>%
  ggplot(aes(x=trait, y=PVE, fill=factor)) + ylim(-0.25,1) + labs(subtitle = "2010") +
  geom_bar(stat="identity", position=position_dodge()) + coord_flip() + theme_minimal()

pve_grape2011 = pve_res %>%
  filter(year=="harvest2011") %>%
  gather(c("DGE","IGE","cov"),key="factor",value="PVE") %>%
  ggplot(aes(x=trait, y=PVE, fill=factor)) + ylim(-0.25,1) + labs(subtitle = "2011") +
  geom_bar(stat="identity", position=position_dodge()) + coord_flip() + theme_minimal()

pve_grape2012 = pve_res %>%
  filter(year=="harvest2012") %>%
  gather(c("DGE","IGE","cov"),key="factor",value="PVE") %>%
  ggplot(aes(x=trait, y=PVE, fill=factor)) + ylim(-0.25,1) + labs(subtitle = "2012") +
  geom_bar(stat="identity", position=position_dodge()) + coord_flip() + theme_minimal()

pve_grape = pve_grape2010+pve_grape2011+pve_grape2012+plot_annotation(tag_levels = "a")
ggsave(pve_grape,filename="../grape_pve.pdf",width=12,height=4)

#############
# for GWAS
for(trait in colnames(pheno)[8:26]){ print(trait)
  covneilmm_res = c()
  EM3res12 = try(readRDS(paste0("./output/grape/EM3res12_",trait,".rds")))
  for(i in c(2010,2011,2012)){ print(i)
    
    obj = paste0("harvest",i)
    if(is.na(EM3res12[[obj]])[1]) {
      covneilmm_res = c(covneilmm_res,list(NA))
      next
    } else if(class(EM3res12[[obj]])=="try-error") {
      covneilmm_res = c(covneilmm_res,list(NA))
      next
    }
    
    #subset into years
    pheno2 = pheno[pheno$harvest==i,]
    g_self2 = g_self[pheno$harvest==i,]
    g_nei2 = g_nei[pheno$harvest==i,]
    
    y = as.matrix(pheno2[,trait])
    
    # calc. K
    selfCross = tcrossprod(g_self2)
    neiCross = tcrossprod(g_nei2)
    
    q = ncol(g_self2)
    selfK = (q+selfCross)/(2*(q-1))
    neiK = neiCross/(q-1)
    
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
    
    K12 = tcrossprod(g_self2,g_nei2) / (q-1)
    K21 = t(K12)
    rownames(K12) = ind_names
    colnames(K12) = ind_names
    rownames(K21) = ind_names
    colnames(K21) = ind_names
    
    covList1 <- rep(list(rep(list(NULL), 2)), 2)
    covList1[[1]][[2]] <- K12
    covList1[[2]][[1]] <- K21
    
    X0 = model.matrix(~factor(pheno2$block)+poly(pheno2$rank,2)+poly(pheno2$location,2))
    rownames(X0) = ind_names
    res = covnei_lmm(geno=g_self2,g_nei=g_nei2,pheno=y[,1],X0=X0,EM3.res12=EM3res12[[obj]],K1=selfK,K2=neiK,K12=K12,K21=K21,n_core=6L)
    covneilmm_res = c(covneilmm_res,list(res))
  }
  names(covneilmm_res) = c("harvest2010","harvest2011","harvest2012")
  saveRDS(covneilmm_res,file=paste0("./output/grape/covnei_lmm_",trait,".rds"))
}
  

# create Manhattan and QQ plots
pdf()
nsnps = c()
gifs = c()
n_ind = c()
for(trait in colnames(pheno)[8:26]){
  n_ind = c(n_ind,
  sum(is.na(pheno[pheno$harvest=="2010",trait])==FALSE),
  sum(is.na(pheno[pheno$harvest=="2011",trait])==FALSE),
  sum(is.na(pheno[pheno$harvest=="2012",trait])==FALSE)
  )
  
  res = try(readRDS(paste0("./output/grape/covnei_lmm_",trait,".rds")))
  for(i in 1:3){
    if(is.na(res[[i]])[1]!=TRUE) {
      gaston::manhattan(x=data.frame(chr=geno$gmap$chr,pos=geno$gmap$pos,p=res[[i]]$p2),las=1,main=paste(trait,i))
      abline(h=-log10(0.05/length(res[[i]]$p2)))
      gaston::qqplot.pvalues(p=res[[i]]$p2,main=trait,las=1)
      nsnps = c(nsnps,sum(p.adjust(res[[i]]$p2,method="BH")<0.05))
      gifs = c(gifs,median(log10(res[[i]]$p2))/median(log10(ppoints(length(res[[i]]$p2)))))
    } else {
      nsnps = c(nsnps,NA)
      gifs = c(gifs,NA)
    }
  }
}
dev.off()

# export for Table S3
write.csv(cbind(nsnps,gifs,n_ind),"../tmp.csv")
