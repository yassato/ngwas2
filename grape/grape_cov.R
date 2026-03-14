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
group = paste0(pheno$block,"_",pheno$rank) #location 1m < rank 2.5m
g_nei = nei_coval(geno=g_self,smap=smap,grouping=group,scale=1+0.01) #location 1m < rank 2.5m

# convert date to days
mm = as.numeric(substr(gsub("-","",pheno$datever),6,6))*31 - 7*31
dd = as.numeric(substr(gsub("-","",pheno$datever),7,8))
pheno$datever = mm + dd

########
# PVE including cov
for(trait in colnames(pheno)[8:26]) {
  y = as.matrix(pheno[,trait])
  
  # calc. K
  selfCross = tcrossprod(g_self)
  neiCross = tcrossprod(g_nei)
  
  q = ncol(g_self)
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
  
  K12 = tcrossprod(g_self,g_nei) / (q-1)
  K21 = t(K12)
  rownames(K12) = ind_names
  colnames(K12) = ind_names
  rownames(K21) = ind_names
  colnames(K21) = ind_names
  
  covList1 <- rep(list(rep(list(NULL), 2)), 2)
  covList1[[1]][[2]] <- K12
  covList1[[2]][[1]] <- K21
  
  X0 = model.matrix(~factor(pheno$block))
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
                       nIterOptimization = 10^5, traceInside = 1)
  
  saveRDS(EM3.res11,file=paste0("./output/grape/EM3res11_",trait,".rds"))
  saveRDS(EM3.res12,file=paste0("./output/grape/EM3res12_",trait,".rds"))
}

# check PVE
pve_res = c()
for(trait in colnames(pheno)[8:26]){
  EM3.res12 = try(readRDS(paste0("./output/grape/EM3res12_",trait,".rds")))
  
  if(class(EM3.res12)!="try-error") {
    EM3.res11 = try(readRDS(paste0("./output/grape/EM3res11_",trait,".rds")))
    LL0 = EM3.res11$LL
    chi = 2*(EM3.res12$LL - LL0)
    pchi = pchisq(chi,1,lower.tail=FALSE)
    pve_res = rbind(pve_res,c(EM3.res12$Vu,EM3.res12$Ve,EM3.res12$weights,EM3.res12$rhosMat[1,2],chi,pchi))
  }
}

# show results for Table S3
colnames(pve_res) = c("Vu","Ve","w1","w2","rho12","chisq","pval")
rownames(pve_res) = colnames(pheno)[8:26]
pve_res = as.data.frame(pve_res)

# plot Figure 5
h2 = pve_res$Vu / (pve_res$Vu+pve_res$Ve)
pve_res = data.frame(trait=factor(rownames(pve_res),levels=rownames(pve_res)),pve_res,DGE=h2*pve_res$w1,IGE=h2*pve_res$w2,cov=h2*pve_res$rho12*sqrt(pve_res$w1*pve_res$w2))

pve_grape = pve_res %>%
  gather(c("DGE","IGE","cov"),key="factor",value="PVE") %>%
  ggplot(aes(x=trait, y=PVE, fill=factor)) + ylim(-0.25,1) +
  geom_bar(stat="identity", position=position_dodge()) + coord_flip() + theme_minimal()

ggsave(pve_grape,filename="../grape_pve.pdf",width=4,height=5)

#############
# for GWAS
for(trait in colnames(pheno)[8:26][1:2]){
  EM3.res12 = try(readRDS(paste0("./output/grape/EM3res12_",trait,".rds")))
  
  if(class(EM3.res12)!="try-error") {
    X0 = model.matrix(~factor(pheno$block))
    res = covnei_lmm(geno=g_self,g_nei=g_nei,pheno=pheno[,trait],X0=X0,EM3.res12=EM3.res12,K1=selfK,K2=neiK,K12=K12,K21=K21,n_core=8L)
    res = data.frame(geno$gmap,res)
    saveRDS(res,file=paste0("./output/grape/covnei_lmm_",trait,".rds"))
  }
}

# create Manhattan and QQ plots
pdf()
nsnps = c()
for(trait in colnames(pheno)[8:26]){
  res = try(readRDS(paste0("./output/grape/covnei_lmm_",trait,".rds")))
  if(res$p2[1]!="NaN") {
    gaston::manhattan(x=data.frame(chr=res$chr,pos=res$pos,p=res$p2),las=1,main=trait)
    abline(h=-log10(0.05/length(res$p2)))
    gaston::qqplot.pvalues(p=res$p2,main=trait,las=1)
    nsnps = c(nsnps,sum(p.adjust(res$p2,method="BH")<0.05))
  } else {
    nsnps = c(nsnps,NaN)
  }
}
dev.off()

# for Figure S6
pdf(width=12,height=9)
par(mfrow=c(2,1))
gaston::manhattan(x=data.frame(chr=res$chr,pos=res$pos,p=res$p2),las=1,main="(a)")
th = CalcThreshold(input=data.frame(name=c(1:length(res$p2)),chr=res$chr,pos=res$pos,p=-log10(res$p2)))
abline(h=th,lty=2)
gaston::qqplot.pvalues(p=res$p2,main="(b)",las=1)
dev.off()
