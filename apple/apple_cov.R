#####################
# apple REFPOP data

library(rNeighborGWAS)
library(gaston)
library(RAINBOWR)
library(lme4)
library(tidyverse)
library(patchwork)
source("coord.R")
source("covnei_lmm.R")

#############
# imputed genotypes provided by Michaela Jung
pheno = read.csv("./apple/REFPOP_Phenotypic_data_3years_non_adjusted_reduced_traits_paper_2_reduced_M1.txt", sep="\t")
pheno$Trunk_increment = log(pheno$Trunk_increment+0.01) #log transformation for growth
g3 = read.bed.matrix("./apple/SNPs_final_2022/SNPs_final_2022")
g3 = g3[,g3@snps$maf>0.05] # MAF cut-off at 5%
length(intersect(g3@ped$id,pheno$Genotype))
dim(g3@ped)

geno = gaston2neiGWAS(g3)
pheno$Genotype %in% rownames(geno$geno)
sum(length(intersect(pheno$Genotype,g3@ped$id))==FALSE)

#######
# PVE including cov
for(k in colnames(pheno)[7:30]){ #w/o rows and positions for Color_ground due to singularity
  EM3res11 = c()
  EM3res12 = c()
  for(i in c(2018,2019,2020)){
    for(j in levels(factor(pheno$Country))) {
      print(i)
      print(j)
      pheno2 = filter(pheno,Country==j&Year==i)
      geno2 = geno$geno[pheno2$Genotype,]
      smap2 = cbind(pheno2$Row,pheno2$Position)
      gnei2 = nei_coval(geno=geno2,smap=smap2,scale=1,grouping=pheno2$Row) #only neighbors within "Position"
      
      selfCross = tcrossprod(geno2)
      neiCross = tcrossprod(gnei2)
      
      q = ncol(geno2)
      selfK = (q+selfCross)/(2*(q-1))
      neiK = neiCross/(q-1)
      
      if(sum(is.na(pheno2[,k]))==length(pheno2[,k])) {
        EM3res11 = c(EM3res11,list(NA))
        EM3res12 = c(EM3res12,list(NA)) 
        next
      }
      
      y = as.matrix(pheno2[,k])
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
      
      K12 = tcrossprod(geno2,gnei2) / (q-1)
      K21 = t(K12)
      rownames(K12) = ind_names
      colnames(K12) = ind_names
      rownames(K21) = ind_names
      colnames(K21) = ind_names
      
      covList1 <- rep(list(rep(list(NULL), 2)), 2)
      covList1[[1]][[2]] <- K12
      covList1[[2]][[1]] <- K21
      
      # X0 = matrix(1,nrow=nrow(pheno2))
      X0 = model.matrix(~poly(pheno2$Row,2)+poly(pheno2$Position,2))
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
                           nIterOptimization=10^5, traceInside=1))
      
      EM3res11 = c(EM3res11,list(EM3.res11))
      EM3res12 = c(EM3res12,list(EM3.res12))  
    }
  }
  
  names(EM3res11) = paste0(rep(levels(factor(pheno$Country)),3),rep(c(2018,2019,2020),each=6))
  names(EM3res12) = paste0(rep(levels(factor(pheno$Country)),3),rep(c(2018,2019,2020),each=6))
  
  saveRDS(EM3res11,file=paste0("./output/apple/",k,"_EM3res11.rds"))
  saveRDS(EM3res12,file=paste0("./output/apple/",k,"_EM3res12.rds"))
}

#########
#SE calc
for(k in colnames(pheno)[7:30]){ #w/o rows and positions for Color_ground due to singularity
  EM3res11 = readRDS(file=paste0("./output/apple/",k,"_EM3res11.rds"))
  EM3res12 = readRDS(file=paste0("./output/apple/",k,"_EM3res12.rds"))
  EM3res11se = c()
  EM3res12se = c()
  for(i in c(2018,2019,2020)){
    for(j in levels(factor(pheno$Country))) {
      print(i)
      print(j)
      pheno2 = filter(pheno,Country==j&Year==i)
      geno2 = geno$geno[pheno2$Genotype,]
      smap2 = cbind(pheno2$Row,pheno2$Position)
      gnei2 = nei_coval(geno=geno2,smap=smap2,scale=1,grouping=pheno2$Row) #only neighbors within "Position"
      
      selfCross = tcrossprod(geno2)
      neiCross = tcrossprod(gnei2)
      
      q = ncol(geno2)
      selfK = (q+selfCross)/(2*(q-1))
      neiK = neiCross/(q-1)
      
      y = as.matrix(pheno2[,k])
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
      
      K12 = tcrossprod(geno2,gnei2) / (q-1)
      K21 = t(K12)
      rownames(K12) = ind_names
      colnames(K12) = ind_names
      rownames(K21) = ind_names
      colnames(K21) = ind_names
      
      covList1 <- rep(list(rep(list(NULL), 2)), 2)
      covList1[[1]][[2]] <- K12
      covList1[[2]][[1]] <- K21
      
      # X0 = matrix(1,nrow=nrow(pheno2))
      X0 = model.matrix(~poly(pheno2$Row,2)+poly(pheno2$Position,2))
      rownames(X0) = ind_names
      
      EM3.res = EM3res11[[paste0(j,i)]]
      EM3.SE.res <- try(calcParamSE(
        y = y,
        MM.res = EM3.res,
        ZETA = ZETA1,
        X0 = X0
        ))
      
      EM3.cov.res = EM3res12[[paste0(j,i)]]
      EM3cov.SE.res = try(calcParamSE.cov(y = y,
                               X0 = X0,
                               EM3.cov.res = EM3.cov.res,
                               ZETA = ZETA1,
                               covList = covList1
                               ))
      
      EM3res11se = c(EM3res11se,list(EM3.SE.res))
      EM3res12se = c(EM3res12se,list(EM3cov.SE.res))  
    }
  }
  names(EM3res11se) = paste0(rep(levels(factor(pheno$Country)),3),rep(c(2018,2019,2020),each=6))
  names(EM3res12se) = paste0(rep(levels(factor(pheno$Country)),3),rep(c(2018,2019,2020),each=6))
  
  saveRDS(EM3res11se,file=paste0("./output/apple/",k,"_EM3res11se.rds"))
  saveRDS(EM3res12se,file=paste0("./output/apple/",k,"_EM3res12se.rds"))
}



#####
# for IGE model
for(k in colnames(pheno)[7:30]){ #w/o rows and positions for Color_ground due to singularity
  EM3res21 = c()
  EM3res22 = c()
  for(i in c(2018,2019,2020)){
    for(j in levels(factor(pheno$Country))) {
      print(i)
      print(j)
      pheno2 = filter(pheno,Country==j&Year==i)
      geno2 = geno$geno[pheno2$Genotype,]
      smap2 = cbind(pheno2$Row,pheno2$Position)

      selfCross = tcrossprod(geno2)
      q = ncol(geno2)
      selfK = (q+selfCross)/(2*(q-1))

      if(sum(is.na(pheno2[,k]))==length(pheno2[,k])) {
        EM3res21 = c(EM3res21,list(NA))
        EM3res22 = c(EM3res22,list(NA)) 
        next
      }
      
      ZAdj = matrix(rep(0,nrow(pheno2)^2),nrow(pheno2),nrow(pheno2))
      for(l in levels(factor(pheno2$Row))){
        sub = pheno2[pheno2$Row==l,]
        m = as.matrix(dist(sub$Position))
        Zsub = mapply(FUN=function(x){as.numeric(m[,x]<=1.0&m[,x]>0)},1:nrow(m))
        ZAdj[pheno2$Row==l,pheno2$Row==l] = Zsub
      }
      
      y = as.matrix(pheno2[,k])
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
      
      if(k=="Color_ground") {
        X0 = matrix(1,nrow=nrow(pheno2))
      } else {
        X0 = model.matrix(~poly(pheno2$Row,2)+poly(pheno2$Position,2))
      }
      
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
  }
  
  names(EM3res21) = paste0(rep(levels(factor(pheno$Country)),3),rep(c(2018,2019,2020),each=6))
  names(EM3res22) = paste0(rep(levels(factor(pheno$Country)),3),rep(c(2018,2019,2020),each=6))
  
  saveRDS(EM3res21,file=paste0("./output/apple/",k,"_EM3res21.rds"))
  saveRDS(EM3res22,file=paste0("./output/apple/",k,"_EM3res22.rds"))
}



#####################
# GWAS for all traits
for(k in colnames(pheno)[7:30]){
  covneilmm_res = c()
  for(i in c(2018,2019,2020)){
    for(j in levels(factor(pheno$Country))) {
      
      EM3res12 = readRDS(file=paste0("./output/apple/apple/",k,"_EM3res12.rds"))
      print(k)
      print(i)
      print(j)
      
      if(is.na(EM3res12[[paste0(j,i)]])[1]) {
        covneilmm_res = c(covneilmm_res,list(NA))
        next
      } else if(class(EM3res12[[paste0(j,i)]])=="try-error") {
        covneilmm_res = c(covneilmm_res,list(NA))
        next
      }

      
      pheno2 = filter(pheno,Country==j&Year==i)
      geno2 = geno$geno[pheno2$Genotype,]
      smap2 = cbind(pheno2$Row,pheno2$Position)
      gnei2 = nei_coval(geno=geno2,smap=smap2,scale=1,grouping=pheno2$Row) #only neighbors within "Position"

      selfCross = tcrossprod(geno2)
      neiCross = tcrossprod(gnei2)
      
      q = ncol(geno2)
      selfK = (q+selfCross)/(2*(q-1))
      neiK = neiCross/(q-1)
      
      y = as.matrix(pheno2[,k])
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
      
      K12 = tcrossprod(geno2,gnei2) / (q-1)
      K21 = t(K12)
      rownames(K12) = ind_names
      colnames(K12) = ind_names
      rownames(K21) = ind_names
      colnames(K21) = ind_names
      
      covList1 <- rep(list(rep(list(NULL), 2)), 2)
      covList1[[1]][[2]] <- K12
      covList1[[2]][[1]] <- K21
      
      X0 = model.matrix(~poly(pheno2$Row,2)+poly(pheno2$Position,2))
      
      rownames(X0) = ind_names
      
      res = covnei_lmm(geno=geno2,g_nei=gnei2,pheno=pheno2[,k],X0=X0,EM3.res12=EM3res12[[paste0(j,i)]],K1=selfK,K2=neiK,K12=K12,K21=K21,n_core=8L)
      covneilmm_res = c(covneilmm_res,list(res))
    }
  }
  
  names(covneilmm_res) = paste0(rep(levels(factor(pheno$Country)),3),rep(c(2018,2019,2020),each=6))
  saveRDS(covneilmm_res,file=paste0("./output/apple/covneilmm_",k,".rds"))
}

#######
# output analysis

## LL comparison
ll12=c(); ll21=c(); ll22=c()
rho22 = c(); w22 = c()
trait_names = colnames(pheno)[7:30]
for(j in trait_names){
  EM3res12 = readRDS(paste0("./output/apple/",j,"_EM3res12.rds"))
  EM3res21 = readRDS(paste0("./output/apple/",j,"_EM3res21.rds"))
  EM3res22 = readRDS(paste0("./output/apple/",j,"_EM3res22.rds"))
  
  for(i in 1:18){
    if(is.na(EM3res12[[i]])[1]==FALSE){
      ll12 = c(ll12,EM3res12[[i]]$LL)
      if(class(EM3res22[[i]])!="try-error") {
        ll21 = c(ll21,EM3res21[[i]]$LL)
        ll22 = c(ll22,EM3res22[[i]]$LL) 
        rho22 = c(rho22,EM3res22[[i]]$rhosMat[1,2])
        w22 = c(w22,EM3res22[[i]]$weights[2])
      } else { 
        ll21 = c(ll21,NA)
        ll22 = c(ll22,NA) 
        rho22 = c(rho22,NA)
        w22 = c(w22,NA)
          }
    } else {
      ll12 = c(ll12,NA)
      ll21 = c(ll21,NA)
      ll22 = c(ll22,NA)
      rho22 = c(rho22,NA)
      w22 = c(w22,NA)
    }
  }
}
wilcox.test(ll12,ll22,paired=TRUE)
appleLLp = ggplot(data=data.frame(diffLL=ll12-ll22),aes(x=diffLL)) + geom_histogram(fill=grey(0.5,0.5)) + geom_boxplot(colour="black",size=0.5) + theme_bw() + xlab("delta_LL")
ggsave(appleLLp,filename="../figs/appleLLp.pdf",height=6,width=6)

# export for Table S2
which(pchisq(-2*(ll21-ll22),1,lower.tail=FALSE)<0.05)
which(rho22[which(pchisq(-2*(ll21-ll22),1,lower.tail=FALSE)<0.05)]<0)
var_res[which(pchisq(-2*(ll21-ll22),1,lower.tail=FALSE)<0.05)[20],]

chi=-2*(ll21-ll22)
pval = pchisq(-2*(ll21-ll22),1,lower.tail=FALSE)
write.csv(data.frame(w22,rho22,ll22,chi,pval),"../tmp.csv")


## SE for Table S2
cov_res = c()
trait_names = colnames(pheno)[7:30]
for(j in trait_names){
  EM3res11 = readRDS(paste0("./output/apple/",j,"_EM3res11.rds"))
  EM3res12 = readRDS(paste0("./output/apple/",j,"_EM3res12.rds"))
  EM3res12se = readRDS(paste0("./output/apple/",j,"_EM3res12se.rds"))
  
  covlist = c()
  for(i in 1:18) {
    if(is.na(EM3res12[[i]])[1]==FALSE){
      chi2 = -2*(EM3res11[[i]]$LL - EM3res12[[i]]$LL)
      pchi = pchisq(chi2,1,lower.tail=FALSE)
      if(class(EM3res12se[[i]])!="try-error"){
        covlist = rbind(covlist,c(substr(names(EM3res12)[i],1,3),
                                  substr(names(EM3res12)[i],4,7),
                                  EM3res12se[[i]]$parameterTable[1,3],EM3res12se[[i]]$parameterTable[1,4],
                                  EM3res12se[[i]]$parameterTable[3,3],EM3res12se[[i]]$parameterTable[3,4],
                                  EM3res12se[[i]]$parameterTable[2,3],EM3res12se[[i]]$parameterTable[2,4],
                                  EM3res12se[[i]]$parameterTable[4,3],EM3res12se[[i]]$parameterTable[4,4],
                                  chi2,pchi))
      } else {
        covlist = rbind(covlist,c(substr(names(EM3res12)[i],1,3),
                                  substr(names(EM3res12)[i],4,7),
                                  EM3res12[[i]]$Vu*EM3res12[[i]]$weights[1],NA,
                                  EM3res12[[i]]$Vu*EM3res12[[i]]$weights[2],NA,
                                  EM3res12[[i]]$rhosMat[1,2],NA,
                                  EM3res12[[i]]$Ve,NA,
                                  chi2,pchi))
      }
      
    } else {
      covlist = rbind(covlist,c(substr(names(EM3res12)[i],1,3),
                                substr(names(EM3res12)[i],4,7),
                                rep(NA,8),
                                NA,NA))
    }
  }
  colnames(covlist) = c("site","year","sigma1","sigma1_se","sigma2","sigma2_se","sigma12","sigma12_se","sigmae","sigmae_se","chisq","pval")
  covlist = as.data.frame(covlist)
  
  cov_res = rbind(cov_res,covlist)
}
var_res = data.frame(trait=rep(trait_names,each=18),cov_res)
write.csv(var_res,"../apple_tableS2.csv")

## PVEs for figures
cov_res = c()
trait_names = colnames(pheno)[7:30]
trait_names = trait_names
for(j in trait_names){
  EM3res11 = readRDS(paste0("./output/apple/",j,"_EM3res11.rds"))
  EM3res12 = readRDS(paste0("./output/apple/",j,"_EM3res12.rds"))
  
  covlist = c()
  for(i in 1:18) {
    if(is.na(EM3res12[[i]])[1]==FALSE){
      chi2 = -2*(EM3res11[[i]]$LL - EM3res12[[i]]$LL)
      #print(pchisq(chi2,1,lower.tail=FALSE))
      covlist = rbind(covlist,c(substr(names(EM3res12)[i],1,3),
                                substr(names(EM3res12)[i],4,7),
                                EM3res12[[i]]$Vu,
                                EM3res12[[i]]$Ve,
                                EM3res12[[i]]$weights,
                                EM3res12[[i]]$rhosMat[1,2],
                                chi2))
    } else {
      covlist = rbind(covlist,c(substr(names(EM3res12)[i],1,3),
                                substr(names(EM3res12)[i],4,7),
                                NA,
                                NA,
                                NA,NA,
                                NA,
                                NA))
    }
  }
  colnames(covlist) = c("site","year","Vu","Ve","w1","w2","rho12","chisq")
  covlist = as.data.frame(covlist)
  
  cov_res = c(cov_res,list(covlist))
}
names(cov_res) = trait_names

for(i in 1:length(cov_res)){
  h2 = as.numeric(cov_res[[i]]$Vu)/(as.numeric(cov_res[[i]]$Vu) + as.numeric(cov_res[[i]]$Ve))
  DGE=h2*as.numeric(cov_res[[i]]$w1)
  IGE=h2*as.numeric(cov_res[[i]]$w2)
  cov=h2*as.numeric(cov_res[[i]]$rho12)*sqrt(as.numeric(cov_res[[i]]$w1)*as.numeric(cov_res[[i]]$w2))
  site_year = paste0(cov_res[[i]]$site,"_",cov_res[[i]]$year)
  cov_res[[i]] = data.frame(cov_res[[i]],h2,DGE,IGE,cov,site_year)
}

pvep = c()
for(i in 1:length(cov_res)){
  p = cov_res[[i]] %>%
    gather(c("DGE","IGE","cov"),key="factor",value="PVE") %>%
    ggplot(aes(x=site_year, y=PVE, fill=factor)) + ylim(-0.25,1) + labs(subtitle=names(cov_res)[[i]]) +
    geom_bar(stat="identity", position=position_dodge()) + coord_flip() + theme_minimal()
  pvep = c(pvep,list(p))
}

m = (pvep[[1]] | pvep[[2]] | pvep[[3]] | pvep[[4]]) / (pvep[[6]] | pvep[[7]] | pvep[[8]] | pvep[[9]]) / (pvep[[10]] | pvep[[11]] | pvep[[12]] | pvep[[13]]) / (pvep[[14]] | pvep[[15]] | pvep[[16]] | pvep[[17]]) / (pvep[[18]] | pvep[[19]] | pvep[[20]] | pvep[[21]]) / (pvep[[22]] | plot_spacer() | plot_spacer() | plot_spacer()) + plot_annotation(tag_levels = "a")
ggsave(m,filename="../apple_pve.pdf",width=12,height=20) # export a supplementary figure

# LMM regressing PVE~cov~ on the mean trunk diameter
mean_diameter = aggregate(Trunk_diameter~Country+Year,pheno,mean)[,3]
summary(lme4::lmer(cov~mean_diameter+(1|site)+(1|year),data=cov_res$Trunk_increment))
summary(lme4::lmer(cov~mean_diameter+(1|site)+(1|year),data=cov_res$Trunk_diameter))
summary(lme4::lmer(cov~mean_diameter+(1|site)+(1|year),data=cov_res$Flowering_intensity))

pchi1 = pchisq(as.numeric(cov_res$Trunk_increment$chisq),1,lower.tail=FALSE)
p11 = data.frame(cov_res$Trunk_increment,mean_diameter=mean_diameter,pLL=pchi1) %>%
  ggplot(aes(x=mean_diameter,y=cov)) + geom_point(color=1+as.numeric(pchi1<0.05)) + 
  geom_abline(intercept=0.173584,slope=-0.007585) + 
  labs(subtitle="Trunk_increment") + xlab("Mean trunk diameter (cm)") + ylab(expression(PVE[cov])) +
  theme_bw()

pchi2 = pchisq(as.numeric(cov_res$Trunk_diameter$chisq),1,lower.tail=FALSE)
p21 = data.frame(cov_res$Trunk_diameter,mean_diameter=mean_diameter,pLL=pchi2) %>%
  ggplot(aes(x=mean_diameter,y=cov)) + geom_point(color=1+as.numeric(pchi2<0.05)) + 
  geom_abline(intercept=0.052754,slope=-0.004215) + 
  labs(subtitle="Trunk_diameter") + xlab("Mean trunk diameter (cm)") + ylab(expression(PVE[cov])) +
  theme_bw()

pchi3 = pchisq(as.numeric(cov_res$Flowering_intensity$chisq),1,lower.tail=FALSE)
p31 = data.frame(cov_res$Flowering_intensity,mean_diameter=mean_diameter,pLL=pchi3) %>%
  ggplot(aes(x=mean_diameter,y=cov)) + geom_point(color=1+as.numeric(pchi3<0.05)) + 
  geom_abline(intercept=0.0288035,slope=0.0002315,lty=2) + 
  labs(subtitle="Flowering_intensity") + xlab("Mean trunk diameter (cm)") + ylab(expression(PVE[cov])) +
  theme_bw()

# create all Manhattan and QQ plots
# count no. of significant SNPs for Table S2
pdf()
nsnps = c()
gifs = c()
for(k in trait_names) { print(k)
  out = readRDS(file=paste0("./output/apple/covneilmm_",k,".rds"))
  TF = is.na(out)
  for(l in 1:18){
    if(TF[l]!=TRUE){
      nsnps = c(nsnps,sum(p.adjust(out[[l]]$p2,method="BH")<0.05))
      gifs = c(gifs,median(log10(out[[l]]$p2))/median(log10(ppoints(length(out[[l]]$p2)))))
      gaston::manhattan(data.frame(chr=geno$gmap$chr,pos=geno$gmap$pos,p=out[[l]]$p2),
                        main=k,sub=names(out)[[l]],las=1)
      abline(h=-log10(0.05/length(out[[l]]$p2)))
      gaston::qqplot.pvalues(p=out[[l]]$p2,main=k,sub=names(out)[[l]],las=1)
    } else {
      nsnps = c(nsnps,NA)
      gifs = c(gifs,NA)
      plot(1,1,type="n",main=k,sub=names(out)[[l]],las=1)
    }
  }
}
dev.off()

# no. of effective samples for Table S2
n_ind = c()
for(k in colnames(pheno)[7:30]){
  for(i in c(2018,2019,2020)){
    for(j in levels(factor(pheno$Country))){
      subpheno = filter(pheno,Country==j&Year==i)
      n_kij = sum(is.na(subpheno[,k])==FALSE)
      n_ind = c(n_ind,n_kij)
    }
  }
}
write.csv(cbind(nsnps,gifs,n_ind),"../tmp.csv") # export for Table S2

# trunk increment, BEL2019, significant on chr7
gwasnps = data.frame(trait = rep(trait_names,each=length(out)),trial = rep(names(out),length(trait_names)),nsnps)
na.omit(gwasnps)

ggMan = function(chr,pos,p) {
  chr_rep = table(chr)
  alphas = 1-(chr/2 - chr %/% 2)
  cols = rgb(0,0,alphas,alphas)
  x = coord(chr,pos)
  y = -log10(p)
  fdr = CalcThreshold(input=data.frame(name=c(1:length(p)),chr=chr,pos=pos,p=-log10(p)))
  
  man = ggplot(NULL,aes(x=x$coord,y=y)) + geom_point(colour=cols) + theme_classic() + 
    scale_x_continuous(name="Chromosomes", breaks=x$tic, labels=names(chr_rep)) +
    ylab(expression(-log[10]*(italic(p)[2]))) + geom_hline(yintercept=fdr,lty=1,colour="red")
  return(man)
}

pm = ggMan(chr=geno$gmap$chr,pos=geno$gmap$pos,p=out$BEL2019$p2)
pm = pm + labs(subtitle="BEL2019: Trunk_increment")

# export Figure 4
p1234 = (pvep[[5]] | p31) / (pvep[[23]] | p21) / (pvep[[24]] | p11) / pm + plot_layout(widths = c(2,1), heights = c(1,1,1,0.5)) + plot_annotation(tag_levels = "a")
ggsave(p1234,filename="../fig4_apple.jpg",width=9,height=10,dpi=600)
# ggsave(p1234,filename="../fig4_apple.pdf",width=9,height=10)

