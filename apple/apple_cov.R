#####################
# apple REFPOP data

library(rNeighborGWAS)
library(gaston)
library(RAINBOWR)
library(lme4)
library(tidyverse)
library(patchwork)
source("coord.R")

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
for(k in colnames(pheno)[7:29]){ #error in "Russet_freq_cheek CHR2019, optimizer = "nlminb", fine"
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
      
      X0 = matrix(1,nrow=nrow(pheno2))
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

#####################
# GWAS for all traits
for(k in colnames(pheno)[22]){
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
      
      X0 = matrix(1,nrow=nrow(pheno2))
      rownames(X0) = ind_names
      
      res = covnei_lmm(geno=geno2,g_nei=gnei2,pheno=pheno2[,k],X0=X0,EM3.res12=EM3res12[[paste0(j,i)]],K1=selfK,K2=neiK,K12=K12,K21=K21,n_core=6L)
      covneilmm_res = c(covneilmm_res,list(res))
    }
  }
  
  names(covneilmm_res) = paste0(rep(levels(factor(pheno$Country)),3),rep(c(2018,2019,2020),each=6))
  saveRDS(covneilmm_res,file=paste0("./output/apple/covneilmm_",k,".rds"))
}

#######
# output analysis

## PVEs for Table S2
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
ggsave(m,filename="../apple_pve.pdf",width=12,height=20)

# LMM regressing PVE~cov~ on the mean trunk diameter
mean_diameter = aggregate(Trunk_diameter~Country+Year,pheno,mean)[,3]
summary(lme4::lmer(cov~mean_diameter+(1|site)+(1|year),data=cov_res$Trunk_increment))
summary(lme4::lmer(cov~mean_diameter+(1|site)+(1|year),data=cov_res$Trunk_diameter))
summary(lme4::lmer(cov~mean_diameter+(1|site)+(1|year),data=cov_res$Flowering_intensity))

pchi1 = pchisq(as.numeric(cov_res$Trunk_increment$chisq),1,lower.tail=FALSE)
p11 = data.frame(cov_res$Trunk_increment,mean_diameter=mean_diameter,pLL=pchi1) %>%
  ggplot(aes(x=mean_diameter,y=cov)) + geom_point(color=1+as.numeric(pchi1<0.05)) + 
  geom_abline(intercept=0.132052,slope=-0.006452) + 
  labs(subtitle="Trunk_increment") + xlab("Mean trunk diameter (cm)") + ylab(expression(PVE[cov])) +
  theme_bw()

pchi2 = pchisq(as.numeric(cov_res$Trunk_diameter$chisq),1,lower.tail=FALSE)
p21 = data.frame(cov_res$Trunk_diameter,mean_diameter=mean_diameter,pLL=pchi2) %>%
  ggplot(aes(x=mean_diameter,y=cov)) + geom_point(color=1+as.numeric(pchi2<0.05)) + 
  geom_abline(intercept=-0.034339,slope=-0.002247,lty=2) + 
  labs(subtitle="Trunk_diameter") + xlab("Mean trunk diameter (cm)") + ylab(expression(PVE[cov])) +
  theme_bw()

pchi3 = pchisq(as.numeric(cov_res$Flowering_intensity$chisq),1,lower.tail=FALSE)
p31 = data.frame(cov_res$Flowering_intensity,mean_diameter=mean_diameter,pLL=pchi3) %>%
  ggplot(aes(x=mean_diameter,y=cov)) + geom_point(color=1+as.numeric(pchi3<0.05)) + 
  geom_abline(intercept=0.005555,slope=0.000791,lty=2) + 
  labs(subtitle="Flowering_intensity") + xlab("Mean trunk diameter (cm)") + ylab(expression(PVE[cov])) +
  theme_bw()

# create all Manhattan and QQ plots
# count no. of sifnificant SNPs for Table S2
pdf()
nsnps = c()
for(k in trait_names) { print(k)
  out = readRDS(file=paste0("./output/apple/covneilmm_",k,".rds"))
  TF = is.na(out)
  for(l in 1:18){
    if(TF[l]!=TRUE){
      nsnps = c(nsnps,sum(p.adjust(out[[l]]$p2,method="BH")<0.05))
      gaston::manhattan(data.frame(chr=geno$gmap$chr,pos=geno$gmap$pos,p=out[[l]]$p2),
                        main=k,sub=names(out)[[l]],las=1)
      abline(h=-log10(0.05/length(out[[l]]$p2)))
      gaston::qqplot.pvalues(p=out[[l]]$p2,main=k,sub=names(out)[[l]],las=1)
    } else {
      nsnps = c(nsnps,NA)
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
    ylab(expression(-log[10]*(italic(p)[2]))) + geom_hline(yintercept=fdr,lty=2,colour="black")
  return(man)
}

pm = ggMan(chr=geno$gmap$chr,pos=geno$gmap$pos,p=out$BEL2019$p2)
pm = pm + labs(subtitle="BEL2019: Trunk_increment")

p1234 = (pvep[[5]] | p31) / (pvep[[23]] | p21) / (pvep[[24]] | p11) / pm + plot_layout(widths = c(2,1), heights = c(1,1,1,0.5)) + plot_annotation(tag_levels = "a")
ggsave(p1234,filename="../fig4_apple.jpg",width=9,height=10,dpi=300)

