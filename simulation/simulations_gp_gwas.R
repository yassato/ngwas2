set.seed(1234)

library(gaston)
library(vcfR)
library(rNeighborGWAS)
library(pROC)
library(RAINBOWR)
library(tidyverse)
library(patchwork)
source("covnei_lmm.R")
source("simulation/simPhenoWithCov.R")

#################################
# for GP and variance simulations
iter_i = function(i){
  f_path = paste0("./simulation/",i,"_STVS.vcf.gz")
  d0 = read.vcf(f_path)
  d = select.snps(d0, d0@snps$maf>0.01)
  
  ann = read.vcfR(f_path)
  ann = ann[which(d0@snps$maf>0.01),]
  
  geno = gaston2neiGWAS(d)
  
  ngrid = 40
  g = geno$geno[sample(1:nrow(geno$geno),ngrid^2,replace=FALSE),]
  g[g==0] = 1 # assume complete dominance
  
  smap = cbind(rep(1:(nrow(g)/ngrid),each=nrow(g)/ngrid),
               rep(1:(nrow(g)/ngrid),nrow(g)/ngrid)) # 40 x 40 map
  g_nei = nei_coval(g,smap,scale=1.42)
  
  selfCross = tcrossprod(g)
  neiCross = tcrossprod(g_nei)
  
  q = ncol(g)
  selfK = (q+selfCross)/(2*(q-1))
  neiK = neiCross/(q-1)
  
  ind_names = paste0("ind",1:nrow(g))
  Z = diag(nrow(g))
  rownames(Z) = ind_names
  colnames(Z) = ind_names
  rownames(selfK) = ind_names
  colnames(selfK) = ind_names
  rownames(neiK) = ind_names
  colnames(neiK) = ind_names
  ZETA1 <- list(A = list(Z = Z, K = selfK),
                Int = list(Z = Z, K = neiK))
  
  K12 = tcrossprod(g,g_nei) / (q-1)
  K21 = tcrossprod(g_nei,g) / (q-1)
  rownames(K12) = ind_names
  colnames(K12) = ind_names
  rownames(K21) = ind_names
  colnames(K21) = ind_names
  
  covList1 <- rep(list(rep(list(NULL), 2)), 2)
  covList1[[1]][[2]] <- K12
  covList1[[2]][[1]] <- K21
  
  X0 = as.matrix(rep(1,nrow(g)))
  rownames(X0) = ind_names
  
  ### Parameters
  herit <- 0.6
  weights <- c(0.5, 0.5)
  rhosMat <- matrix(
    data = c(1, 0.6, # changing 0.6 simulates the other cases of \rho_{12}
             0.6, 1),
    nrow = 2, ncol = 2
  )
  X <- NULL
  beta <- NULL
  
  simRes <- simPhenoWithCov(ZETA = ZETA1, covList = covList1,
                            herit = herit, weights = weights,
                            rhosMat = rhosMat, X = X, beta = beta)
  y <- simRes$y
  
  # keep 1/10 for test data
  test_y = sample(1:length(y),ngrid^2/10)
  org_y = y
  y[test_y] = NA
  
  PVEnei.res = calc_PVEnei(pheno=org_y,geno=g,smap=smap,scale_seq=1.42)
  
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
  
  ysim = data.frame(org_y,y)
  res = list(PVEnei=PVEnei.res,EM3res11=EM3.res11,EM3res12=EM3.res12,ysim=ysim)
  saveRDS(res,file=paste0("./output/simulation/",i,"_EM3cov_pos.rds"))
  return(NULL)
}

set.seed(1234)
res = parallel::mcmapply(iter_i,1:2,mc.cores=2L)

# variance component check
rho0 = c()
w01 = c()
w02 = c()
chi0 = c()
for(i in 1:30){
  out = readRDS(paste0("./output/simulation/",i,"_EM3cov_null.rds"))
  rho0 = c(rho0,out$EM3res12$rhosMat[1,2])
  w01 = rbind(w01,out$EM3res11$weights)
  w02 = rbind(w02,out$EM3res12$weights)
  chi0 = c(chi0,-2*(out$EM3res11$LL-out$EM3res12$LL))
}
pchisq(chi0,1,lower.tail=FALSE) # all non-significant

mean((abs(rho0+0)))
p01 = ggplot(data=NULL,aes(x="",y=rho0)) + geom_boxplot() + geom_jitter() + 
  geom_hline(yintercept=0,lwd=1,lty=2) + ylab(expression(hat(rho)[12])) + xlab("cov") + theme_bw()

w0 = data.frame(noncov=w01[,2],cov=w02[,2])
mean(abs(w0[,2]-0.5))
p02 = w0 %>%
  gather(c("noncov","cov"),key="model",value="w2") %>%
  ggplot(aes(x=model,y=w2)) + geom_boxplot() + geom_jitter() + ylim(0,1) + 
  geom_hline(yintercept=0.5,lwd=1,lty=2) + ylab(expression(hat(w)[2])) + theme_bw()


for(i in 1:30){
  out = readRDS(paste0("./output/simulation/",i,"_EM3cov_pos.rds"))
  PVE0 = out$PVEnei$PVEself[1]
  PVE11 = out$EM3res11$Vu / (out$EM3res11$Vu + out$EM3res11$Ve)
  PVE12 = out$EM3res12$Vu / (out$EM3res12$Vu + out$EM3res12$Ve)
  print(c(PVE0,PVE11,PVE12))
}

rhon = c()
wn11 = c()
wn12 = c()
for(i in 1:30){
  out = readRDS(paste0("./output/simulation/",i,"_EM3cov_neg.rds"))
  rhon = c(rhon,out$EM3res12$rhosMat[1,2])
  wn11 = rbind(wn11,out$EM3res11$weights)
  wn12 = rbind(wn12,out$EM3res12$weights)
}
mean((abs(rhon+0.6)))
p1 = ggplot(data=NULL,aes(x="",y=rhon)) + geom_boxplot() + geom_jitter() + ylim(NA,0) +
  geom_hline(yintercept=-0.6,lwd=1,lty=2) + ylab(expression(hat(rho)[12])) + xlab("cov") + theme_bw()

mean(abs(wn[,2]-0.5))
wn = data.frame(noncov=wn11[,2],cov=wn12[,2])
p2 = wn %>%
  gather(c("noncov","cov"),key="model",value="w2") %>%
  ggplot(aes(x=model,y=w2)) + geom_boxplot() + geom_jitter() + ylim(0,1) + 
  geom_hline(yintercept=0.5,lwd=1,lty=2) + ylab(expression(hat(w)[2])) + theme_bw()

rhop = c()
wp11 = c()
wp12 = c()
for(i in 1:30){
  out = readRDS(paste0("./output/simulation/",i,"_EM3cov_pos.rds"))
  rhop = c(rhop,out$EM3res12$rhosMat[1,2])
  wp11 = rbind(wp11,out$EM3res12$weights)
  wp12 = rbind(wp12,out$EM3res12$weights)
}
mean((abs(rhop-0.6)))
p3 = ggplot(data=NULL,aes(x="",y=rhop)) + geom_boxplot() + geom_jitter() + ylim(0,NA) +
  geom_hline(yintercept=0.6,lwd=1,lty=2) + ylab(expression(hat(rho)[12])) + xlab("cov") + theme_bw()

mean(abs(wp[,2]-0.5))
wp = data.frame(noncov=wp11[,2],cov=wp12[,2])
p4 = wp %>%
  gather(c("noncov","cov"),key="model",value="w2") %>%
  ggplot(aes(x=model,y=w2)) + geom_boxplot() + geom_jitter() + ylim(0,1) + 
  geom_hline(yintercept=0.5,lwd=1,lty=2) + ylab(expression(hat(w)[2])) + theme_bw()


# GP output
ps = c()
r2s = c()
r11s = c(); r12s = c()
mae11s = c(); mae12s = c()
for(i in 1:30){
  out = readRDS(paste0("./output/simulation/",i,"_EM3cov_null.rds"))
  ps = c(ps,pchisq(-2*(out$EM3res11$LL-out$EM3res12$LL),1,lower.tail=FALSE))
  r2s = c(r2s,cor(out$EM3res12$y.pred,out$EM3res11$y.pred)^2)
  r11s = c(r11s,cor(out$ysim$org_y[is.na(out$ysim$y)==FALSE],out$EM3res11$y.pred[is.na(out$ysim$y)==FALSE])^2)
  r12s = c(r12s,cor(out$ysim$org_y[is.na(out$ysim$y)==FALSE],out$EM3res12$y.pred[is.na(out$ysim$y)==FALSE])^2)
  mae11s = c(mae11s,mean(abs(out$ysim$org_y[is.na(out$ysim$y)==FALSE]-out$EM3res11$y.pred[is.na(out$ysim$y)==FALSE])))
  mae12s = c(mae12s,mean(abs(out$ysim$org_y[is.na(out$ysim$y)==FALSE]-out$EM3res12$y.pred[is.na(out$ysim$y)==FALSE])))
  print(i)
}

r2_0 = data.frame(noncov=r11s,cov=r12s)
p05 = r2_0 %>%
  gather(c("noncov","cov"),key="model",value="R2") %>%
  ggplot(aes(x=model,y=R2)) + geom_boxplot() + geom_jitter() + ylab(expression(R^2)) + ylim(0.5,1) + theme_bw()

mae_0 = data.frame(noncov=mae11s,cov=mae12s)
p06 = mae_0 %>%
  gather(c("noncov","cov"),key="model",value="MAE") %>%
  ggplot(aes(x=model,y=MAE)) + geom_boxplot() + geom_jitter() + ylim(0,0.45) + theme_bw()

ps = c()
r2s = c()
r11s = c(); r12s = c()
mae11s = c(); mae12s = c()
for(i in 1:30){
  out = readRDS(paste0("./output/simulation/",i,"_EM3cov_neg.rds"))
  ps = c(ps,pchisq(-2*(out$EM3res11$LL-out$EM3res12$LL),1,lower.tail=FALSE))
  r2s = c(r2s,cor(out$EM3res12$y.pred,out$EM3res11$y.pred)^2)
  r11s = c(r11s,cor(out$ysim$org_y[is.na(out$ysim$y)==FALSE],out$EM3res11$y.pred[is.na(out$ysim$y)==FALSE])^2)
  r12s = c(r12s,cor(out$ysim$org_y[is.na(out$ysim$y)==FALSE],out$EM3res12$y.pred[is.na(out$ysim$y)==FALSE])^2)
  mae11s = c(mae11s,mean(abs(out$ysim$org_y[is.na(out$ysim$y)==FALSE]-out$EM3res11$y.pred[is.na(out$ysim$y)==FALSE])))
  mae12s = c(mae12s,mean(abs(out$ysim$org_y[is.na(out$ysim$y)==FALSE]-out$EM3res12$y.pred[is.na(out$ysim$y)==FALSE])))
  print(i)
}

r2_neg = data.frame(noncov=r11s,cov=r12s)
p5 = r2_neg %>%
  gather(c("noncov","cov"),key="model",value="R2") %>%
  ggplot(aes(x=model,y=R2)) + geom_boxplot() + geom_jitter() + ylab(expression(R^2)) + ylim(0.5,1) + theme_bw()

mae_neg = data.frame(noncov=mae11s,cov=mae12s)
p6 = mae_neg %>%
  gather(c("noncov","cov"),key="model",value="MAE") %>%
  ggplot(aes(x=model,y=MAE)) + geom_boxplot() + geom_jitter() + ylim(0,0.45) + theme_bw()


ps = c()
r2s = c()
r11s = c(); r12s = c()
mae11s = c(); mae12s = c()
for(i in 1:30){
  out = readRDS(paste0("./output/simulation/",i,"_EM3cov_pos.rds"))
  ps = c(ps,pchisq(-2*(out$EM3res11$LL-out$EM3res12$LL),1,lower.tail=FALSE))
  r2s = c(r2s,cor(out$EM3res12$y.pred,out$EM3res11$y.pred)^2)
  r11s = c(r11s,cor(out$ysim$org_y[is.na(out$ysim$y)==FALSE],out$EM3res11$y.pred[is.na(out$ysim$y)==FALSE])^2)
  r12s = c(r12s,cor(out$ysim$org_y[is.na(out$ysim$y)==FALSE],out$EM3res12$y.pred[is.na(out$ysim$y)==FALSE])^2)
  mae11s = c(mae11s,mean(abs(out$ysim$org_y[is.na(out$ysim$y)==FALSE]-out$EM3res11$y.pred[is.na(out$ysim$y)==FALSE])))
  mae12s = c(mae12s,mean(abs(out$ysim$org_y[is.na(out$ysim$y)==FALSE]-out$EM3res12$y.pred[is.na(out$ysim$y)==FALSE])))
  print(i)
}

r2_pos = data.frame(noncov=r11s,cov=r12s)
p7 = r2_pos %>%
  gather(c("noncov","cov"),key="model",value="R2") %>%
  ggplot(aes(x=model,y=R2)) + geom_boxplot() + geom_jitter() + ylab(expression(R^2)) + ylim(0.5,1) + theme_bw()

mae_pos = data.frame(noncov=mae11s,cov=mae12s)
p8 = mae_pos %>%
  gather(c("noncov","cov"),key="model",value="MAE") %>%
  ggplot(aes(x=model,y=MAE)) + geom_boxplot() + geom_jitter() + ylim(0,0.45) + theme_bw()

# R2 between two models
r2s = c()
for(i in 1:30){
  out = readRDS(paste0("./output/simulation/",i,"_EM3cov_pos.rds"))
  r2s = c(r2s,cor(out$EM3res12$y.pred,out$EM3res11$y.pred)^2)
  
  out = readRDS(paste0("./output/simulation/",i,"_EM3cov_pos03.rds"))
  r2s = c(r2s,cor(out$EM3res12$y.pred,out$EM3res11$y.pred)^2)
  
  out = readRDS(paste0("./output/simulation/",i,"_EM3cov_null.rds"))
  r2s = c(r2s,cor(out$EM3res12$y.pred,out$EM3res11$y.pred)^2)
  
  out = readRDS(paste0("./output/simulation/",i,"_EM3cov_neg03.rds"))
  r2s = c(r2s,cor(out$EM3res12$y.pred,out$EM3res11$y.pred)^2)
  
  out = readRDS(paste0("./output/simulation/",i,"_EM3cov_neg.rds"))
  r2s = c(r2s,cor(out$EM3res12$y.pred,out$EM3res11$y.pred)^2)
  
  print(i)
}

r2 = data.frame(rho=rep(c(0.6,0.3,0,-0.3,-0.6), 30), R2 = r2s)
ps1 = r2 %>%
  ggplot(aes(x=factor(rho),y=R2)) + geom_boxplot() + geom_jitter() + 
  xlab(expression(rho[12])) + ylab(expression(R^2)) + ylim(0.9,1) + theme_bw()

ggsave(ps1,filename="../BLUP_R2.pdf",width = 4,height = 6)


##########################
#for GWAS simulations
for(i in 1:30){
  f_path = paste0("./simulation/",i,"_STVS.vcf.gz")
  d0 = read.vcf(f_path)
  d = select.snps(d0, d0@snps$maf>0.01)
  
  ann = read.vcfR(f_path)
  ann = ann[which(d0@snps$maf>0.01),]
  
  geno = gaston2neiGWAS(d)
  
  ngrid = 40
  g = geno$geno[sample(1:nrow(geno$geno),ngrid^2,replace=FALSE),]
  g[g==0] = 1 # assume complete dominance
  
  smap = cbind(rep(1:(nrow(g)/ngrid),each=nrow(g)/ngrid),
               rep(1:(nrow(g)/ngrid),nrow(g)/ngrid)
  )#40 x 40 map
  g_nei = nei_coval(g,smap,scale=1.42)
  
  selfCross = tcrossprod(g)
  neiCross = tcrossprod(g_nei)
  
  q = ncol(g)
  selfK = (q+selfCross)/(2*(q-1))
  neiK = neiCross/(q-1)
  
  ind_names = paste0("ind",1:nrow(g))
  Z = diag(nrow(g))
  rownames(Z) = ind_names
  colnames(Z) = ind_names
  rownames(selfK) = ind_names
  colnames(selfK) = ind_names
  rownames(neiK) = ind_names
  colnames(neiK) = ind_names
  ZETA1 <- list(A = list(Z = Z, K = selfK),
                Int = list(Z = Z, K = neiK))
  
  K12 = tcrossprod(g,g_nei) / (q-1)
  K21 = tcrossprod(g_nei,g) / (q-1)
  rownames(K12) = ind_names
  colnames(K12) = ind_names
  rownames(K21) = ind_names
  colnames(K21) = ind_names
  
  covList1 <- rep(list(rep(list(NULL), 2)), 2)
  covList1[[1]][[2]] <- K12
  covList1[[2]][[1]] <- K21
  
  X0 = as.matrix(rep(1,nrow(g)))
  rownames(X0) = ind_names
  
  ### Parameters
  herit <- 0.6
  weights <- c(0.5, 0.5)
  rhosMat <- matrix(
    data = c(1, 0.6, # changing 0.6 into -0.6 simulates the case of negative \rho_12
             0.6, 1),
    nrow = 2, ncol = 2
  )
  
  beta1 = rep(0,nrow(d@snps))
  n_causal1 = length(grep("MT=3",ann@fix[,"INFO"]))
  beta1[grep("MT=3",ann@fix[,"INFO"])] = 0.1
  
  beta2 = rep(0,nrow(d@snps))
  n_causal2 = 10
  s1 = sample(grep("MT=3",ann@fix[,"INFO"]),n_causal2) 
  s2 = sample(grep("MT=3",ann@fix[,"INFO"]),n_causal2) 
  beta2[s2] = 0.1 #PFDS
  beta2[s1] = - 0.1 # NFDS
  beta1[s1] = -0.1 # make them correlated
  
  g_12 = g*g_nei
  beta12 = rep(0,nrow(d@snps))
  n_causal12 = 5
  s12 = sample(s1,n_causal12)
  s11 = sample(s2,n_causal12)
  beta12[s12] = 0.1 #PFDS
  beta12[s11] = - 0.1 # NFDS
  
  X <- cbind(g,g_nei,g_12)
  beta <- c(beta1,beta2,beta12)
  
  simRes <- simPhenoWithCov(ZETA = ZETA1, covList = covList1,
                            herit = herit, weights = weights,
                            rhosMat = rhosMat, X = X, beta = beta)
  y <- simRes$y
  
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
  
  
  res11 = covnei_lmm11(g,g_nei,pheno=y,X0=X0,EM3.res11=EM3.res11,n_core=1L)
  res12 = covnei_lmm12(g,g_nei,pheno=y,X0=X0,EM3.res12=EM3.res12,n_core=1L)
  
  res.roc = list(beta2=beta2,beta12=beta12,res11=res11,res12=res12,EM3res11=EM3.res11,EM3res12=EM3.res12)
  saveRDS(res.roc,file=paste0("./output/simulation/",i,"_gwas_pos.rds"))
}

# GWAS output
aucs = c(); scores = c(); pval = c()
pred11 = c(); resp11 = c()
pred12 = c(); resp12 = c()
for(i in 1:30){ print(i)
  out = readRDS(file=paste0("./output/simulation/",i,"_gwas_neg.rds"))
  print(out$EM3res12$rhosMat[1,2])
  aucs = rbind(aucs,c(out$roc11$auc,out$roc12$auc))
  scores = rbind(scores,c(mean(out$roc11$predictor[out$roc11$response==1]),
                   mean(out$roc12$predictor[out$roc12$response==1])))
  
  pval = c(pval,pchisq(-2*(out$EM3res11$LL-out$EM3res12$LL),1,lower.tail=FALSE))
  
  pred11 = c(pred11,out$roc11$predictor)
  resp11 = c(resp11,out$roc11$response)
  
  pred12 = c(pred12,out$roc12$predictor)
  resp12 = c(resp12,out$roc12$response)
}


aucs = c(); scores = c(); power = c()
for(i in 1:30){
  out = readRDS(file=paste0("./output/simulation/",i,"_gwas_pos.rds"))
  ans = as.numeric(out$beta12==0&out$beta2!=0)
  p11 = -log10(out$res11$p2)
  p12 = -log10(out$res12$p2)
  roc11 = roc(ans~p11,smooth=TRUE)
  roc12 = roc(ans~p12,smooth=TRUE)
  aucs = rbind(aucs,c(roc11$auc,roc12$auc))
  power = rbind(aucs,c(max(roc11$sensitivities[roc11$specificities>0.95]),max(roc12$sensitivities[roc12$specificities>0.95])))
  scores = rbind(scores,c(median(p11[ans!=0]),median(p12[ans!=0])))
}


colnames(aucs) = c("noncov","cov")
auc_neg = as.data.frame(aucs)
p9 = auc_neg %>%
  gather(c("noncov","cov"),key="model",value="AUC") %>%
  ggplot(aes(x=model,y=AUC)) + geom_boxplot() + geom_jitter() + ylim(NA,1.0) + theme_bw()

colnames(power) = c("noncov","cov")
power_neg = as.data.frame(power)
p10 = power_neg %>%
  gather(c("noncov","cov"),key="model",value="sensitivity") %>%
  ggplot(aes(x=model,y=sensitivity)) + geom_boxplot() + geom_jitter() + ylim(NA,1.0) + ylab("sensitivity at FDR=0.05") + theme_bw()

colnames(aucs) = c("noncov","cov")
auc_pos = as.data.frame(aucs)
p11 = auc_pos %>%
  gather(c("noncov","cov"),key="model",value="AUC") %>%
  ggplot(aes(x=model,y=AUC)) + geom_boxplot() + geom_jitter() + ylim(NA,1.0) + theme_bw()

colnames(power) = c("noncov","cov")
power_pos = as.data.frame(power)
p12 = power_pos %>%
  gather(c("noncov","cov"),key="model",value="sensitivity") %>%
  ggplot(aes(x=model,y=sensitivity)) + geom_boxplot() + geom_jitter() + ylim(NA,1.0) + ylab("sensitivity at FDR=0.05") + theme_bw()


allp=((p3+labs(title="cov+") | p4 | p7 | p8 | p11 | p12) / (p1+labs(title="cov-") | p2 | p5 | p6 | p9 | p10)) + plot_annotation(tag_levels = "a")
ggsave(allp,filename="../fig2_simulation.pdf",width=12,height=6)
