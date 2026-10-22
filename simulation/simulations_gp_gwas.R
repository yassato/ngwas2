set.seed(1234)

# load packages and functions
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
  message(paste("======iteration",i,"======"))
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
  g_nei = nei_coval(g,smap,scale=1)
  
  # m = as.matrix(dist(smap))
  # ZAdj = mapply(FUN=function(x){as.numeric(m[,x]<=1.42&m[,x]>0)},1:nrow(m))
  
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
  
  X0 = model.matrix(~poly(smap[,1],2)+poly(smap[,2],2))
  rownames(X0) = ind_names
  
  ### Parameters
  herit <- 0.6
  weights <- c(0.5, 0.5)　# changing c(0.9, 0.1) simulates DGE:IGE = 9:1
  rhosMat <- matrix(
    data = c(1, -0.6, # changing 0.6 simulates the other cases of \rho_{12}
             -0.6, 1),
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
  
  PVEnei.res = calc_PVEnei(pheno=org_y,geno=g,smap=smap,scale_seq=1)
  
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
  
  
  EM3cov.SE.res <- try(calcParamSE.cov(
    y = y,
    EM3.cov.res = EM3.res12,
    ZETA = ZETA1,
    covList = covList1,
    X0 = X0,
    REML = TRUE
  ))
  
  ysim = data.frame(org_y,y)
  res = list(PVEnei=PVEnei.res,EM3res11=EM3.res11,EM3res12=EM3.res12,EM3res12se=EM3cov.SE.res,ysim=ysim,simRes=simRes)
  saveRDS(res,file=paste0("../output/",i,"_EM3cov_neg_rc1.rds"))
  return(NULL)
}

set.seed(1234)
res = parallel::mcmapply(iter_i,1:30,mc.cores=8L)

# variance component check
# null, pos, neg for \rho_12 = 0, +0.6, -0.6
# pos03, neg03 for \rho_12 = +0.3, -0.3
# 91pos, 91neg for w1:w2=0.9:0.1 with \rho_12 = +0.6, -0.6
rho0 = c()
w01 = c()
w02 = c()
chi0 = c()
for(i in 1:30){
  out = readRDS(paste0("./output/simulation/",i,"_EM3cov_null_rc1.rds"))
  resid = out$ysim$org_y-out$EM3res12$y.pred
  print(cor(resid,out$EM3res12$u)[1,1])
  
  rho0 = c(rho0,out$EM3res12$rhosMat[1,2])
  w01 = rbind(w01,out$EM3res11$weights)
  w02 = rbind(w02,out$EM3res12$weights)
  chi0 = c(chi0,-2*(out$EM3res11$LL-out$EM3res12$LL))
}
pchisq(chi0,1,lower.tail=FALSE) # all non-significant

p01 = ggplot(data=NULL,aes(x="",y=rho0)) + geom_boxplot() + geom_jitter() + ylim(-0.5,0.5) +
  geom_hline(yintercept=0,lwd=1,lty=2) + ylab(expression(hat(rho)[12])) + xlab("cov") + theme_bw()

w0 = data.frame(noncov=w01[,2],cov=w02[,2])
p02 = w0 %>%
  gather(c("noncov","cov"),key="model",value="w2") %>%
  ggplot(aes(x=model,y=w2)) + geom_boxplot() + geom_jitter() + ylim(0,1) + 
  geom_hline(yintercept=0.5,lwd=1,lty=2) + ylab(expression(hat(w)[2])) + theme_bw()

rhon = c()
wn11 = c()
wn12 = c()
for(i in 1:30){
  out = readRDS(paste0("./output/simulation/",i,"_EM3cov_neg_rc1.rds"))
  rhon = c(rhon,out$EM3res12$rhosMat[1,2])
  wn11 = rbind(wn11,out$EM3res11$weights)
  wn12 = rbind(wn12,out$EM3res12$weights)
}
p1 = ggplot(data=NULL,aes(x="",y=rhon)) + geom_boxplot() + geom_jitter() + ylim(NA,0.0) +
  geom_hline(yintercept=-0.6,lwd=1,lty=2) + ylab(expression(hat(rho)[12])) + xlab("cov") + theme_bw()

wn = data.frame(noncov=wn11[,2],cov=wn12[,2])
p2 = wn %>%
  gather(c("noncov","cov"),key="model",value="w2") %>%
  ggplot(aes(x=model,y=w2)) + geom_boxplot() + geom_jitter() + ylim(0,1) + 
  geom_hline(yintercept=0.5,lwd=1,lty=2) + ylab(expression(hat(w)[2])) + theme_bw()

rhop = c()
wp11 = c()
wp12 = c()
for(i in 1:30){
  out = readRDS(paste0("./output/simulation/",i,"_EM3cov_pos_rc1.rds"))
  rhop = c(rhop,out$EM3res12$rhosMat[1,2])
  wp11 = rbind(wp11,out$EM3res12$weights)
  wp12 = rbind(wp12,out$EM3res12$weights)
}
p3 = ggplot(data=NULL,aes(x="",y=rhop)) + geom_boxplot() + geom_jitter() + ylim(0,NA) +
  geom_hline(yintercept=0.6,lwd=1,lty=2) + ylab(expression(hat(rho)[12])) + xlab("cov") + theme_bw()

wp = data.frame(noncov=wp11[,2],cov=wp12[,2])
p4 = wp %>%
  gather(c("noncov","cov"),key="model",value="w2") %>%
  ggplot(aes(x=model,y=w2)) + geom_boxplot() + geom_jitter() + ylim(0,1) + 
  geom_hline(yintercept=0.5,lwd=1,lty=2) + ylab(expression(hat(w)[2])) + theme_bw()


# GP output
r2s = c()
r11s = c(); r12s = c()
mae11s = c(); mae12s = c()
for(i in 1:30){
  out = readRDS(paste0("./output/simulation/",i,"_EM3cov_null_rc1.rds"))
  r11s = c(r11s,cor(out$simRes$uAll[is.na(out$ysim$y)==FALSE],out$EM3res11$u[is.na(out$ysim$y)==FALSE])^2)
  r12s = c(r12s,cor(out$simRes$uAll[is.na(out$ysim$y)==FALSE],out$EM3res12$u[is.na(out$ysim$y)==FALSE])^2)
  mae11s = c(mae11s,mean(abs(out$simRes$uAll[is.na(out$ysim$y)==FALSE]-out$EM3res11$u[is.na(out$ysim$y)==FALSE])))
  mae12s = c(mae12s,mean(abs(out$simRes$uAll[is.na(out$ysim$y)==FALSE]-out$EM3res12$u[is.na(out$ysim$y)==FALSE])))
  print(i)
}

r2_0 = data.frame(noncov=r11s,cov=r12s)
p05 = r2_0 %>%
  gather(c("noncov","cov"),key="model",value="R2") %>%
  ggplot(aes(x=model,y=R2)) + geom_boxplot() + geom_jitter() + ylab(expression(R^2)) + ylim(NA,1) + theme_bw()

mae_0 = data.frame(noncov=mae11s,cov=mae12s)
p06 = mae_0 %>%
  gather(c("noncov","cov"),key="model",value="MAE") %>%
  ggplot(aes(x=model,y=MAE)) + geom_boxplot() + geom_jitter() + ylim(0,NA) + theme_bw()

r2s = c()
r11s = c(); r12s = c()
mae11s = c(); mae12s = c()
for(i in 1:30){
  out = readRDS(paste0("./output/simulation/",i,"_EM3cov_neg_rc1.rds"))
  r11s = c(r11s,cor(out$simRes$uAll[is.na(out$ysim$y)==FALSE],out$EM3res11$u[is.na(out$ysim$y)==FALSE])^2)
  r12s = c(r12s,cor(out$simRes$uAll[is.na(out$ysim$y)==FALSE],out$EM3res12$u[is.na(out$ysim$y)==FALSE])^2)
  mae11s = c(mae11s,mean(abs(out$simRes$uAll[is.na(out$ysim$y)==FALSE]-out$EM3res11$u[is.na(out$ysim$y)==FALSE])))
  mae12s = c(mae12s,mean(abs(out$simRes$uAll[is.na(out$ysim$y)==FALSE]-out$EM3res12$u[is.na(out$ysim$y)==FALSE])))
  print(i)
}

r2_neg = data.frame(noncov=r11s,cov=r12s)
p5 = r2_neg %>%
  gather(c("noncov","cov"),key="model",value="R2") %>%
  ggplot(aes(x=model,y=R2)) + geom_boxplot() + geom_jitter() + ylab(expression(R^2)) + ylim(NA,1) + theme_bw()

mae_neg = data.frame(noncov=mae11s,cov=mae12s)
p6 = mae_neg %>%
  gather(c("noncov","cov"),key="model",value="MAE") %>%
  ggplot(aes(x=model,y=MAE)) + geom_boxplot() + geom_jitter() + ylim(0,NA) + theme_bw()


r2s = c()
r11s = c(); r12s = c()
mae11s = c(); mae12s = c()
for(i in 1:30){
  out = readRDS(paste0("./output/simulation/",i,"_EM3cov_pos_rc1.rds"))
  r11s = c(r11s,cor(out$simRes$uAll[is.na(out$ysim$y)==FALSE],out$EM3res11$u[is.na(out$ysim$y)==FALSE])^2)
  r12s = c(r12s,cor(out$simRes$uAll[is.na(out$ysim$y)==FALSE],out$EM3res12$u[is.na(out$ysim$y)==FALSE])^2)
  mae11s = c(mae11s,mean(abs(out$simRes$uAll[is.na(out$ysim$y)==FALSE]-out$EM3res11$u[is.na(out$ysim$y)==FALSE])))
  mae12s = c(mae12s,mean(abs(out$simRes$uAll[is.na(out$ysim$y)==FALSE]-out$EM3res12$u[is.na(out$ysim$y)==FALSE])))
  print(i)
}

r2_pos = data.frame(noncov=r11s,cov=r12s)
p7 = r2_pos %>%
  gather(c("noncov","cov"),key="model",value="R2") %>%
  ggplot(aes(x=model,y=R2)) + geom_boxplot() + geom_jitter() + ylab(expression(R^2)) + ylim(NA,1) + theme_bw()

mae_pos = data.frame(noncov=mae11s,cov=mae12s)
p8 = mae_pos %>%
  gather(c("noncov","cov"),key="model",value="MAE") %>%
  ggplot(aes(x=model,y=MAE)) + geom_boxplot() + geom_jitter() + ylim(0,NA) + theme_bw()

# export a supplementary figure
gpp = (p01 + labs(subtitle=expression(rho[12]*" = 0"))| p02 | p05 | p06) /  (p3 + labs(subtitle=expression(rho[12]*" = 0.3"))| p4 | p7 | p8) / (p1 + labs(subtitle=expression(rho[12]*" = −0.3"))| p2 | p5 | p6) + plot_annotation(tag_levels = "a")
ggsave(gpp,filename="../figs/figS_simulation003.pdf",width=8, height=8)


# SE check
est = c()
se = c()
for(i in 1:30){ print(i)
  out = readRDS(paste0("./output/simulation/",i,"_EM3cov_neg_rc1se.rds"))
  est = rbind(est,out$EM3res12se$estimate[2:3])
  se = rbind(se,out$EM3res12se$SE[2:3])
}

se0 = data.frame(iteration=c(1:30),est=est,se=se)
dot01 = ggplot(data=se0,aes(x=iteration,y=est.rho.A.Int)) + geom_point() + geom_errorbar(aes(ymin=est.rho.A.Int-se.rho.A.Int,ymax=est.rho.A.Int+se.rho.A.Int)) + ylim(-1.5,1.5) + ylab(expression(hat(rho)[12])) +
  geom_hline(yintercept=mean(se0$est.rho.A.Int),color="red") + geom_hline(yintercept=mean(se0$est.rho.A.Int) + sd(se0$est.rho.A.Int),lty=2,color="red") + geom_hline(yintercept=mean(se0$est.rho.A.Int) - sd(se0$est.rho.A.Int),lty=2,color="red") + theme_bw()
dot02 = ggplot(data=se0,aes(x=iteration,y=est.sigma2.Int)) + geom_point() + geom_errorbar(aes(ymin=est.sigma2.Int-se.sigma2.Int,ymax=est.sigma2.Int+se.sigma2.Int)) + ylim(0,1.0) + ylab(expression(hat(sigma)[2]^2)) +
  geom_hline(yintercept=mean(se0$est.sigma2.Int),color="red") + geom_hline(yintercept=mean(se0$est.sigma2.Int) + sd(se0$est.sigma2.Int),lty=2,color="red") + geom_hline(yintercept=mean(se0$est.sigma2.Int) - sd(se0$est.sigma2.Int),lty=2,color="red") + theme_bw()

nse = data.frame(iteration=c(1:30),est=est,se=se)
ndot1 = ggplot(data=nse,aes(x=iteration,y=est.rho.A.Int)) + geom_point() + geom_errorbar(aes(ymin=est.rho.A.Int-se.rho.A.Int,ymax=est.rho.A.Int+se.rho.A.Int)) + ylim(-1.1,0) + ylab(expression(hat(rho)[12])) +
  geom_hline(yintercept=mean(nse$est.rho.A.Int),color="red") + geom_hline(yintercept=mean(nse$est.rho.A.Int) + sd(nse$est.rho.A.Int),lty=2,color="red") + geom_hline(yintercept=mean(nse$est.rho.A.Int) - sd(nse$est.rho.A.Int),lty=2,color="red") + theme_bw()
ndot2 = ggplot(data=nse,aes(x=iteration,y=est.sigma2.Int)) + geom_point() + geom_errorbar(aes(ymin=est.sigma2.Int-se.sigma2.Int,ymax=est.sigma2.Int+se.sigma2.Int)) + ylim(0,1.0) + ylab(expression(hat(sigma)[2]^2)) +
  geom_hline(yintercept=mean(nse$est.sigma2.Int),color="red") + geom_hline(yintercept=mean(nse$est.sigma2.Int) + sd(nse$est.sigma2.Int),lty=2,color="red") + geom_hline(yintercept=mean(nse$est.sigma2.Int) - sd(nse$est.sigma2.Int),lty=2,color="red") + theme_bw()

pse = data.frame(iteration=c(1:30),est=est,se=se)
pdot1 = ggplot(data=pse,aes(x=iteration,y=est.rho.A.Int)) + geom_point() + geom_errorbar(aes(ymin=est.rho.A.Int-se.rho.A.Int,ymax=est.rho.A.Int+se.rho.A.Int)) + ylim(0,1.1) + ylab(expression(hat(rho)[12])) +
  geom_hline(yintercept=mean(pse$est.rho.A.Int),color="red") + geom_hline(yintercept=mean(pse$est.rho.A.Int) + sd(pse$est.rho.A.Int),lty=2,color="red") + geom_hline(yintercept=mean(pse$est.rho.A.Int) - sd(pse$est.rho.A.Int),lty=2,color="red") + theme_bw()
pdot2 = ggplot(data=pse,aes(x=iteration,y=est.sigma2.Int)) + geom_point() + geom_errorbar(aes(ymin=est.sigma2.Int-se.sigma2.Int,ymax=est.sigma2.Int+se.sigma2.Int)) + ylim(0,1.0) + ylab(expression(hat(sigma)[2]^2)) +
  geom_hline(yintercept=mean(pse$est.sigma2.Int),color="red") + geom_hline(yintercept=mean(pse$est.sigma2.Int) + sd(pse$est.sigma2.Int),lty=2,color="red") + geom_hline(yintercept=mean(pse$est.sigma2.Int) - sd(pse$est.sigma2.Int),lty=2,color="red") + theme_bw()

# export a supplementary figure
ppse = ((pdot1+labs(subtitle="cov+")) | (ndot1+labs(subtitle="cov-")) | (dot01+labs(subtitle="cov0"))) / (pdot2 | ndot2 | dot02) + plot_annotation(tag_levels = "a")
ggsave(ppse,filename="../figs/rho_sim_se.pdf",height=6,width=8)

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
  g_nei = nei_coval(g,smap,scale=1)
  
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
  
  X0 = model.matrix(~poly(smap[,1],2)+poly(smap[,2],2))
  rownames(X0) = ind_names
  
  ### Parameters
  herit <- 0.6
  weights <- c(0.5, 0.5)
  rhosMat <- matrix(
    data = c(1, -0.6, # changing 0.6 into -0.6 simulates the case of negative \rho_12
             -0.6, 1),
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
  saveRDS(res.roc,file=paste0("./output/simulation/",i,"_gwas_neg_rc1.rds"))
}

# GWAS output
aucs = c(); scores = c(); power = c(); gifs = c()
for(i in 1:30){
  out = readRDS(file=paste0("./output/simulation/",i,"_gwas_pos_rc1.rds")) # ncor for not correlated
  ans = as.numeric(out$beta12==0&out$beta2!=0)
  p_11 = -log10(out$res11$p2)
  p_12 = -log10(out$res12$p2)
  roc11 = roc(ans~p_11,smooth=TRUE)
  roc12 = roc(ans~p_12,smooth=TRUE)
  gif11 = median(log10(out$res11$p2))/median(log10(ppoints(length(out$res11$p2))))
  gif12 = median(log10(out$res12$p2))/median(log10(ppoints(length(out$res12$p2))))
  aucs = rbind(aucs,c(roc11$auc,roc12$auc))
  power = rbind(aucs,c(max(roc11$sensitivities[roc11$specificities>0.95]),max(roc12$sensitivities[roc12$specificities>0.95])))
  scores = rbind(scores,c(median(p_11[ans!=0]),median(p_12[ans!=0])))
  gifs = rbind(gifs,c(gif11,gif12))
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

colnames(gifs) = c("noncov","cov")
gif_neg = as.data.frame(gifs)
p13 = gif_neg %>%
  gather(c("noncov","cov"),key="model",value="GIF") %>%
  ggplot(aes(x=model,y=GIF)) + geom_boxplot() + geom_jitter() + ylab("GIF") + theme_bw()


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

colnames(gifs) = c("noncov","cov")
gif_pos = as.data.frame(gifs)
p14 = gif_pos %>%
  gather(c("noncov","cov"),key="model",value="GIF") %>%
  ggplot(aes(x=model,y=GIF)) + geom_boxplot() + geom_jitter() + ylab("GIF") + theme_bw()


# export a main figure
allp=((p3+labs(title="cov+") | p4 | p7 | p8 | p11 | p12 | p14) / (p1+labs(title="cov-") | p2 | p5 | p6 | p9 | p10 | p13)) + plot_annotation(tag_levels = "a")
ggsave(allp,filename="../figs/fig2_simulation.pdf",width=12,height=5)

# export a supplementary figure
pncor=(p11+labs(title="cov+") | p12 | p14) /  (p9+labs(title="cov-") | p10 | p13) + plot_annotation(tag_levels = "a")
ggsave(pncor,filename="../figs/figS_sim_gwas_ncor.pdf",width=6,height=5)
