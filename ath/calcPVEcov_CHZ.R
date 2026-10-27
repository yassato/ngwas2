########################
#PVEcov for Zurich data#
########################

library(rNeighborGWAS)
library(gaston)
library(RAINBOWR)

# load data
# input data are available at Zenodo at https://doi.org/10.5281/zenodo.7945317 
geno_d = readRDS(file="./genoData/sub_snpMAF5.rds")
geno_d[geno_d==0]=-1 #replace 0 into -1

position = readRDS(file="./genoData/positionsMAF5.rds")

pheno_d = read.csv("./insectData/SurveyCombined4GWAS_max.csv",header=TRUE)
pheno_d = subset(pheno_d, Site=="ZH")
pheno_d = pheno_d[-which(is.na(pheno_d$gwasID)),]
pheno_d$gwasID = paste0("X",pheno_d$gwasID)

# reshape pheno.data
pheno_d = subset(pheno_d, gwasID!="X7329")
pheno_d = subset(pheno_d, gwasID!="X3")
naID = which(is.na(pheno_d$InitLeafLen))
pheno_d = pheno_d[-naID,]

n_marker = nrow(geno_d)
n_plants = nrow(pheno_d)

geno = geno_d[,as.character(pheno_d$gwasID)]
geno = t(geno)

# set spatial map
smap = cbind(pheno_d$position_X,pheno_d$position_Y)

# spatial scale sqrt(2)+0.01 for the first nearest neighbors s
scale = sqrt(2)+0.01

g_nei = nei_coval(geno,smap=smap,scale=scale,grouping=pheno_d$Block,n_core=1L)

q = ncol(geno)
K_self = tcrossprod(geno)
K_self = ((q-1)/2+K_self/2)/(q-1)
# K_self = as.matrix(Matrix::nearPD(K_self,maxit=10^6)$mat)

K_nei = tcrossprod(g_nei)/(q-1)
# K_nei = as.matrix(Matrix::nearPD(K_nei,maxit=10^6)$mat)

ind_names = paste0("ind",1:nrow(geno))
Z = diag(nrow(geno))
rownames(Z) = ind_names
colnames(Z) = ind_names
rownames(K_self) = ind_names
colnames(K_self) = ind_names
rownames(K_nei) = ind_names
colnames(K_nei) = ind_names
ZETA1 <- list(A = list(Z = Z, K = K_self),
              Int = list(Z = Z, K = K_nei))

K12 = tcrossprod(geno,g_nei) / (q-1)
K21 = t(K12)
rownames(K12) = ind_names
colnames(K12) = ind_names
rownames(K21) = ind_names
colnames(K21) = ind_names

covList1 <- rep(list(rep(list(NULL), 2)), 2)
covList1[[1]][[2]] <- K12
covList1[[2]][[1]] <- K21

X = as.matrix(model.matrix(~factor(pheno_d$Block)+scale(pheno_d$InitLeafLen)+pheno_d$Bolting+pheno_d$edge))
rownames(X) = ind_names

Y = as.matrix(scale(log(pheno_d$Holes+1)))
EM3.res11 <- EM3.cpp(y = Y,
                     X0 = X,
                     ZETA = ZETA1)
w = EM3.res11$weights
if(sum(w==TRUE)==1) {
  w[w==1] = 0.99
  w[w==0] = 0.01
}

EM3.res12 <- EM3.cov(y = Y,
                     X0 = X,
                     ZETA = ZETA1,
                     covList = covList1,
                     optimizeWeights = TRUE,
                     parInitForWeights = w,
                     forceApproxK = TRUE,
                     nIterOptimization=10^5, traceInside=1)

EM3.res12se = calcParamSE.cov(
  y = Y,
  EM3.cov.res = EM3.res12,
  ZETA = ZETA1,
  covList = covList1,
  X0 = X,
  REML = TRUE,
  check.Vinv = TRUE,
  Vinv.tol = 1e-06,
  suppressWarning = TRUE
)

saveRDS(list(EM3res11=EM3.res11,EM3res12=EM3.res12,EM3res12se=EM3.res12se),"../EM3outHolesZH.rds")


Y = as.matrix(scale(log(pheno_d$richness+1)))
EM3.res11 <- EM3.cpp(y = Y,
                     X0 = X,
                     ZETA = ZETA1)
w = EM3.res11$weights
if(sum(w==TRUE)==1) {
  w[w==1] = 0.99
  w[w==0] = 0.01
}
EM3.res12 <- EM3.cov(y = Y,
                     X0 = X,
                     ZETA = ZETA1,
                     covList = covList1,
                     optimizeWeights = TRUE,
                     parInitForWeights = w,
                     forceApproxK = TRUE,
                     nIterOptimization=10^5, traceInside=1)

saveRDS(list(EM3res11=EM3.res11,EM3res12=EM3.res12),"../EM3outRichness.rds")

Y = as.matrix(scale(log(pheno_d$sucker+1)))
EM3.res11 <- EM3.cpp(y = Y,
                     X0 = X,
                     ZETA = ZETA1)
w = EM3.res11$weights
if(sum(w==TRUE)==1) {
  w[w==1] = 0.99
  w[w==0] = 0.01
}
EM3.res12 <- EM3.cov(y = Y,
                     X0 = X,
                     ZETA = ZETA1,
                     covList = covList1,
                     optimizeWeights = TRUE,
                     parInitForWeights = w,
                     forceApproxK = TRUE,
                     nIterOptimization=10^5, traceInside=1)
saveRDS(list(EM3res11=EM3.res11,EM3res12=EM3.res12),"../EM3outSucker.rds")

Y = as.matrix(scale(log(pheno_d$chewer+1)))
EM3.res11 <- EM3.cpp(y = Y,
                     X0 = X,
                     ZETA = ZETA1)
w = EM3.res11$weights
if(sum(w==TRUE)==1) {
  w[w==1] = 0.99
  w[w==0] = 0.01
}
EM3.res12 <- EM3.cov(y = Y,
                     X0 = X,
                     ZETA = ZETA1,
                     covList = covList1,
                     optimizeWeights = TRUE,
                     parInitForWeights = w,
                     forceApproxK = TRUE,
                     nIterOptimization=10^5, traceInside=1)
saveRDS(list(EM3res11=EM3.res11,EM3res12=EM3.res12),"../EM3outChewer.rds")

