library(tidyverse)
library(patchwork)

##############################################
# Figure S1: spatial simulations based on the Ising model

# define functions
XiXj1 = function(id,dmat,range) {
  xi = Xi[id]
  if(xi==0) { xi = 1 } else { xi = xi }
  j = as.numeric(which((dmat[,id]<=range)&(dmat[,id]>0)))
  xj = Xi[j]
  xj[xj==0] = 1
  return(xi*sum(xj))
}

# set coefficients
J = -0.2 #= beta_2; -0.2 and 0.2 for Fig. S1a and b, respectively
h = 0.0 #= beta_1

# set a field
rect = 40
N = rect*rect
Xi = sample(c(-1,1),N,replace=T)

smap = cbind(rep(1:rect,each=rect),
             rep(1:rect,rect)
)
dmat = as.matrix(dist(smap))

group = rep(1,rect*rect)


xixj = c()
for(i in 1:N) xixj = c(xixj, XiXj1(i,dmat=dmat,range=sqrt(2)))
Ei_t0 = (J*xixj+h*Xi)

# MCMC by Gibbs sampling
for(j in 1:100) {
  Xi_t0 = Xi
  
  perm = sample(1:N,N)
  for(i in perm) {
    Xi[i] = sample(c(-1,1),1)
    xi = Xi[i]
    if(xi==0) { xi = 1 }
    Ei = (J*XiXj1(i,dmat=dmat,range=sqrt(2))+h*xi) # set the energy positive to maximize itself
    
    # Metropolis algorithm
    if(Ei_t0[i]<=Ei) {
      Xi[i] = Xi[i]; Ei_t0[i] = Ei
    } else if(runif(1,0,1)<exp(Ei-Ei_t0[i])) {
      Xi[i] = Xi[i]; Ei_t0[i] = Ei
    } else {
      Xi[i] = Xi_t0[i]
    }
  }
  E = mean(Ei_t0)
  image(matrix(Xi,rect,rect),main=paste(j,E),col=c("black","white"))
}


pdf(file=paste0("../IsingJ",J,"_h",h,"cnt_d1.pdf"), width=5, height=5)
image(t(matrix(Xi,rect,rect)),main=paste(j,E),col=c("black","white"),las=1,xaxt="n",yaxt="n")
dev.off()

############
# Figure S2
W_AA = function(x,b1,b2,b12) (b12+b2)*(2*x-1)+b1+1
W_Aa = function(x,b1,b2,b12) 1
W_aa = function(x,b1,b2,b12) (b12-b2)*(2*x-1)-b1+1

W_A = function(x,b1,b2,b12) x*W_AA(x,b1,b2,b12)+(1-x)*W_Aa(x,b1,b2,b12)
W_a = function(x,b1,b2,b12) x*W_Aa(x,b1,b2,b12)+(1-x)*W_aa(x,b1,b2,b12)
W_m = function(x,b1,b2,b12) x*W_A(x,b1,b2,b12)+(1-x)*W_a(x,b1,b2,b12)

wplot = function(b1,b2,b12) {
  p2 = 0.5*(b12-b2)/b12
  w2 = p2*W_A(p2,b1,b2,b12)+(1-p2)*W_a(p2,b1,b2,b12)
  p = ggplot(NULL,aes(x=0.5,y=0.5*W_A(0.5,b1,b2,b12)+(1-0.5)*W_a(0.5,b1,b2,b12))) + 
    ylim(0.5,1.5) + xlim(0,1) + theme_classic() +
    geom_function(fun=W_A, col=grey(0.0,1.0),args=list(b1,b2,b12)) +
    geom_function(fun=W_a, col=grey(0.0,0.5),args=list(b1,b2,b12)) + 
    geom_function(fun=W_m, lty=2,args=list(b1,b2,b12)) + 
    ylab(expression(y[k[i]])) + xlab(expression(f["q,A"])) +
    labs(subtitle=bquote(beta["q,s"]*"="*.(b1)*"; "*beta["q,n"]*"="*.(b12)*"; "*beta["q,sxn"]*"="*.(b2)))
  return(p)
} 

pa = wplot(b1=0.1,b2=-0.2,b12=0.0)
pa = pa + geom_text(aes(x=0.025,y=1.5),label="(black): A alleles",size=3,hjust=0) + 
  geom_text(aes(x=0.025,y=1.4),label="(gray): a alleles",colour=grey(0.5,1),size=3,hjust=0) +
  geom_text(aes(x=0.025,y=1.3),label="(dashed): weighted mean",size=3,hjust=0)
pb = wplot(b1=0.0,b2=-0.2,b12=0.3)
pc = wplot(b1=0.1,b2=-0.2,b12=-0.3)
pd = wplot(b1=0.1,b2=0.2,b12=0)
pe = wplot(b1=0.0,b2=0.2,b12=0.3)
pf = wplot(b1=0.1,b2=0.2,b12=-0.3)

p = (pa | pb | pc) / (pd | pe | pf) + plot_annotation(tag_levels = "a")
ggsave(p,filename="../figs_numerical.pdf",width=9,height=6)

