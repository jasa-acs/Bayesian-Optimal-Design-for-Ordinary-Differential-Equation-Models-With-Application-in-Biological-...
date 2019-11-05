##################################################################################################

### PLACENTA #####################################################################################

##################################################################################################

PlacentaSolver<-function(grid, x2u2, B, model.uncertainty = FALSE){

X<-x2u2	
	
aa<-c(80,0.02,80,80)
bb<-c(120,0.08,120,120)

mods<-rep(2,B)
if(model.uncertainty){
mods<-rbinom(n=B,size=1,prob=0.5)+1}

THETA2<-cbind(rtriangle(n=B,a=aa[1],b=bb[1]),rtriangle(n=B,a=aa[2],b=bb[2]),rtriangle(n=B,a=aa[3],b=bb[3]),rtriangle(n=B,a=aa[4],b=bb[4]))
phi<-matrix(runif(n=4*B,min=0,max=0.05),nrow=B)

THETA<-0*THETA2
for(k in 1:4){
u<-runif(B,min=1-phi[,k],max=1+phi[,k])
THETA[,k]<-u*THETA2[,k]}

THETA2[mods==1,4]<-THETA2[mods==1,3]
THETA[mods==1,4]<-THETA[mods==1,3]

n0<-length(grid)

qrbeval<-placentaqr(grid,placenta$S[1:placenta$N0])%*%(placenta$Binv/placenta$alpha)
qqqq<-placentaqq(grid,grid)/placenta$alpha
rqrq<-placentarq(placenta$S[1:placenta$N0],grid)/placenta$alpha
sCCC<-diag(qqqq-qrbeval%*%rqrq)
sCCC<-ifelse(is.na(sCCC),0,sCCC)
sCCC<-diag(sqrt(ifelse(sCCC<0,0,sCCC)))

mu<-placsolver1.cpp(X=matrix(X,nrow=1),S=placenta$S,THETA=THETA,qrbeval=qrbeval,QRB=placenta$QRB,SAPV=placenta$SAPV,NM=c(n0,1), randoms1 = matrix(rnorm(B*length(placenta$SAPV)),ncol=B),randoms2 = matrix(rnorm(B*length(placenta$SAPV)),ncol=B))
mu<-mu+matrix(as.vector(sCCC%*%matrix(rnorm(n=n0*B),nrow=n0)),byrow=TRUE,nrow=B)

list(mu=mu,mods=mods)}

PlacentaDefaultEvals<-function(M, utility = "NSEL"){

if(M!=7){
stop("Only implemented for M = 7")  
}
    
if(utility=="EST01"){
out<-placentaopt$dist.def}
if(utility=="NSEL"){
out<-placentaopt$nsel.def}
if(utility=="NAEL"){
out<-placentaopt$nael.def}
if(utility=="MOD01"){
out<-placentaopt$mod.def}
	
out}

PlacentaEvals<-function(M, rep = NULL, utility = "NSEL"){

if(!is.null(rep)){	
if(rep<0|rep>20){
stop("rep should be 1, 2, ..., 20")}} 
	
if(utility=="EST01"){
out<-placentaopt$dist.out[M-1,,]}
if(utility=="NSEL"){
out<-placentaopt$nsel.out[M-1,,]}
if(utility=="NAEL"){
out<-placentaopt$nael.out[M-1,,]}
if(utility=="MOD01"){
out<-placentaopt$mod.out[M-1,,]}
	
if(!is.null(rep)){
out<-out[rep,]}
	
out}
	
PlacentaTerminal<-function(M, rep, utility = "NSEL", scaled = TRUE){

if(utility=="EST01"){
out<-placentaopt$dist.designs[placentaopt$dist.designs[,1]==M & placentaopt$dist.designs[,2]==rep,-c(1,2)]}
if(utility=="NSEL"){
out<-placentaopt$nsel.designs[placentaopt$nsel.designs[,1]==M & placentaopt$nsel.designs[,2]==rep,-c(1,2)]}
if(utility=="NAEL"){
out<-placentaopt$nael.designs[placentaopt$nael.designs[,1]==M & placentaopt$nael.designs[,2]==rep,-c(1,2)]}
if(utility=="MOD01"){
out<-placentaopt$mod.designs[placentaopt$mod.designs[,1]==M & placentaopt$mod.designs[,2]==rep,-c(1,2)]}

if(!scaled){
out<-list(concentrations = matrix((out[1:M,]+1)*500,ncol=2),times = as.vector(out[-(1:M),]+1)*300)
conc.ord <- do.call(order, data.frame(out$concentrations))
time.ord <- do.call(order, data.frame(out$times))
out$concentrations <- out$concentrations[conc.ord, ]
out$times <- out$times[time.ord]
} 

out}	

PlacentaFinal<-function(M, utility = "NSEL", scaled = TRUE){

if(utility=="EST01"){
out<-placentaopt$dist.design
out<-out[out[,1]==M,-1]}
if(utility=="NSEL"){
out<-placentaopt$nsel.design
out<-out[out[,1]==M,-1]}
if(utility=="NAEL"){
out<-placentaopt$nael.design
out<-out[out[,1]==M,-1]}
if(utility=="MOD01"){
out<-placentaopt$mod.design
out<-out[out[,1]==M,-1]}

if(!scaled){
out<-list(concentrations = matrix((out[1:M,]+1)*500,ncol=2),times = as.vector(out[-(1:M),]+1)*300)
conc.ord <- do.call(order, data.frame(out$concentrations))
time.ord <- do.call(order, data.frame(out$times))
out$concentrations <- out$concentrations[conc.ord, ]
out$times <- out$times[time.ord]
}
out}

placentalimits2<-function(d){
  x<-c(-1,sort(runif(9998))*2-1,1)
  n<-length(d)
  for(t in (1:n)){
    int<-d[t]+c(-1,1)*(5/297.5)
    x<-x[x<int[1] | x>int[2]]}
  x}

PlacentaStart<-function(M){
start.d<-runif(1)*2-1
for(j in 2:8){
  cand<-placentalimits2(d=start.d)
  start.d<-c(start.d,cand[sample(x=1:length(cand),size=1)])
  start.d<-sort(as.vector(start.d))}
start.d<-rbind(randomLHS(n=M,k=2)*2-1,matrix(start.d,ncol=2))
start.d}

placsolver1.cpp <- function(X, S, THETA, qrbeval, QRB, SAPV, NM, randoms1, randoms2){
	.Call( "placsolver1cpp", X, S, THETA, qrbeval, QRB, SAPV, NM, randoms1, randoms2, PACKAGE = "aceodes" )
}

nsel4.cpp <- function(THETA, sig, mu, y){
	.Call( "nsel4cpp", THETA, sig, mu, y, PACKAGE = "aceodes" )
}

# sig.cpp <- function(sig2, mu1, mu2, y){
# 	.Call( "sigcpp", sig2, mu1, mu2, y, PACKAGE = "aceodes" )
# }

nael.cpp <- function(THETA, sig, mu, y, order){
	.Call( "naelcpp", THETA, sig, mu, y, order, PACKAGE = "aceodes" )
}

mod.cpp <- function(sig, mu1, mu2, y){
	.Call( "modcpp", sig, mu1, mu2, y, PACKAGE = "aceodes" )
}

erf<-function(x){
2*pnorm(sqrt(2)*x) - 1}

placentarr<-function(u1,v1){
N1<-length(u1)
N2<-length(v1)
u<-matrix(rep(u1,N2),ncol=N2,nrow=N1,byrow=FALSE)
v<-matrix(rep(v1,N1),ncol=N2,nrow=N1,byrow=TRUE)
exp(-0.25*((u-v)^2)/placenta$lam2)*spi*placenta$lam}

placentarq<-function(u1,v1){
N1<-length(u1)
N2<-length(v1)
u<-matrix(rep(u1,N2),ncol=N2,nrow=N1,byrow=FALSE)
v<-matrix(rep(v1,N1),ncol=N2,nrow=N1,byrow=TRUE)
pi*placenta$lam2*(erf(0.5*(v-u)/placenta$lam)+erf(0.5*(u-placenta$a)/placenta$lam))}

placentaqr<-function(u1,v1){
t(placentarq(v1,u1))}

placentaqq<-function(u1,v1){
N1<-length(u1)
N2<-length(v1)
u<-matrix(rep(u1,N2),ncol=N2,nrow=N1,byrow=FALSE)
v<-matrix(rep(v1,N1),ncol=N2,nrow=N1,byrow=TRUE)
pi*placenta$lam2*(u-placenta$a)*erf(0.5*(u-placenta$a)/placenta$lam)+2*spi*placenta$lam3*exp(-0.25*((u-placenta$a)^2)/placenta$lam2)-pi*placenta$lam2*(v-u)*erf(0.5*(v-u)/placenta$lam)-2*spi*placenta$lam3*exp(-0.25*((v-u)^2)/placenta$lam2)+pi*placenta$lam2*(v-placenta$a)*erf(0.5*(v-placenta$a)/placenta$lam)+2*spi*placenta$lam3*exp(-0.25*((v-placenta$a)^2)/placenta$lam2)-2*spi*placenta$lam3}

PlacentaMod01<-function(d, B){

aa<-c(80,0.02,80,80)
bb<-c(120,0.08,120,120)
P<-1
N<-8
M<-dim(d)[1]-0.5*N
n<-N*P*M
  
BB<-B+B
BBB<-BB+B

mods<-rbinom(n=B,size=1,prob=0.5)+1

THETA2<-cbind(rtriangle(n=B,a=aa[1],b=bb[1]),rtriangle(n=B,a=aa[2],b=bb[2]),rtriangle(n=B,a=aa[3],b=bb[3]),rtriangle(n=B,a=aa[4],b=bb[4]))
phi<-matrix(runif(n=4*B,min=0,max=0.05),nrow=B)

THETA<-matrix(0,nrow=B,ncol=4*M)
for(j in 1:M){
for(k in 1:4){
u<-runif(B,min=1-phi[,k],max=1+phi[,k])
THETA[,(j-1)*4+k]<-u*THETA2[,k]}}

THETA2[mods==1,4]<-THETA2[mods==1,3]
THETA[mods==1,4*(1:M)]<-THETA[mods==1,4*(1:M)-1]

theta2<-cbind(rtriangle(n=B,a=aa[1],b=bb[1]),rtriangle(n=B,a=aa[2],b=bb[2]),rtriangle(n=B,a=aa[3],b=bb[3]),rtriangle(n=B,a=aa[4],b=bb[4]))
phi<-matrix(runif(n=4*B,min=0,max=0.05),nrow=B)

theta<-matrix(0,nrow=B,ncol=4*M)
for(j in 1:M){
for(k in 1:4){
u<-runif(B,min=1-phi[,k],max=1+phi[,k])
theta[,(j-1)*4+k]<-u*theta2[,k]}}
theta2[,4]<-theta2[,3]
theta[,4*(1:M)]<-theta[,4*(1:M)-1]

THETA<-rbind(THETA,theta)
THETA2<-rbind(THETA2,theta2)

theta2<-cbind(rtriangle(n=B,a=aa[1],b=bb[1]),rtriangle(n=B,a=aa[2],b=bb[2]),rtriangle(n=B,a=aa[3],b=bb[3]),rtriangle(n=B,a=aa[4],b=bb[4]))
phi<-matrix(runif(n=4*B,min=0,max=0.05),nrow=B)

theta<-matrix(0,nrow=B,ncol=4*M)
for(j in 1:M){
for(k in 1:4){
u<-runif(B,min=1-phi[,k],max=1+phi[,k])
theta[,(j-1)*4+k]<-u*theta2[,k]}}

THETA<-rbind(THETA,theta)
THETA2<-rbind(THETA2,theta2)

sig<-runif(n=BB,min=0,max=1)

d0<-d[1:M,]
X<-1000*0.5*(matrix(d0,nrow=M)+1)
grid<-sort((as.vector(d[-(1:M),])+1)*0.5*595+5)

qrbeval<-placentaqr(grid,placenta$S[1:placenta$N0])%*%(placenta$Binv/placenta$alpha)
qqqq<-placentaqq(grid,grid)/placenta$alpha
rqrq<-placentarq(placenta$S[1:placenta$N0],grid)/placenta$alpha
sCCC<-diag(qqqq-qrbeval%*%rqrq)
sCCC<-ifelse(is.na(sCCC),0,sCCC)
sCCC<-diag(sqrt(ifelse(sCCC<0,0,sCCC)))

mu<-placsolver1.cpp(X=X,S=placenta$S,THETA=THETA,qrbeval=qrbeval,QRB=placenta$QRB,SAPV=placenta$SAPV,NM=c(N,M), randoms1 = matrix(rnorm(BBB*length(placenta$SAPV)*M),ncol=BBB),randoms2 = matrix(rnorm(BBB*length(placenta$SAPV)*M),ncol=BBB))
mu<-mu+matrix(as.vector(sCCC%*%matrix(rnorm(n=M*N*BBB),nrow=N)),byrow=TRUE,nrow=BBB)

y<-mu[1:B,]+matrix(rep(sqrt(sig[1:B]),n),byrow=FALSE,ncol=n)*matrix(rnorm(n*B),ncol=n)

mu<-mu[-(1:B),]
sig<-sig[-(1:B)]

out<-mod.cpp(sig=sig,mu1=mu[1:B,],mu2=mu[-(1:B),],y=y) 

out2<-as.numeric(as.vector(out)==mods)
out2[is.na(out2)]<-0  # Set any NAs to 0

out2
}

PlacentaNael<-function(d, B){

aa<-c(80,0.02,80,80)
bb<-c(120,0.08,120,120)
P<-1
N<-8
M<-dim(d)[1]-0.5*N
n<-N*P*M  
  
BB<-B+B

THETA2<-cbind(rtriangle(n=BB,a=aa[1],b=bb[1]),rtriangle(n=BB,a=aa[2],b=bb[2]),rtriangle(n=BB,a=aa[3],b=bb[3]),rtriangle(n=BB,a=aa[4],b=bb[4]))
phi<-matrix(runif(n=4*BB,min=0,max=0.05),nrow=BB)

##
THETA<-matrix(0,nrow=BB,ncol=4*M)
for(j in 1:M){
for(k in 1:4){
u<-runif(BB,min=1-phi[,k],max=1+phi[,k])
THETA[,(j-1)*4+k]<-u*THETA2[,k]}}

sig<-runif(n=BB,min=0,max=1)

d0<-d[1:M,]
X<-1000*0.5*(matrix(d0,nrow=M)+1)
grid<-sort((as.vector(d[-(1:M),])+1)*0.5*595+5)

qrbeval<-placentaqr(grid,placenta$S[1:placenta$N0])%*%(placenta$Binv/placenta$alpha)
qqqq<-placentaqq(grid,grid)/placenta$alpha
rqrq<-placentarq(placenta$S[1:placenta$N0],grid)/placenta$alpha
sCCC<-diag(qqqq-qrbeval%*%rqrq)
sCCC<-ifelse(is.na(sCCC),0,sCCC)
sCCC<-diag(sqrt(ifelse(sCCC<0,0,sCCC)))

mu<-placsolver1.cpp(X=X,S=placenta$S,THETA=THETA,qrbeval=qrbeval,QRB=placenta$QRB,SAPV=placenta$SAPV,NM=c(N,M), randoms1 = matrix(rnorm(BB*length(placenta$SAPV)*M),ncol=BB),randoms2 = matrix(rnorm(BB*length(placenta$SAPV)*M),ncol=BB))
mu<-mu+matrix(as.vector(sCCC%*%matrix(rnorm(n=M*N*BB),nrow=N)),byrow=TRUE,nrow=BB)

y<-mu[1:B,]+matrix(rep(sqrt(sig[1:B]),n),byrow=FALSE,ncol=n)*matrix(rnorm(n*B),ncol=n)

mu<-matrix(mu[-(1:B),],nrow=B)
sig<-sig[-(1:B)]

out<-nael.cpp(THETA=THETA2[-(1:B),],sig=sig,mu=mu,y=y,order=apply(THETA2[-(1:B),],2,order)) 

out2<-apply(abs(THETA2[1:B,]-out),1,sum)
out2[is.na(out2)]<-12.5  # Set any NAs to approximate performance of original design

-out2}

PlacentaNsel<-function(d, B){

aa<-c(80,0.02,80,80)
bb<-c(120,0.08,120,120)
P<-1
N<-8
M<-dim(d)[1]-0.5*N
n<-N*P*M  
  
BB<-B+B

THETA2<-cbind(rtriangle(n=BB,a=aa[1],b=bb[1]),rtriangle(n=BB,a=aa[2],b=bb[2]),rtriangle(n=BB,a=aa[3],b=bb[3]),rtriangle(n=BB,a=aa[4],b=bb[4]))
phi<-matrix(runif(n=4*BB,min=0,max=0.05),nrow=BB)

##
THETA<-matrix(0,nrow=BB,ncol=4*M)
for(j in 1:M){
for(k in 1:4){
u<-runif(BB,min=1-phi[,k],max=1+phi[,k])
THETA[,(j-1)*4+k]<-u*THETA2[,k]}}

sig<-runif(n=BB,min=0,max=1)

d0<-d[1:M,]
X<-1000*0.5*(matrix(d0,nrow=M)+1)
grid<-sort((as.vector(d[-(1:M),])+1)*0.5*595+5)

qrbeval<-placentaqr(grid,placenta$S[1:placenta$N0])%*%(placenta$Binv/placenta$alpha)
qqqq<-placentaqq(grid,grid)/placenta$alpha
rqrq<-placentarq(placenta$S[1:placenta$N0],grid)/placenta$alpha
sCCC<-diag(qqqq-qrbeval%*%rqrq)
sCCC<-ifelse(is.na(sCCC),0,sCCC)
sCCC<-diag(sqrt(ifelse(sCCC<0,0,sCCC)))

mu<-placsolver1.cpp(X=X,S=placenta$S,THETA=THETA,qrbeval=qrbeval,QRB=placenta$QRB,SAPV=placenta$SAPV,NM=c(N,M), randoms1 = matrix(rnorm(BB*length(placenta$SAPV)*M),ncol=BB),randoms2 = matrix(rnorm(BB*length(placenta$SAPV)*M),ncol=BB))
mu<-mu+matrix(as.vector(sCCC%*%matrix(rnorm(n=M*N*BB),nrow=N)),byrow=TRUE,nrow=BB)

y<-mu[1:B,]+matrix(rep(sqrt(sig[1:B]),n),byrow=FALSE,ncol=n)*matrix(rnorm(n*B),ncol=n)

mu<-matrix(mu[-(1:B),],nrow=B)
sig<-sig[-(1:B)]

out<-nsel4.cpp(THETA=THETA2[-(1:B),],sig=sig,mu=mu,y=y) 

out2<-apply((THETA2[1:B,]-out)^2,1,sum)
out2[is.na(out2)]<-80  # Set any NAs to approximate performance of original design

-out2}

PlacentaEst01<-function(d, B){

aa<-c(80,0.02,80,80)
bb<-c(120,0.08,120,120)
P<-1
N<-8
M<-dim(d)[1]-0.5*N
n<-N*P*M  
  
BB<-B+B

THETA2<-cbind(rtriangle(n=BB,a=aa[1],b=bb[1]),rtriangle(n=BB,a=aa[2],b=bb[2]),rtriangle(n=BB,a=aa[3],b=bb[3]),rtriangle(n=BB,a=aa[4],b=bb[4]))
phi<-matrix(runif(n=4*BB,min=0,max=0.05),nrow=BB)

##
THETA<-matrix(0,nrow=BB,ncol=4*M)
for(j in 1:M){
for(k in 1:4){
u<-runif(BB,min=1-phi[,k],max=1+phi[,k])
THETA[,(j-1)*4+k]<-u*THETA2[,k]}}

sig<-runif(n=BB,min=0,max=1)

d0<-d[1:M,]
X<-1000*0.5*(matrix(d0,nrow=M)+1)
grid<-sort((as.vector(d[-(1:M),])+1)*0.5*595+5)

qrbeval<-placentaqr(grid,placenta$S[1:placenta$N0])%*%(placenta$Binv/placenta$alpha)
qqqq<-placentaqq(grid,grid)/placenta$alpha
rqrq<-placentarq(placenta$S[1:placenta$N0],grid)/placenta$alpha
sCCC<-diag(qqqq-qrbeval%*%rqrq)
sCCC<-ifelse(is.na(sCCC),0,sCCC)
sCCC<-diag(sqrt(ifelse(sCCC<0,0,sCCC)))

mu<-placsolver1.cpp(X=X,S=placenta$S,THETA=THETA,qrbeval=qrbeval,QRB=placenta$QRB,SAPV=placenta$SAPV,NM=c(N,M), randoms1 = matrix(rnorm(BB*length(placenta$SAPV)*M),ncol=BB),randoms2 = matrix(rnorm(BB*length(placenta$SAPV)*M),ncol=BB))
mu<-mu+matrix(as.vector(sCCC%*%matrix(rnorm(n=M*N*BB),nrow=N)),byrow=TRUE,nrow=BB)

y<-mu[1:B,]+matrix(rep(sqrt(sig[1:B]),n),byrow=FALSE,ncol=n)*matrix(rnorm(n*B),ncol=n)

mu<-matrix(mu[-(1:B),],nrow=B)
sig<-sig[-(1:B)]

out<-nsel4.cpp(THETA=THETA2[-(1:B),],sig=sig,mu=mu,y=y) 

delta<-c(5,0.01,5,5)
Delta<-matrix(rep(delta,each=B),nrow=B)

out2<-(out>(THETA2[1:B,]-Delta))&(out<(THETA2[1:B,]+Delta))

out3<-as.numeric(apply(out2,1,all))
out3[is.na(out3)]<-0  # Set any NAs to 0

out3}

PlacentaLimits<-function(i,j,d){

N<-8
M<-dim(d)[1]-0.5*N

hered<-d
if(i>M){
x<-c(-1,sort(runif(9998))*2-1,1)
exc<-as.vector(hered[-(1:M),])
n<-length(exc)
for(t in (1:n)[exc!=hered[i,j]]){
int<-exc[t]+c(-1,1)*(5/297.5)
x<-x[x<=int[1] | x>=int[2]]}} else{
x<-c(-1,sort(runif(9998))*2-1,1)}
x}

##################################################################################################

### JAKSTAT #####################################################################################

##################################################################################################

jakstatlimits2<-function(d){
  x<-c(-1,sort(runif(9998))*2-1,1)
  n<-dim(d)[1]
  for(t in (1:n)){
    int<-d[t,]+c(-1,1)*(1/30)
    x<-x[x<int[1] | x>int[2]]}
  x}

JakstatStart<-function(){
start.d<-matrix(runif(1)*2-1,ncol=1)
for(j in 2:16){
  cand<-jakstatlimits2(d=start.d)
  start.d<-rbind(start.d,cand[sample(x=1:length(cand),size=1)])
  start.d<-matrix(sort(as.vector(start.d)),ncol=1)}
rbind(runif(1)*2-1,start.d)}

jakstatrr<-function(u1,v1){
N1<-length(u1)
N2<-length(v1)
u<-matrix(rep(u1,N2),ncol=N2,nrow=N1,byrow=FALSE)
v<-matrix(rep(v1,N1),ncol=N2,nrow=N1,byrow=TRUE)
ret<-ifelse(u>=v,v,u)-ifelse(u<=v,v,u)
(ret+2*jakstat$lam)*ifelse(ret>=(-2*jakstat$lam),1,0)}

jakstatrq<-function(u1,v1){
N1<-length(u1)
N2<-length(v1)
u<-matrix(rep(u1,N2),ncol=N2,nrow=N1,byrow=FALSE)
v<-matrix(rep(v1,N1),ncol=N2,nrow=N1,byrow=TRUE)
A<-matrix(jakstat$a,nrow=N1,ncol=N2)
2*jakstat$lam*(pmin(v-jakstat$lam,u+jakstat$lam)-pmax(A+jakstat$lam,u-jakstat$lam))*as.numeric(pmin(v-jakstat$lam,u+jakstat$lam)>pmax(A+jakstat$lam,u-jakstat$lam))+(0.5*pmin(A+jakstat$lam,pmin(v-jakstat$lam,u+jakstat$lam))^2+(jakstat$lam-jakstat$a)*pmin(A+jakstat$lam,pmin(v-jakstat$lam,u+jakstat$lam))-0.5*pmax(A-jakstat$lam,u-jakstat$lam)^2 -(jakstat$lam-jakstat$a)*pmax(A-jakstat$lam,u-jakstat$lam))*as.numeric(pmin(A+jakstat$lam,pmin(v-jakstat$lam,u+jakstat$lam))>pmax(A-jakstat$lam,u-jakstat$lam))+((v+jakstat$lam)*pmin(v+jakstat$lam,u+jakstat$lam) - 0.5*pmin(v+jakstat$lam,u+jakstat$lam)^2-(v+jakstat$lam)*pmax(A+jakstat$lam,pmax(v-jakstat$lam,u-jakstat$lam))+ 0.5*pmax(A+jakstat$lam,pmax(v-jakstat$lam,u-jakstat$lam))^2)*as.numeric(pmin(v+jakstat$lam,u+jakstat$lam)>pmax(A+jakstat$lam,pmax(v-jakstat$lam,u-jakstat$lam)))+ ((v-jakstat$a)*pmin(A+jakstat$lam,u+jakstat$lam) - (v-jakstat$a)*pmax(v-jakstat$lam,u-jakstat$lam))*as.numeric(pmin(A+jakstat$lam,u+jakstat$lam)>pmax(v-jakstat$lam,u-jakstat$lam))
}

jakstatqr<-function(u1,v1){
t(jakstatrq(v1,u1))}

jakstatqq<-function(u1,v1){
N1<-length(u1)
N2<-length(v1)
u<-matrix(rep(u1,N2),ncol=N2,nrow=N1,byrow=FALSE)
v<-matrix(rep(v1,N1),ncol=N2,nrow=N1,byrow=TRUE)
A<-matrix(jakstat$a,nrow=N1,ncol=N2)
((4*jakstat$lam^2)*(pmin(u-jakstat$lam,v-jakstat$lam)-(A+jakstat$lam))*as.numeric(pmin(u-jakstat$lam,v-jakstat$lam)>(A+jakstat$lam))+(2*jakstat$lam)*((v+jakstat$lam)*pmin(u-jakstat$lam,v+jakstat$lam) - 0.5*pmin(u-jakstat$lam,v+jakstat$lam)^2-(v+jakstat$lam)*pmax(A+jakstat$lam,v-jakstat$lam) + 0.5*pmax(A+jakstat$lam,v-jakstat$lam)^2)*as.numeric(pmin(u-jakstat$lam,v+jakstat$lam)>pmax(A+jakstat$lam,v-jakstat$lam))+((1/3)*pmin(A+jakstat$lam,pmin(u-jakstat$lam,v-jakstat$lam))^3 + (jakstat$lam-jakstat$a)*pmin(A+jakstat$lam,pmin(u-jakstat$lam,v-jakstat$lam))^2 + (jakstat$lam-jakstat$a)^2*pmin(A+jakstat$lam,pmin(u-jakstat$lam,v-jakstat$lam))-(1/3)*(A-jakstat$lam)^3 - (jakstat$lam-jakstat$a)*(jakstat$a-jakstat$lam)^2 - (jakstat$lam-jakstat$a)^2*(jakstat$a-jakstat$lam))*as.numeric(pmin(A+jakstat$lam,pmin(u-jakstat$lam,v-jakstat$lam))>(jakstat$a-jakstat$lam))+(v-jakstat$a)*(0.5*pmin(A+jakstat$lam,u-jakstat$lam)^2 + (jakstat$lam-jakstat$a)*pmin(A+jakstat$lam,u-jakstat$lam)- 0.5*pmax(A-jakstat$lam,v-jakstat$lam)^2 - (jakstat$lam-jakstat$a)*pmax(A-jakstat$lam,v-jakstat$lam))*as.numeric(pmin(A+jakstat$lam,u-jakstat$lam)>pmax(A-jakstat$lam,v-jakstat$lam))+(2*jakstat$lam)*((u+jakstat$lam)*pmin(u+jakstat$lam,v-jakstat$lam)-0.5*pmin(u+jakstat$lam,v-jakstat$lam)^2-(u+jakstat$lam)*pmax(A+jakstat$lam,u-jakstat$lam) + 0.5*pmax(A+jakstat$lam,u-jakstat$lam)^2)*as.numeric(pmin(u+jakstat$lam,v-jakstat$lam)>pmax(A+jakstat$lam,u-jakstat$lam))+((u+jakstat$lam)*(v+jakstat$lam)*pmin(u+jakstat$lam,v+jakstat$lam) - 0.5*(u+v+2*jakstat$lam)*pmin(u+jakstat$lam,v+jakstat$lam)^2 + (1/3)*pmin(u+jakstat$lam,v+jakstat$lam)^3 - (u+jakstat$lam)*(v+jakstat$lam)*pmax(A+jakstat$lam,pmax(u-jakstat$lam,v-jakstat$lam)) + 0.5*(u+v+2*jakstat$lam)*pmax(A+jakstat$lam,pmax(u-jakstat$lam,v-jakstat$lam))^2 - (1/3)*pmax(A+jakstat$lam,pmax(u-jakstat$lam,v-jakstat$lam))^3)*as.numeric(pmin(u+jakstat$lam,v+jakstat$lam)>pmax(A+jakstat$lam,pmax(u-jakstat$lam,v-jakstat$lam)))+(u-jakstat$a)*(0.5*pmin(A+jakstat$lam,v-jakstat$lam)^2 + (jakstat$lam-jakstat$a)*pmin(A+jakstat$lam,v-jakstat$lam)- 0.5*pmax(A-jakstat$lam,u-jakstat$lam)^2 - (jakstat$lam-jakstat$a)*pmax(A-jakstat$lam,u-jakstat$lam))*as.numeric(pmin(jakstat$a+jakstat$lam,v-jakstat$lam)>pmax(jakstat$a-jakstat$lam,u-jakstat$lam))+(u-jakstat$a)*(v-jakstat$a)*((jakstat$a+jakstat$lam)-pmax(u-jakstat$lam,v-jakstat$lam))*as.numeric((A+jakstat$lam)>pmax(u-jakstat$lam,v-jakstat$lam)))}


jakstatsolver2.cpp <- function(S, THETA2, U0, TAU, qrbeval, QRB, SAPV, epo, randoms1, randoms2, randoms3, randoms4){
	.Call( "jakstatsolver2cpp", S, THETA2, U0, TAU, qrbeval, QRB, SAPV, epo, randoms1, randoms2, randoms3, randoms4, PACKAGE = "aceodes" )
}

jakstatnsel4.cpp <- function(THETA, mu, v, y){
	.Call( "jakstatnsel4cpp", THETA, mu, v, y, PACKAGE = "aceodes" )
}

jakstatsig.cpp <- function(mu1, mu2, v, y){
	.Call( "jakstatsigcpp", mu1, mu2, v, y, PACKAGE = "aceodes" )
}

jakstatnael.cpp <- function(THETA, mu, v, y, order){
	.Call( "jakstatnaelcpp", THETA, mu, v, y, order, PACKAGE = "aceodes" )
}

JakstatNael<-function(d, B){

n0<-dim(d)[1]-1  
BB<-B+B

sam<-sample(x=1:dim(jakstat$mcmc)[1],size=BB,replace=TRUE)

sig<-cbind(runif(n=BB,min=0,max=0.1)^2,runif(n=BB,min=0,max=0.1)^2,runif(n=BB,min=0,max=20)^2,runif(n=BB,min=0,max=0.1)^2)
THETA2<-jakstat$mcmc[sam,1:6]
TAU<-jakstat$mcmc[sam,7]
U0<-jakstat$mcmc[sam,8]

d0<-d[-1,]
grid<-sort((as.vector(d0)+1)*30)
extra<-(d[1,]+1)*30
grid<-c(grid,extra)

epo<-exp(matrix(rnorm(n=BB*jakstat$N0,mean=rep(jakstat$q,BB),sd=rep(sqrt(jakstat$R),BB)),nrow=BB,byrow=TRUE))

qrbeval<-jakstatqr(grid,jakstat$S[1:jakstat$N0])%*%(jakstat$Binv/jakstat$alpha)
qqqq<-jakstatqq(grid,grid)/jakstat$alpha
rqrq<-jakstatrq(jakstat$S[1:jakstat$N0],grid)/jakstat$alpha
sCCC<-diag(qqqq-qrbeval%*%rqrq)
sCCC<-ifelse(is.na(sCCC),0,sCCC)
sCCC<-diag(sqrt(ifelse(sCCC<0,0,sCCC)))

mu<-jakstatsolver2.cpp(S=jakstat$S,THETA2=THETA2,U0=U0,TAU=TAU,qrbeval=qrbeval,QRB=jakstat$QRB,SAPV=jakstat$SAPV,epo=epo,randoms1=matrix(rnorm(n=length(jakstat$SAPV)*BB),ncol=BB),randoms2=matrix(rnorm(n=length(jakstat$SAPV)*BB),ncol=BB),randoms3=matrix(rnorm(n=length(jakstat$SAPV)*BB),ncol=BB),randoms4=matrix(rnorm(n=length(jakstat$SAPV)*BB),ncol=BB))
mu[is.na(mu)]<-0
mu<-mu+matrix(as.vector(sCCC%*%matrix(rnorm(n=3*(n0+1)*BB),nrow=n0+1)),byrow=TRUE,nrow=BB)

g1<-matrix(rep(THETA2[,5],n0),ncol=n0,byrow=FALSE)*(mu[,(n0+2):(2*n0+1)]+2*mu[,(2*n0+3):(3*n0+2)])
g2<-matrix(rep(THETA2[,6],n0),ncol=n0,byrow=FALSE)*(mu[,1:n0]+mu[,(n0+2):(2*n0+1)]+2*mu[,(2*n0+3):(3*n0+2)])
g3<-U0
g4<-mu[,3*n0+3]/(mu[,(2*n0+2)]+mu[,3*n0+3])

MU<-cbind(g1,g2,g3,g4)
V<-cbind(matrix(rep(sig[,1],n0),ncol=n0,byrow=FALSE),matrix(rep(sig[,2],n0),ncol=n0,byrow=FALSE),sig[,3],sig[,4])

y<-MU[1:B,]+sqrt(V[1:B,])*matrix(rnorm(n=B*(2*n0+2)),ncol=2*n0+2)

#out<-nsel4.cpp(THETA=THETA2[-(1:B),],mu=MU[-(1:B),],v=V[-(1:B),],y=y) 
out<-jakstatnael.cpp(THETA=THETA2[-(1:B),],mu=MU[-(1:B),],v=V[-(1:B),],y=y,order=apply(THETA2[-(1:B),],2,order))

#out2<-apply((THETA2[1:B,]-out)^2,1,sum)
out2<-apply(abs(THETA2[1:B,]-out),1,sum)

#out2<-sig.cpp(mu1=MU[1:B,],mu2=MU[-(1:B),],v=V[-(1:B),],y=y) 
out2[is.na(out2)]<-0.48 # Set any NAs to approximate performance of original design

-out2}

JakstatNsel<-function(d, B){

n0<-dim(d)[1]-1  
BB<-B+B

sam<-sample(x=1:dim(jakstat$mcmc)[1],size=BB,replace=TRUE)

sig<-cbind(runif(n=BB,min=0,max=0.1)^2,runif(n=BB,min=0,max=0.1)^2,runif(n=BB,min=0,max=20)^2,runif(n=BB,min=0,max=0.1)^2)
THETA2<-jakstat$mcmc[sam,1:6]
TAU<-jakstat$mcmc[sam,7]
U0<-jakstat$mcmc[sam,8]

d0<-d[-1,]
grid<-sort((as.vector(d0)+1)*30)
extra<-(d[1,]+1)*30
grid<-c(grid,extra)

epo<-exp(matrix(rnorm(n=BB*jakstat$N0,mean=rep(jakstat$q,BB),sd=rep(sqrt(jakstat$R),BB)),nrow=BB,byrow=TRUE))

qrbeval<-jakstatqr(grid,jakstat$S[1:jakstat$N0])%*%(jakstat$Binv/jakstat$alpha)
qqqq<-jakstatqq(grid,grid)/jakstat$alpha
rqrq<-jakstatrq(jakstat$S[1:jakstat$N0],grid)/jakstat$alpha
sCCC<-diag(qqqq-qrbeval%*%rqrq)
sCCC<-ifelse(is.na(sCCC),0,sCCC)
sCCC<-diag(sqrt(ifelse(sCCC<0,0,sCCC)))

mu<-jakstatsolver2.cpp(S=jakstat$S,THETA2=THETA2,U0=U0,TAU=TAU,qrbeval=qrbeval,QRB=jakstat$QRB,SAPV=jakstat$SAPV,epo=epo,randoms1=matrix(rnorm(n=length(jakstat$SAPV)*BB),ncol=BB),randoms2=matrix(rnorm(n=length(jakstat$SAPV)*BB),ncol=BB),randoms3=matrix(rnorm(n=length(jakstat$SAPV)*BB),ncol=BB),randoms4=matrix(rnorm(n=length(jakstat$SAPV)*BB),ncol=BB))
mu[is.na(mu)]<-0
mu<-mu+matrix(as.vector(sCCC%*%matrix(rnorm(n=3*(n0+1)*BB),nrow=n0+1)),byrow=TRUE,nrow=BB)

g1<-matrix(rep(THETA2[,5],n0),ncol=n0,byrow=FALSE)*(mu[,(n0+2):(2*n0+1)]+2*mu[,(2*n0+3):(3*n0+2)])
g2<-matrix(rep(THETA2[,6],n0),ncol=n0,byrow=FALSE)*(mu[,1:n0]+mu[,(n0+2):(2*n0+1)]+2*mu[,(2*n0+3):(3*n0+2)])
g3<-U0
g4<-mu[,3*n0+3]/(mu[,(2*n0+2)]+mu[,3*n0+3])

MU<-cbind(g1,g2,g3,g4)
V<-cbind(matrix(rep(sig[,1],n0),ncol=n0,byrow=FALSE),matrix(rep(sig[,2],n0),ncol=n0,byrow=FALSE),sig[,3],sig[,4])

y<-MU[1:B,]+sqrt(V[1:B,])*matrix(rnorm(n=B*(2*n0+2)),ncol=2*n0+2)

out<-jakstatnsel4.cpp(THETA=THETA2[-(1:B),],mu=MU[-(1:B),],v=V[-(1:B),],y=y) 
#out<-nael.cpp(THETA=THETA2[-(1:B),],mu=MU[-(1:B),],v=V[-(1:B),],y=y,order=apply(THETA2[-(1:B),],2,order))

out2<-apply((THETA2[1:B,]-out)^2,1,sum)
#out2<-apply(abs(THETA2[1:B,]-out),1,sum)

#out2<-sig.cpp(mu1=MU[1:B,],mu2=MU[-(1:B),],v=V[-(1:B),],y=y) 
out2[is.na(out2)]<-0.25 # Set any NAs to approximate performance of original design

-out2}

JakstatSig<-function(d, B){

n0<-dim(d)[1]-1  
BB<-B+B

sam<-sample(x=1:dim(jakstat$mcmc)[1],size=BB,replace=TRUE)

sig<-cbind(runif(n=BB,min=0,max=0.1)^2,runif(n=BB,min=0,max=0.1)^2,runif(n=BB,min=0,max=20)^2,runif(n=BB,min=0,max=0.1)^2)
THETA2<-jakstat$mcmc[sam,1:6]
TAU<-jakstat$mcmc[sam,7]
U0<-jakstat$mcmc[sam,8]

d0<-d[-1,]
grid<-sort((as.vector(d0)+1)*30)
extra<-(d[1,]+1)*30
grid<-c(grid,extra)

epo<-exp(matrix(rnorm(n=BB*jakstat$N0,mean=rep(jakstat$q,BB),sd=rep(sqrt(jakstat$R),BB)),nrow=BB,byrow=TRUE))

qrbeval<-jakstatqr(grid,jakstat$S[1:jakstat$N0])%*%(jakstat$Binv/jakstat$alpha)
qqqq<-jakstatqq(grid,grid)/jakstat$alpha
rqrq<-jakstatrq(jakstat$S[1:jakstat$N0],grid)/jakstat$alpha
sCCC<-diag(qqqq-qrbeval%*%rqrq)
sCCC<-ifelse(is.na(sCCC),0,sCCC)
sCCC<-diag(sqrt(ifelse(sCCC<0,0,sCCC)))

mu<-jakstatsolver2.cpp(S=jakstat$S,THETA2=THETA2,U0=U0,TAU=TAU,qrbeval=qrbeval,QRB=jakstat$QRB,SAPV=jakstat$SAPV,epo=epo,randoms1=matrix(rnorm(n=length(jakstat$SAPV)*BB),ncol=BB),randoms2=matrix(rnorm(n=length(jakstat$SAPV)*BB),ncol=BB),randoms3=matrix(rnorm(n=length(jakstat$SAPV)*BB),ncol=BB),randoms4=matrix(rnorm(n=length(jakstat$SAPV)*BB),ncol=BB))
mu[is.na(mu)]<-0
mu<-mu+matrix(as.vector(sCCC%*%matrix(rnorm(n=3*(n0+1)*BB),nrow=n0+1)),byrow=TRUE,nrow=BB)

g1<-matrix(rep(THETA2[,5],n0),ncol=n0,byrow=FALSE)*(mu[,(n0+2):(2*n0+1)]+2*mu[,(2*n0+3):(3*n0+2)])
g2<-matrix(rep(THETA2[,6],n0),ncol=n0,byrow=FALSE)*(mu[,1:n0]+mu[,(n0+2):(2*n0+1)]+2*mu[,(2*n0+3):(3*n0+2)])
g3<-U0
g4<-mu[,3*n0+3]/(mu[,(2*n0+2)]+mu[,3*n0+3])

MU<-cbind(g1,g2,g3,g4)
V<-cbind(matrix(rep(sig[,1],n0),ncol=n0,byrow=FALSE),matrix(rep(sig[,2],n0),ncol=n0,byrow=FALSE),sig[,3],sig[,4])

y<-MU[1:B,]+sqrt(V[1:B,])*matrix(rnorm(n=B*(2*n0+2)),ncol=2*n0+2)

#out<-nsel4.cpp(THETA=THETA2[-(1:B),],mu=MU[-(1:B),],v=V[-(1:B),],y=y) 
#out<-nael.cpp(THETA=THETA2[-(1:B),],mu=MU[-(1:B),],v=V[-(1:B),],y=y,order=apply(THETA2[-(1:B),],2,order))

#out2<-apply((THETA2[1:B,]-out)^2,1,sum)
#out2<-apply(abs(THETA2[1:B,]-out),1,sum)

out2<-jakstatsig.cpp(mu1=MU[1:B,],mu2=MU[-(1:B),],v=V[-(1:B),],y=y) 
out2[is.na(out2)]<-5.4 # Set any NAs to approximate performance of original design

out2}

JakstatLimits<-function(i,j,d){
hered<-d
if(i>1){
x<-c(-1,sort(runif(9998))*2-1,1)
exc<-as.vector(hered[-1,])
n<-length(exc)
for(t in (1:n)[exc!=hered[i,j]]){
int<-exc[t]+c(-1,1)*(1/30)
x<-x[x<=int[1] | x>=int[2]]}} else{
x<-c(-1,sort(runif(9998))*2-1,1)}
x}

JakstatTerminal<-function(rep, utility = "SIG", scaled = TRUE){

if(utility=="SIG"){
des<-as.vector(jakstatopt$sig.designs[jakstatopt$sig.designs[,1]==rep,-1])}
if(utility=="NSEL"){
des<-as.vector(jakstatopt$nsel.designs[jakstatopt$nsel.designs[,1]==rep,-1])}
if(utility=="NAEL"){
des<-as.vector(jakstatopt$nael.designs[jakstatopt$nael.designs[,1]==rep,-1])}

if(!scaled){
des<-30*(des+1)} 

des} 

JakstatEvals<-function(rep = NULL, utility = "SIG"){

if(!is.null(rep)){	
if(rep<0|rep>20){
stop("rep should be 1, 2, ..., 20")}} 
	
if(utility=="SIG"){
out<-jakstatopt$sig.out}
if(utility=="NSEL"){
out<-jakstatopt$nsel.out}
if(utility=="NAEL"){
out<-jakstatopt$nael.out}

if(!is.null(rep)){
out<-out[rep,]}
	
out}

JakstatDefaultEvals<-function(utility = "SIG"){

if(utility=="SIG"){
out<-jakstatopt$sig.def}
if(utility=="NSEL"){
out<-jakstatopt$nsel.def}
if(utility=="NAEL"){
out<-jakstatopt$nael.def}

out}


JakstatFinal<-function(utility = "SIG", scaled = TRUE){

if(utility=="SIG"){
des<-as.vector(jakstatopt$sig.design)}
if(utility=="NSEL"){
des<-as.vector(jakstatopt$nsel.design)}
if(utility=="NAEL"){
des<-as.vector(jakstatopt$nael.design)}

if(!scaled){
des<-30*(des+1)} 

des}  

JakstatSolver<-function(grid, B){

BB<-B  
n0<-length(grid)
sam<-sample(x=1:dim(jakstat$mcmc)[1],size=B,replace=TRUE)
THETA2<-jakstat$mcmc[sam,1:6]
TAU<-jakstat$mcmc[sam,7]
U0<-jakstat$mcmc[sam,8]

epo<-exp(matrix(rnorm(n=B*jakstat$N0,mean=rep(jakstat$q,B),sd=rep(sqrt(jakstat$R),B)),nrow=B,byrow=TRUE))

qrbeval<-jakstatqr(grid,jakstat$S[1:jakstat$N0])%*%(jakstat$Binv/jakstat$alpha)
qqqq<-jakstatqq(grid,grid)/jakstat$alpha
rqrq<-jakstatrq(jakstat$S[1:jakstat$N0],grid)/jakstat$alpha
sCCC<-diag(qqqq-qrbeval%*%rqrq)
sCCC<-ifelse(is.na(sCCC),0,sCCC)
sCCC<-diag(sqrt(ifelse(sCCC<0,0,sCCC)))

mu<-jakstatsolver2.cpp(S=jakstat$S,THETA2=THETA2,U0=U0,TAU=TAU,qrbeval=qrbeval,QRB=jakstat$QRB,SAPV=jakstat$SAPV,epo=epo,randoms1=matrix(rnorm(n=length(jakstat$SAPV)*BB),ncol=BB),randoms2=matrix(rnorm(n=length(jakstat$SAPV)*BB),ncol=BB),randoms3=matrix(rnorm(n=length(jakstat$SAPV)*BB),ncol=BB),randoms4=matrix(rnorm(n=length(jakstat$SAPV)*BB),ncol=BB))
mu[is.na(mu)]<-0
mu<-mu+matrix(as.vector(sCCC%*%matrix(rnorm(n=3*n0*B),nrow=n0)),byrow=TRUE,nrow=B)

g1<-matrix(rep(THETA2[,5],n0),ncol=n0,byrow=FALSE)*(mu[,(n0+1):(2*n0)]+2*mu[,(2*n0+1):(3*n0)])
g2<-matrix(rep(THETA2[,6],n0),ncol=n0,byrow=FALSE)*(mu[,1:n0]+mu[,(n0+1):(2*n0)]+2*mu[,(2*n0+1):(3*n0)])
g3<-U0
g4<-mu[,(2*n0+1):(3*n0)]/(mu[,(n0+1):(2*n0)]+mu[,(2*n0+1):(3*n0)])

list(G1=g1,G2=g2,G3=g3,G4=g4)
}	

##################################################################################################

### FITZHUGH-NAGUMO #####################################################################################

##################################################################################################

fitzhughlimits2<-function(d){
  x<-c(-1,sort(runif(9998))*2-1,1)
  n<-dim(d)[1]
  for(t in (1:n)){
    int<-d[t,]+c(-1,1)*(0.25/10)
    x<-x[x<int[1] | x>int[2]]}
  x}

FitzhughStart<-function(){
  start.d<-matrix(runif(1)*2-1,ncol=1)
  for(j in 2:21){
    cand<-fitzhughlimits2(d=start.d)
    start.d<-rbind(start.d,cand[sample(x=1:length(cand),size=1)])
    start.d<-matrix(sort(as.vector(start.d)),ncol=1)}
start.d}

fitzhughrr<-function(u1,v1){
N1<-length(u1)
N2<-length(v1)
u<-matrix(rep(u1,N2),ncol=N2,nrow=N1,byrow=FALSE)
v<-matrix(rep(v1,N1),ncol=N2,nrow=N1,byrow=TRUE)
ret<-ifelse(u>=v,v,u)-ifelse(u<=v,v,u)
(ret+2*fitzhugh$lam)*ifelse(ret>=(-2*fitzhugh$lam),1,0)}

fitzhughrq<-function(u1,v1){
N1<-length(u1)
N2<-length(v1)
u<-matrix(rep(u1,N2),ncol=N2,nrow=N1,byrow=FALSE)
v<-matrix(rep(v1,N1),ncol=N2,nrow=N1,byrow=TRUE)
A<-matrix(fitzhugh$a,nrow=N1,ncol=N2)
2*fitzhugh$lam*(pmin(v-fitzhugh$lam,u+fitzhugh$lam)-pmax(A+fitzhugh$lam,u-fitzhugh$lam))*as.numeric(pmin(v-fitzhugh$lam,u+fitzhugh$lam)>pmax(A+fitzhugh$lam,u-fitzhugh$lam))+(0.5*pmin(A+fitzhugh$lam,pmin(v-fitzhugh$lam,u+fitzhugh$lam))^2+(fitzhugh$lam-fitzhugh$a)*pmin(A+fitzhugh$lam,pmin(v-fitzhugh$lam,u+fitzhugh$lam))-0.5*pmax(A-fitzhugh$lam,u-fitzhugh$lam)^2 -(fitzhugh$lam-fitzhugh$a)*pmax(A-fitzhugh$lam,u-fitzhugh$lam))*as.numeric(pmin(A+fitzhugh$lam,pmin(v-fitzhugh$lam,u+fitzhugh$lam))>pmax(A-fitzhugh$lam,u-fitzhugh$lam))+((v+fitzhugh$lam)*pmin(v+fitzhugh$lam,u+fitzhugh$lam) - 0.5*pmin(v+fitzhugh$lam,u+fitzhugh$lam)^2-(v+fitzhugh$lam)*pmax(A+fitzhugh$lam,pmax(v-fitzhugh$lam,u-fitzhugh$lam))+ 0.5*pmax(A+fitzhugh$lam,pmax(v-fitzhugh$lam,u-fitzhugh$lam))^2)*as.numeric(pmin(v+fitzhugh$lam,u+fitzhugh$lam)>pmax(A+fitzhugh$lam,pmax(v-fitzhugh$lam,u-fitzhugh$lam)))+ ((v-fitzhugh$a)*pmin(A+fitzhugh$lam,u+fitzhugh$lam) - (v-fitzhugh$a)*pmax(v-fitzhugh$lam,u-fitzhugh$lam))*as.numeric(pmin(A+fitzhugh$lam,u+fitzhugh$lam)>pmax(v-fitzhugh$lam,u-fitzhugh$lam))
}

fitzhughqr<-function(u1,v1){
t(fitzhughrq(v1,u1))}

fitzhughqq<-function(u1,v1){
N1<-length(u1)
N2<-length(v1)
u<-matrix(rep(u1,N2),ncol=N2,nrow=N1,byrow=FALSE)
v<-matrix(rep(v1,N1),ncol=N2,nrow=N1,byrow=TRUE)
A<-matrix(fitzhugh$a,nrow=N1,ncol=N2)
((4*fitzhugh$lam^2)*(pmin(u-fitzhugh$lam,v-fitzhugh$lam)-(A+fitzhugh$lam))*as.numeric(pmin(u-fitzhugh$lam,v-fitzhugh$lam)>(A+fitzhugh$lam))+(2*fitzhugh$lam)*((v+fitzhugh$lam)*pmin(u-fitzhugh$lam,v+fitzhugh$lam) - 0.5*pmin(u-fitzhugh$lam,v+fitzhugh$lam)^2-(v+fitzhugh$lam)*pmax(A+fitzhugh$lam,v-fitzhugh$lam) + 0.5*pmax(A+fitzhugh$lam,v-fitzhugh$lam)^2)*as.numeric(pmin(u-fitzhugh$lam,v+fitzhugh$lam)>pmax(A+fitzhugh$lam,v-fitzhugh$lam))    +((1/3)*pmin(A+fitzhugh$lam,pmin(u-fitzhugh$lam,v-fitzhugh$lam))^3 + (fitzhugh$lam-fitzhugh$a)*pmin(A+fitzhugh$lam,pmin(u-fitzhugh$lam,v-fitzhugh$lam))^2 + (fitzhugh$lam-fitzhugh$a)^2*pmin(A+fitzhugh$lam,pmin(u-fitzhugh$lam,v-fitzhugh$lam))-(1/3)*(A-fitzhugh$lam)^3 - (fitzhugh$lam-fitzhugh$a)*(fitzhugh$a-fitzhugh$lam)^2 - (fitzhugh$lam-fitzhugh$a)^2*(fitzhugh$a-fitzhugh$lam))*as.numeric(pmin(A+fitzhugh$lam,pmin(u-fitzhugh$lam,v-fitzhugh$lam))>(fitzhugh$a-fitzhugh$lam))+(v-fitzhugh$a)*(0.5*pmin(A+fitzhugh$lam,u-fitzhugh$lam)^2 + (fitzhugh$lam-fitzhugh$a)*pmin(A+fitzhugh$lam,u-fitzhugh$lam)- 0.5*pmax(A-fitzhugh$lam,v-fitzhugh$lam)^2 - (fitzhugh$lam-fitzhugh$a)*pmax(A-fitzhugh$lam,v-fitzhugh$lam))*as.numeric(pmin(A+fitzhugh$lam,u-fitzhugh$lam)>pmax(A-fitzhugh$lam,v-fitzhugh$lam))+(2*fitzhugh$lam)*((u+fitzhugh$lam)*pmin(u+fitzhugh$lam,v-fitzhugh$lam)-0.5*pmin(u+fitzhugh$lam,v-fitzhugh$lam)^2-(u+fitzhugh$lam)*pmax(A+fitzhugh$lam,u-fitzhugh$lam) + 0.5*pmax(A+fitzhugh$lam,u-fitzhugh$lam)^2)*as.numeric(pmin(u+fitzhugh$lam,v-fitzhugh$lam)>pmax(A+fitzhugh$lam,u-fitzhugh$lam))+((u+fitzhugh$lam)*(v+fitzhugh$lam)*pmin(u+fitzhugh$lam,v+fitzhugh$lam) - 0.5*(u+v+2*fitzhugh$lam)*pmin(u+fitzhugh$lam,v+fitzhugh$lam)^2 + (1/3)*pmin(u+fitzhugh$lam,v+fitzhugh$lam)^3 - (u+fitzhugh$lam)*(v+fitzhugh$lam)*pmax(A+fitzhugh$lam,pmax(u-fitzhugh$lam,v-fitzhugh$lam)) + 0.5*(u+v+2*fitzhugh$lam)*pmax(A+fitzhugh$lam,pmax(u-fitzhugh$lam,v-fitzhugh$lam))^2 - (1/3)*pmax(A+fitzhugh$lam,pmax(u-fitzhugh$lam,v-fitzhugh$lam))^3)*as.numeric(pmin(u+fitzhugh$lam,v+fitzhugh$lam)>pmax(A+fitzhugh$lam,pmax(u-fitzhugh$lam,v-fitzhugh$lam)))+(u-fitzhugh$a)*(0.5*pmin(A+fitzhugh$lam,v-fitzhugh$lam)^2 + (fitzhugh$lam-fitzhugh$a)*pmin(A+fitzhugh$lam,v-fitzhugh$lam)- 0.5*pmax(A-fitzhugh$lam,u-fitzhugh$lam)^2 - (fitzhugh$lam-fitzhugh$a)*pmax(A-fitzhugh$lam,u-fitzhugh$lam))*as.numeric(pmin(fitzhugh$a+fitzhugh$lam,v-fitzhugh$lam)>pmax(fitzhugh$a-fitzhugh$lam,u-fitzhugh$lam))+(u-fitzhugh$a)*(v-fitzhugh$a)*((fitzhugh$a+fitzhugh$lam)-pmax(u-fitzhugh$lam,v-fitzhugh$lam))*as.numeric((A+fitzhugh$lam)>pmax(u-fitzhugh$lam,v-fitzhugh$lam)))}

fitzsolver2.cpp <- function(S, THETA2, qrbeval, QRB, SAPV, randoms1, randoms2){
  .Call( "fitzsolver2cpp", S, THETA2, qrbeval, QRB, SAPV, randoms1, randoms2, PACKAGE = "aceodes" )
}

fitzhughnsel.cpp <- function(THETA, mu, v, y){
  .Call( "fitzhughnsel4cpp", THETA, mu, v, y, PACKAGE = "aceodes" )
}

fitzhughsig.cpp <- function(mu1, mu2, v, y){
  .Call( "fitzhughsigcpp", mu1, mu2, v, y, PACKAGE = "aceodes" )
}

fitzhughnael.cpp <- function(THETA, mu, v, y, order){
  .Call( "fitzhughnaelcpp", THETA, mu, v, y, order, PACKAGE = "aceodes" )
}

FitzhughNael<-function(d, B){

n0<-dim(d)[1]  
BB<-B+B

THETA2<-cbind(runif(n=BB,min=0,max=1),runif(n=BB,min=0,max=1),runif(n=BB,min=1,max=5))
sig<-runif(n=BB,min=0.5,max=1)^2

grid<-(as.vector(d)+1)*0.5*20

qrbeval<-fitzhughqr(grid,fitzhugh$S[1:fitzhugh$N0])%*%(fitzhugh$Binv/fitzhugh$alpha)
qqqq<-fitzhughqq(grid,grid)/fitzhugh$alpha
rqrq<-fitzhughrq(fitzhugh$S[1:fitzhugh$N0],grid)/fitzhugh$alpha
sCCC<-diag(qqqq-qrbeval%*%rqrq)
sCCC<-ifelse(is.na(sCCC),0,sCCC)
sCCC<-diag(sqrt(ifelse(sCCC<0,0,sCCC)))

mu<-fitzsolver2.cpp(S=fitzhugh$S,THETA2=THETA2,qrbeval=qrbeval,QRB=fitzhugh$QRB,SAPV=fitzhugh$SAPV,randoms1 = matrix(rnorm(BB*length(fitzhugh$SAPV)),ncol=BB),randoms2 = matrix(rnorm(BB*length(fitzhugh$SAPV)),ncol=BB))
mu<-mu+matrix(as.vector(sCCC%*%matrix(rnorm(n=n0*BB),nrow=n0)),byrow=TRUE,nrow=BB)
V<-matrix(rep(sig,n0),nrow=BB,byrow=FALSE)
y<-mu[1:B,]+sqrt(V[1:B,])*matrix(rnorm(n=B*n0),ncol=n0)

#out<-nsel4.cpp(THETA=THETA2[-(1:B),],mu=mu[-(1:B),],v=V[-(1:B),],y=y) 
out<-fitzhughnael.cpp(THETA=THETA2[-(1:B),],mu=mu[-(1:B),],v=V[-(1:B),],y=y,order=apply(THETA2[-(1:B),],2,order))

#out2<-apply((THETA2[1:B,]-out)^2,1,sum)
out2<-apply(abs(THETA2[1:B,]-out),1,sum)

#out2<-sig.cpp(mu1=mu[1:B,],mu2=mu[-(1:B),],v=V[-(1:B),],y=y) 
out2[is.na(out2)]<-0.66 # Set any NAs to approximate performance of uniform design

-out2}

FitzhughNsel<-function(d, B){

n0<-dim(d)[1]  
BB<-B+B

THETA2<-cbind(runif(n=BB,min=0,max=1),runif(n=BB,min=0,max=1),runif(n=BB,min=1,max=5))
sig<-runif(n=BB,min=0.5,max=1)^2

grid<-(as.vector(d)+1)*0.5*20

qrbeval<-fitzhughqr(grid,fitzhugh$S[1:fitzhugh$N0])%*%(fitzhugh$Binv/fitzhugh$alpha)
qqqq<-fitzhughqq(grid,grid)/fitzhugh$alpha
rqrq<-fitzhughrq(fitzhugh$S[1:fitzhugh$N0],grid)/fitzhugh$alpha
sCCC<-diag(qqqq-qrbeval%*%rqrq)
sCCC<-ifelse(is.na(sCCC),0,sCCC)
sCCC<-diag(sqrt(ifelse(sCCC<0,0,sCCC)))

mu<-fitzsolver2.cpp(S=fitzhugh$S,THETA2=THETA2,qrbeval=qrbeval,QRB=fitzhugh$QRB,SAPV=fitzhugh$SAPV,randoms1 = matrix(rnorm(BB*length(fitzhugh$SAPV)),ncol=BB),randoms2 = matrix(rnorm(BB*length(fitzhugh$SAPV)),ncol=BB))
mu<-mu+matrix(as.vector(sCCC%*%matrix(rnorm(n=n0*BB),nrow=n0)),byrow=TRUE,nrow=BB)
V<-matrix(rep(sig,n0),nrow=BB,byrow=FALSE)
y<-mu[1:B,]+sqrt(V[1:B,])*matrix(rnorm(n=B*n0),ncol=n0)

out<-fitzhughnsel.cpp(THETA=THETA2[-(1:B),],mu=mu[-(1:B),],v=V[-(1:B),],y=y) 
#out<-nael.cpp(THETA=THETA2[-(1:B),],mu=mu[-(1:B),],v=V[-(1:B),],y=y,order=apply(THETA2[-(1:B),],2,order))

out2<-apply((THETA2[1:B,]-out)^2,1,sum)
#out2<-apply(abs(THETA2[1:B,]-out),1,sum)

#out2<-sig.cpp(mu1=mu[1:B,],mu2=mu[-(1:B),],v=V[-(1:B),],y=y) 
out2[is.na(out2)]<-0.36 # Set any NAs to approximate performance of uniform design

-out2}

FitzhughSig<-function(d, B){

n0<-dim(d)[1]  
BB<-B+B

THETA2<-cbind(runif(n=BB,min=0,max=1),runif(n=BB,min=0,max=1),runif(n=BB,min=1,max=5))
sig<-runif(n=BB,min=0.5,max=1)^2

grid<-(as.vector(d)+1)*0.5*20

qrbeval<-fitzhughqr(grid,fitzhugh$S[1:fitzhugh$N0])%*%(fitzhugh$Binv/fitzhugh$alpha)
qqqq<-fitzhughqq(grid,grid)/fitzhugh$alpha
rqrq<-fitzhughrq(fitzhugh$S[1:fitzhugh$N0],grid)/fitzhugh$alpha
sCCC<-diag(qqqq-qrbeval%*%rqrq)
sCCC<-ifelse(is.na(sCCC),0,sCCC)
sCCC<-diag(sqrt(ifelse(sCCC<0,0,sCCC)))

mu<-fitzsolver2.cpp(S=fitzhugh$S,THETA2=THETA2,qrbeval=qrbeval,QRB=fitzhugh$QRB,SAPV=fitzhugh$SAPV,randoms1 = matrix(rnorm(BB*length(fitzhugh$SAPV)),ncol=BB),randoms2 = matrix(rnorm(BB*length(fitzhugh$SAPV)),ncol=BB))
mu<-mu+matrix(as.vector(sCCC%*%matrix(rnorm(n=n0*BB),nrow=n0)),byrow=TRUE,nrow=BB)
V<-matrix(rep(sig,n0),nrow=BB,byrow=FALSE)
y<-mu[1:B,]+sqrt(V[1:B,])*matrix(rnorm(n=B*n0),ncol=n0)

#out<-nsel4.cpp(THETA=THETA2[-(1:B),],mu=mu[-(1:B),],v=V[-(1:B),],y=y) 
#out<-nael.cpp(THETA=THETA2[-(1:B),],mu=mu[-(1:B),],v=V[-(1:B),],y=y,order=apply(THETA2[-(1:B),],2,order))

#out2<-apply((THETA2[1:B,]-out)^2,1,sum)
#out2<-apply(abs(THETA2[1:B,]-out),1,sum)

out2<-fitzhughsig.cpp(mu1=mu[1:B,],mu2=mu[-(1:B),],v=V[-(1:B),],y=y) 
out2[is.na(out2)]<-3.3 # Set any NAs to approximate performance of uniform design

out2}

FitzhughSolver<-function(grid, B){
theta<-cbind(runif(n=B,min=0,max=1),runif(n=B,min=0,max=1),runif(n=B,min=1,max=5))
qrbeval<-fitzhughqr(grid,fitzhugh$S[1:fitzhugh$N0])%*%(fitzhugh$Binv/fitzhugh$alpha)
qqqq<-fitzhughqq(grid,grid)/fitzhugh$alpha
rqrq<-fitzhughrq(fitzhugh$S[1:fitzhugh$N0],grid)/fitzhugh$alpha
sCCC<-diag(qqqq-qrbeval%*%rqrq)
sCCC<-ifelse(is.na(sCCC),0,sCCC)
sCCC<-diag(sqrt(ifelse(sCCC<0,0,sCCC)))

mu<-fitzsolver2.cpp(S=fitzhugh$S,THETA2=theta,qrbeval=qrbeval,QRB=fitzhugh$QRB,SAPV=fitzhugh$SAPV,randoms1 = matrix(rnorm(B*length(fitzhugh$SAPV)),ncol=B),randoms2 = matrix(rnorm(B*length(fitzhugh$SAPV)),ncol=B))
mu<-mu+matrix(as.vector(sCCC%*%matrix(rnorm(n=length(grid)*B),nrow=length(grid))),byrow=TRUE,nrow=B)

mu}	
	
FitzhughLimits<-function(i, j, d){
x<-c(-1,sort(runif(9998))*2-1,1)
n<-dim(d)[1]
for(t in (1:n)[-i]){
int<-d[t,]+c(-1,1)*(0.25/10)
x<-x[x<int[1] | x>int[2]]}
x}

FitzhughTerminal<-function(rep, utility = "SIG", scaled = TRUE){

if(rep<0|rep>20){
stop("rep should be 1, 2, ..., 20")}  
	
if(utility=="SIG"){
des<-as.vector(fitzopt$sig.designs[fitzopt$sig.designs[,1]==rep,-1])}
if(utility=="NSEL"){
des<-as.vector(fitzopt$nsel.designs[fitzopt$nsel.designs[,1]==rep,-1])}
if(utility=="NAEL"){
des<-as.vector(fitzopt$nael.designs[fitzopt$nael.designs[,1]==rep,-1])}

if(!scaled){
des<-10*(des+1)} 

sort(des)}  
	
FitzhughEvals<-function(rep = NULL, utility = "SIG"){

if(!is.null(rep)){	
if(rep<0|rep>20){
stop("rep should be 1, 2, ..., 20")}} 
	
if(utility=="SIG"){
out<-fitzopt$sig.out}
if(utility=="NSEL"){
out<-fitzopt$nsel.out}
if(utility=="NAEL"){
out<-fitzopt$nael.out}

if(!is.null(rep)){
out<-out[rep,]}
	
out}

FitzhughDefaultEvals<-function(utility = "SIG"){

if(utility=="SIG"){
out<-fitzopt$sig.def}
if(utility=="NSEL"){
out<-fitzopt$nsel.def}
if(utility=="NAEL"){
out<-fitzopt$nael.def}

out}

FitzhughFinal<-function(utility = "SIG", scaled = TRUE){

if(utility=="SIG"){
des<-as.vector(fitzopt$sig.design)}
if(utility=="NSEL"){
des<-as.vector(fitzopt$nsel.design)}
if(utility=="NAEL"){
des<-as.vector(fitzopt$nael.design)}

if(!scaled){
des<-10*(des+1)} 

sort(des)}    
  

##################################################################################################

### COMPARTMENTAL#################################################################################

##################################################################################################

compsolver2.cpp <- function(S, THETA2, qrbeval, QRB, SAPV, randoms1, randoms2){
	.Call( "compsolver2cpp", S, THETA2, qrbeval, QRB, SAPV, randoms1, randoms2, PACKAGE = "aceodes" )
}

compartmentalnsel.cpp <- function(THETA, mu, v, y){
	.Call( "compartmentalnsel4cpp", THETA, mu, v, y, PACKAGE = "aceodes" )
}

compartmentalsig.cpp <- function(mu, v, y){
	.Call( "compartmentalsigcpp", mu, v, y, PACKAGE = "aceodes" )
}

compartmentalnael.cpp <- function(THETA, mu, v, y, order){
	.Call( "compartmentalnaelcpp", THETA, mu, v, y, order, PACKAGE = "aceodes" )
}

compartmentalrr<-function(u1,v1){
N1<-length(u1)
N2<-length(v1)
u<-matrix(rep(u1,N2),ncol=N2,nrow=N1,byrow=FALSE)
v<-matrix(rep(v1,N1),ncol=N2,nrow=N1,byrow=TRUE)
exp(-0.25*((u-v)^2)/compartmental$lam2)*spi*compartmental$lam}

compartmentalrq<-function(u1,v1){
N1<-length(u1)
N2<-length(v1)
u<-matrix(rep(u1,N2),ncol=N2,nrow=N1,byrow=FALSE)
v<-matrix(rep(v1,N1),ncol=N2,nrow=N1,byrow=TRUE)
pi*compartmental$lam2*(erf(0.5*(v-u)/compartmental$lam)+erf(0.5*(u-compartmental$a)/compartmental$lam))}

compartmentalqr<-function(u1,v1){
t(compartmentalrq(v1,u1))}

compartmentalqq<-function(u1,v1){
N1<-length(u1)
N2<-length(v1)
u<-matrix(rep(u1,N2),ncol=N2,nrow=N1,byrow=FALSE)
v<-matrix(rep(v1,N1),ncol=N2,nrow=N1,byrow=TRUE)
pi*compartmental$lam2*(u-compartmental$a)*erf(0.5*(u-compartmental$a)/compartmental$lam)+2*spi*compartmental$lam3*exp(-0.25*((u-compartmental$a)^2)/compartmental$lam2)-pi*compartmental$lam2*(v-u)*erf(0.5*(v-u)/compartmental$lam)-2*spi*compartmental$lam3*exp(-0.25*((v-u)^2)/compartmental$lam2)+pi*compartmental$lam2*(v-compartmental$a)*erf(0.5*(v-compartmental$a)/compartmental$lam)+2*spi*compartmental$lam3*exp(-0.25*((v-compartmental$a)^2)/compartmental$lam2)-2*spi*compartmental$lam3}

CompartmentalNaelProb<-function(d, B){

n0<-dim(d)[1]
D<-400  
BB<-B+B

THETA2<-cbind(rlnorm(n=BB,meanlog=log(0.1),sdlog=sqrt(0.05)),rlnorm(n=BB,meanlog=log(1),sdlog=sqrt(0.05)),rlnorm(n=BB,meanlog=log(20),sdlog=sqrt(0.05)))

grid<-sort((as.vector(d)+1)*0.5*24)

qrbeval<-compartmentalqr(grid,compartmental$S[1:compartmental$N0])%*%(compartmental$Binv/compartmental$alpha)
qqqq<-compartmentalqq(grid,grid)/compartmental$alpha
rqrq<-compartmentalrq(compartmental$S[1:compartmental$N0],grid)/compartmental$alpha
sCCC<-diag(qqqq-qrbeval%*%rqrq)
sCCC<-ifelse(is.na(sCCC),0,sCCC)
sCCC<-diag(sqrt(ifelse(sCCC<0,0,sCCC)))

mu<-compsolver2.cpp(S=compartmental$S,THETA2=THETA2,qrbeval=qrbeval,QRB=compartmental$QRB,SAPV=compartmental$SAPV, randoms1 = matrix(rnorm(BB*length(compartmental$SAPV)),ncol=BB),randoms2 = matrix(rnorm(BB*length(compartmental$SAPV)),ncol=BB))
mu<-mu+matrix(as.vector(sCCC%*%matrix(rnorm(n=n0*BB),nrow=n0)),byrow=TRUE,nrow=BB)
v<-0.1+0.01*mu^2

y<-mu[1:B,]+sqrt(v[1:B,])*matrix(rnorm(n0*B),ncol=n0)

mu<-matrix(mu[-(1:B),],nrow=B)
v<-matrix(v[-(1:B),],nrow=B)

out<-compartmentalnael.cpp(THETA=THETA2[-(1:B),],mu=mu,v=v,y=y,order=apply(THETA2[-(1:B),],2,order)) 

out2<-apply(abs(THETA2[1:B,]-out),1,sum)

out2[is.na(out2)]<-1.10 # Set any NAs to approximate performance of uniform design

-out2}

CompartmentalNaelExact<-function(d, B){
  
D<-400  
n0<-dim(d)[1]
  BB<-B+B
  
  THETA2<-cbind(rlnorm(n=BB,meanlog=log(0.1),sdlog=sqrt(0.05)),rlnorm(n=BB,meanlog=log(1),sdlog=sqrt(0.05)),rlnorm(n=BB,meanlog=log(20),sdlog=sqrt(0.05)))
  
  grid<-sort((as.vector(d)+1)*0.5*24)
  
  mu<-(D*matrix(rep(THETA2[,2],n0),ncol=n0,byrow=FALSE)/matrix(rep(THETA2[,3]*(THETA2[,2]-THETA2[,1]),n0),ncol=n0,byrow=FALSE))*(exp(-matrix(rep(THETA2[,1],n0),ncol=n0,byrow=FALSE)*matrix(rep(grid,BB),ncol=n0,byrow=TRUE))-exp(-matrix(rep(THETA2[,2],n0),ncol=n0,byrow=FALSE)*matrix(rep(grid,BB),ncol=n0,byrow=TRUE)))
  v<-0.1+0.01*mu^2
  
  y<-mu[1:B,]+sqrt(v[1:B,])*matrix(rnorm(n0*B),ncol=n0)
  
  mu<-matrix(mu[-(1:B),],nrow=B)
  v<-matrix(v[-(1:B),],nrow=B)
  
  out<-compartmentalnael.cpp(THETA=THETA2[-(1:B),],mu=mu,v=v,y=y,order=apply(THETA2[-(1:B),],2,order)) 
  
  out2<-apply(abs(THETA2[1:B,]-out),1,sum)
  out2[is.na(out2)]<-1.10 # Set any NAs to approximate performance of uniform design

  -out2}

CompartmentalNselProb<-function(d, B){

n0<-dim(d)[1]  
D<-400
BB<-B+B

THETA2<-cbind(rlnorm(n=BB,meanlog=log(0.1),sdlog=sqrt(0.05)),rlnorm(n=BB,meanlog=log(1),sdlog=sqrt(0.05)),rlnorm(n=BB,meanlog=log(20),sdlog=sqrt(0.05)))

grid<-sort((as.vector(d)+1)*0.5*24)

qrbeval<-compartmentalqr(grid,compartmental$S[1:compartmental$N0])%*%(compartmental$Binv/compartmental$alpha)
qqqq<-compartmentalqq(grid,grid)/compartmental$alpha
rqrq<-compartmentalrq(compartmental$S[1:compartmental$N0],grid)/compartmental$alpha
sCCC<-diag(qqqq-qrbeval%*%rqrq)
sCCC<-ifelse(is.na(sCCC),0,sCCC)
sCCC<-diag(sqrt(ifelse(sCCC<0,0,sCCC)))

mu<-compsolver2.cpp(S=compartmental$S,THETA2=THETA2,qrbeval=qrbeval,QRB=compartmental$QRB,SAPV=compartmental$SAPV, randoms1 = matrix(rnorm(BB*length(compartmental$SAPV)),ncol=BB),randoms2 = matrix(rnorm(BB*length(compartmental$SAPV)),ncol=BB))
mu<-mu+matrix(as.vector(sCCC%*%matrix(rnorm(n=n0*BB),nrow=n0)),byrow=TRUE,nrow=BB)
v<-0.1+0.01*mu^2

y<-mu[1:B,]+sqrt(v[1:B,])*matrix(rnorm(n0*B),ncol=n0)

mu<-matrix(mu[-(1:B),],nrow=B)
v<-matrix(v[-(1:B),],nrow=B)

out<-compartmentalnsel.cpp(THETA=THETA2[-(1:B),],mu=mu,v=v,y=y) 

out2<-apply((THETA2[1:B,]-out)^2,1,sum)
out2[is.na(out2)]<-1.6 # Set any NAs to approximate performance of uniform design

-out2}

CompartmentalNselExact<-function(d, B){
  
D<-400
  n0<-dim(d)[1]
  BB<-B+B
  
  THETA2<-cbind(rlnorm(n=BB,meanlog=log(0.1),sdlog=sqrt(0.05)),rlnorm(n=BB,meanlog=log(1),sdlog=sqrt(0.05)),rlnorm(n=BB,meanlog=log(20),sdlog=sqrt(0.05)))
  
  grid<-sort((as.vector(d)+1)*0.5*24)
  
  mu<-(D*matrix(rep(THETA2[,2],n0),ncol=n0,byrow=FALSE)/matrix(rep(THETA2[,3]*(THETA2[,2]-THETA2[,1]),n0),ncol=n0,byrow=FALSE))*(exp(-matrix(rep(THETA2[,1],n0),ncol=n0,byrow=FALSE)*matrix(rep(grid,BB),ncol=n0,byrow=TRUE))-exp(-matrix(rep(THETA2[,2],n0),ncol=n0,byrow=FALSE)*matrix(rep(grid,BB),ncol=n0,byrow=TRUE)))
  v<-0.1+0.01*mu^2
  
  y<-mu[1:B,]+sqrt(v[1:B,])*matrix(rnorm(n0*B),ncol=n0)
  
  mu<-matrix(mu[-(1:B),],nrow=B)
  v<-matrix(v[-(1:B),],nrow=B)
  
  out<-compartmentalnsel.cpp(THETA=THETA2[-(1:B),],mu=mu,v=v,y=y) 
  
  out2<-apply((THETA2[1:B,]-out)^2,1,sum)
 	out2[is.na(out2)]<-1.6 # Set any NAs to approximate performance of uniform design
 
  -out2}

CompartmentalSigProb<-function(d, B){

D<-400
n0<-dim(d)[1]
BB<-B+B

THETA2<-cbind(rlnorm(n=BB,meanlog=log(0.1),sdlog=sqrt(0.05)),rlnorm(n=BB,meanlog=log(1),sdlog=sqrt(0.05)),rlnorm(n=BB,meanlog=log(20),sdlog=sqrt(0.05)))

grid<-sort((as.vector(d)+1)*0.5*24)

qrbeval<-compartmentalqr(grid,compartmental$S[1:compartmental$N0])%*%(compartmental$Binv/compartmental$alpha)
qqqq<-compartmentalqq(grid,grid)/compartmental$alpha
rqrq<-compartmentalrq(compartmental$S[1:compartmental$N0],grid)/compartmental$alpha
sCCC<-diag(qqqq-qrbeval%*%rqrq)
sCCC<-ifelse(is.na(sCCC),0,sCCC)
sCCC<-diag(sqrt(ifelse(sCCC<0,0,sCCC)))

mu<-compsolver2.cpp(S=compartmental$S,THETA2=THETA2,qrbeval=qrbeval,QRB=compartmental$QRB,SAPV=compartmental$SAPV, randoms1 = matrix(rnorm(BB*length(compartmental$SAPV)),ncol=BB),randoms2 = matrix(rnorm(BB*length(compartmental$SAPV)),ncol=BB))
mu<-mu+matrix(as.vector(sCCC%*%matrix(rnorm(n=n0*BB),nrow=n0)),byrow=TRUE,nrow=BB)
v<-0.1+0.01*mu^2

y<-mu[1:B,]+sqrt(v[1:B,])*matrix(rnorm(n0*B),ncol=n0)

out1<-apply(matrix(dnorm(x=as.vector(y),mean=as.vector(mu[1:B,]),sd=sqrt(as.vector(v[1:B,])),log=TRUE),ncol=n0),1,sum)+0.5*n0*log(2*pi)

mu<-mu[-(1:B),]
v<-v[-(1:B),]

out<-compartmentalsig.cpp(mu=mu,v=v,y=y) 

out2<-out1-as.vector(out)
out2[is.na(out2)]<-3.7 # Set any NAs to approximate performance of uniform design

out2}

CompartmentalSigExact<-function(d, B){
  
D<-400
  n0<-dim(d)[1]
  BB<-B+B
  
  THETA2<-cbind(rlnorm(n=BB,meanlog=log(0.1),sdlog=sqrt(0.05)),rlnorm(n=BB,meanlog=log(1),sdlog=sqrt(0.05)),rlnorm(n=BB,meanlog=log(20),sdlog=sqrt(0.05)))
  
  grid<-sort((as.vector(d)+1)*0.5*24)
  
  mu<-(D*matrix(rep(THETA2[,2],n0),ncol=n0,byrow=FALSE)/matrix(rep(THETA2[,3]*(THETA2[,2]-THETA2[,1]),n0),ncol=n0,byrow=FALSE))*(exp(-matrix(rep(THETA2[,1],n0),ncol=n0,byrow=FALSE)*matrix(rep(grid,BB),ncol=n0,byrow=TRUE))-exp(-matrix(rep(THETA2[,2],n0),ncol=n0,byrow=FALSE)*matrix(rep(grid,BB),ncol=n0,byrow=TRUE)))
  v<-0.1+0.01*mu^2
  
  y<-mu[1:B,]+sqrt(v[1:B,])*matrix(rnorm(n0*B),ncol=n0)
  
  out1<-apply(matrix(dnorm(x=as.vector(y),mean=as.vector(mu[1:B,]),sd=sqrt(as.vector(v[1:B,])),log=TRUE),ncol=n0),1,sum)+0.5*n0*log(2*pi)
  
  mu<-mu[-(1:B),]
  v<-v[-(1:B),]
  
  out<-compartmentalsig.cpp(mu=mu,v=v,y=y) 
  
  out2<-out1-as.vector(out)
  out2[is.na(out2)]<-3.7 # Set any NAs to approximate performance of uniform design
  out2}

CompartmentalLimits<-function(i, j, d){
x<-c(-1,sort(runif(9998))*2-1,1)
n<-dim(d)[1]
for(t in (1:n)[-i]){
int<-d[t,]+c(-1,1)*(0.25/12)
x<-x[x<int[1] | x>int[2]]}
x}

CompartmentalTerminal<-function(rep, utility = "SIG", scaled = TRUE, solver = "probabilistic"){

if(rep<0|rep>20){
stop("rep should be 1, 2, ..., 20")}  
	
if(solver=="probabilistic" & utility=="SIG"){
des<-as.vector(compopt$sig.designs[compopt$sig.designs[,1]==rep,-1])}
if(solver=="probabilistic" & utility=="NSEL"){
des<-as.vector(compopt$nsel.designs[compopt$nsel.designs[,1]==rep,-1])}
if(solver=="probabilistic" & utility=="NAEL"){
des<-as.vector(compopt$nael.designs[compopt$nael.designs[,1]==rep,-1])}

if(solver=="exact" & utility=="SIG"){
des<-as.vector(compopt$sig.designs.exact[compopt$sig.designs.exact[,1]==rep,-1])}
if(solver=="exact" & utility=="NSEL"){
des<-as.vector(compopt$nsel.designs.exact[compopt$nsel.designs.exact[,1]==rep,-1])}
if(solver=="exact" & utility=="NAEL"){
des<-as.vector(compopt$nael.designs.exact[compopt$nael.designs.exact[,1]==rep,-1])}

if(!scaled){
des<-12*(des+1)} 

sort(des)}  	
	
CompartmentalEvals<-function(rep = NULL, utility = "SIG", solver = "probabilistic"){

if(!is.null(rep)){	
if(rep<0|rep>20){
stop("rep should be 1, 2, ..., 20")}} 
	
if(solver=="probabilistic" & utility=="SIG"){
out<-compopt$sig.out}
if(solver=="probabilistic" & utility=="NSEL"){
out<-compopt$nsel.out}
if(solver=="probabilistic" & utility=="NAEL"){
out<-compopt$nael.out}

if(solver=="exact" & utility=="SIG"){
out<-compopt$sig.exact.out}
if(solver=="exact" & utility=="NSEL"){
out<-compopt$nsel.exact.out}
if(solver=="exact" & utility=="NAEL"){
out<-compopt$nael.exact.out}

if(!is.null(rep)){
out<-out[rep,]}
	
out}
	
CompartmentalFinal<-function(utility = "SIG", scaled = TRUE, solver = "probabilistic"){

if(solver=="probabilistic" & utility=="SIG"){
des<-as.vector(compopt$sig.design)}
if(solver=="probabilistic" & utility=="NSEL"){
des<-as.vector(compopt$nsel.design)}
if(solver=="probabilistic" & utility=="NAEL"){
des<-as.vector(compopt$nael.design)}

if(solver=="exact" & utility=="SIG"){
des<-as.vector(compopt$sig.design.exact)}
if(solver=="exact" & utility=="NSEL"){
des<-as.vector(compopt$nsel.design.exact)}
if(solver=="exact" & utility=="NAEL"){
des<-as.vector(compopt$nael.design.exact)}


if(!scaled){
des<-12*(des+1)} 

sort(des)}   
