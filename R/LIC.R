#' The LIC criterion is to determine the most informative subsets so that the subset can retain most of the information contained in the complete data.
#'
#' @param X is a design matrix
#' @param Y is a random response vector of observed values
#' @param alpha is the significance level
#' @param K is the number of subsets
#' @param nk is the sample size of subsets
#'
#' @return MUopt,Bopt,MAEMUopt,MSEMUopt,opt,Yopt
#' @export
#'

#' @examples
#' set.seed(12)
#' X=matrix(data=sample(1:3,1200*5, replace = TRUE) ,nrow=1200,ncol=5)  
#' b=sample(1:3,5, replace = TRUE)         
#' e= rnorm(1200, 0, 1)    
#' Y=X%*%b+e
#' alpha=0.05	
#' K=10
#' nk=1200/K 
#' LIC(X,Y,alpha,K,nk)

LIC=function(X,Y,alpha,K,nk){
n=nrow(X)
p=ncol(X)
nk=n/K
N=L1=c(1:K)  
Rm=matrix(rep(0, nk*K),ncol=K) 
mr=matrix(rep(0,K*nk), ncol=nk)
for(i in 1:K )  {
mr[i, ]=sample(1:n,nk,replace=T)  
r=matrix(c(1:nk,mr[i, ]),ncol=nk,byrow=T)	
Rm[,i]=r[2,]
R=matrix(rep(0, nk*n),ncol=n)
R[t(r)]=1 
X1=R%*%X 	
Y1=R%*%Y 
Hr=X1%*%solve(crossprod(X1))%*%t(X1) 
I1=diag(rep(1,nk)) 
SX=(t(Y1)%*%(I1-Hr)%*%Y1)/(nk-p)  
SY=sqrt(t(Y1)%*%(I1-Hr)%*%Y1)/(nk-p) 
C1=sum(diag(X1%*%solve(crossprod(X1))%*%t(X1)))/nk   
L1[i]=2*SY*C1*qt(1-alpha/2,nk-p)  
N[i]=det(t(X1)%*%X1)  
}  
opt1=Rm[,which.min(L1)]
opt2=Rm[,which.max(N)]
opt=intersect(opt1,opt2)
Yopt=Y[opt]
Xopt=X[opt,]
Bopt=solve(crossprod(Xopt))%*%t(Xopt)%*%Yopt   
MUopt=Xopt%*%Bopt
Nopt=length(Yopt) 
E5=(t(Yopt-MUopt)%*%(Yopt-MUopt))/Nopt    
A5=sum(abs(Yopt-MUopt))/Nopt
return(list(MUopt=MUopt,Bopt=Bopt,MAEMUopt=A5,MSEMUopt=E5,opt=opt,Yopt=Yopt))
}
