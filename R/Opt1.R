#' The Opt1 chooses the optimal index subset based on minimized interval length.
#'
#' @param X is a design matrix
#' @param Y is a random response vector of observed values
#' @param alpha is the significance level
#' @param K is the number of subsets
#' @param nk is the sample size of subsets
#'
#' @return MUopt1,Bopt1,MAEMUopt1,MSEMUopt1,opt1,Yopt1
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
#' Opt1(X,Y,alpha,K,nk)

Opt1=function(X,Y,alpha,K,nk){
n=nrow(X)
p=ncol(X)
nk=n/K  
L1=c(1:K)    
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
}  
opt1=Rm[,which.min(L1)]
Yopt1=Y[opt1]
Xopt1=X[opt1,]
Bopt1=solve(crossprod(Xopt1))%*%t(Xopt1)%*%Yopt1   
MUopt1=Xopt1%*%Bopt1
Nopt1=length(Yopt1) 
E3=(t(Yopt1-MUopt1)%*%(Yopt1-MUopt1))/Nopt1     
A3=sum(abs(Yopt1-MUopt1))/Nopt1
return(list(MUopt1=MUopt1,Bopt1=Bopt1,MAEMUopt1=A3,MSEMUopt1=E3,opt1=opt1,Yopt1=Yopt1))
}
