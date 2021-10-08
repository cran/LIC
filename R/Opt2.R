#' The Opt2 chooses the optimal index subset based on maximized information sub-matrix.
#'
#' @param X is a design matrix
#' @param Y is a random response vector of observed values
#' @param alpha is the significance level
#' @param K is the number of subsets
#' @param nk is the sample size of subsets
#'
#' @return MUopt2,Bopt2,MAEMUopt2,MSEMUopt2,opt2,Yopt2
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
#' Opt2(X,Y,alpha,K,nk)

Opt2=function(X,Y,alpha,K,nk){
n=nrow(X)
p=ncol(X)
nk=n/K  
N=c(1:K)   
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
N[i]=det(t(X1)%*%X1)  
}  
opt2=Rm[,which.max(N)]
Yopt2=Y[opt2]
Xopt2=X[opt2,]
Bopt2=solve(crossprod(Xopt2))%*%t(Xopt2)%*%Yopt2   
MUopt2=Xopt2%*%Bopt2
Nopt2=length(Yopt2) 
E4=(t(Yopt2-MUopt2)%*%(Yopt2-MUopt2))/Nopt2     
A4=sum(abs(Yopt2-MUopt2))/Nopt2
return(list(MUopt2=MUopt2,Bopt2=Bopt2,MAEMUopt2=A4,MSEMUopt2=E4,opt2=opt2,Yopt2=Yopt2))
}
