#' The OSA gives a simple average estimatoris by averaging all these least squares estimators.
#'
#' @param X is a design matrix
#' @param Y is a random response vector of observed values
#' @param alpha is the significance level
#' @param K is the number of subsets
#' @param nk is the sample size of subsets
#'
#' @return MUA,BetaA,MAEMUA,MSEMUA
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
#' OSA(X,Y,alpha,K,nk)

OSA=function(X,Y,alpha,K,nk){
n=nrow(X)
p=ncol(X)
nk=n/K         
BB=matrix(rep(0,p*K), ncol=K)
MU1=Ym=matrix(rep(0, nk*K),ncol=K)   
Xm=matrix(rep(0, (nk*p)*K),ncol=K)     
Rm=matrix(rep(0, nk*K),ncol=K)     
mr=matrix(rep(0,K*nk), ncol=nk)	    
for(i in 1:K ){
mr[i, ]=sample(1:n,nk,replace=T)        
r=matrix(c(1:nk,mr[i, ]),ncol=nk,byrow=T)	
Rm[,i]=r[2,]
R=matrix(rep(0, nk*n),ncol=n)   
R[t(r)]=1   	
X1=R%*%X 
Y1=R%*%Y       
B1=solve(crossprod(X1))%*%t(X1)%*%Y1 
MU1[,i]=X1%*%B1      
BB[,i]=B1    
Ym[,i]= Y1     
Xm[,i]= as.vector(X1)
}         
MUA=matrix(rep(0, nk*1),ncol=1) 
for(i in 1:nk)  {      
MUA[i,]=mean(MU1[i,])
}                        
Ba= matrix(rowMeans(BB),p*1, ncol=1)                     
EA=M1=c(1:K)
for(i in 1:K )  {
XIK= matrix(Xm[,i],nrow = nk,ncol = p)  
YIK=Ym[,i]
EA[i]=(t(YIK-XIK%*%Ba)%*%(YIK-XIK%*%Ba))/nk   
M1[i]=sum(abs(YIK-XIK%*%Ba))/nk    
} 
E1=min(EA)
A1=min(M1)
return(list(MUA=MUA,BetaA=Ba,MAEMUA=A1,MSEMUA=E1))
}
