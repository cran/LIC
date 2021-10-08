#' The OSM is a median processing method for the central processor.
#'
#' @param X is a design matrix
#' @param Y is a random response vector of observed values
#' @param alpha is the significance level
#' @param K is the number of subsets
#' @param nk is the sample size of subsets
#'
#' @return MUM,BetaM,MAEMUM,MSEMUM
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
#' OSM(X,Y,alpha,K,nk)

OSM=function(X,Y,alpha,K,nk){
n=nrow(X)
p=ncol(X) 
nk=n/K         
BB= matrix(rep(0,p*K), ncol=K)  
MU1=Ym=matrix(rep(0, nk*K),ncol=K)    
Xm=matrix(rep(0, (nk*p)*K),ncol=K)     
Rm=matrix(rep(0, nk*K),ncol=K)     
mr=matrix(rep(0,K*nk), ncol=nk)	    
for(i in 1:K ) {
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
MUM=matrix(rep(0, nk*1),ncol=1) 
for(i in 1:nk)  {      
MUM[i,]=median(MU1[i,])
}                    
Bm= matrix(rep(0, p*1), ncol=1)
for(i in 1:p )  {
Bm[i,]=median(BB[i,])
}                          
EM=M2=c(1:K)
for(i in 1:K )  {
XIK= matrix(Xm[,i],nrow = nk,ncol = p)  
YIK=Ym[,i]
EM[i]= (t(YIK-XIK%*%Bm)%*%(YIK-XIK%*%Bm))/nk 
M2[i]= sum(abs(YIK-XIK%*%Bm))/nk
} 
E2=min(EM)
A2=min(M2) 
return(list(MUM=MUM,BetaM=Bm,MAEMUM=A2, MSEMUM=E2))
}
