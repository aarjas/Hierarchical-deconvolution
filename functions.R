#Marginal posterior density
marg_post=function(p,n,data,A,noise,prior,alpha)
{
	P=cbind(matrix(0,n,n/2),diag(n),matrix(0,n,n/2))
	P=as(P,"dgCMatrix")	
	l=exp(p[1:n])
	l=c(rep(l[1],n/2),l,rep(l[n],n/2))	
	a=(-l^2)/sqrt(4*l)
	b=(1+2*l^2)/sqrt(4*l)
	L=bandSparse(2*n,2*n,c(-1,0,1),list(a[2:(2*n)],b,a[1:(2*n-1)]))
	ACAT=crossprod(solve(t(L),t(P)%*%t(A)))
	K=ACAT+noise^2*Diagonal(n)
	du=p[c(2:n,1)]-p[1:n]
	if(prior==1) val=0.5*determinant(K)$modulus+0.5*(t(data)%*%solve(K,data))[1,1]+sum(log(du^2+alpha^2))
	if(prior==2) val=0.5*determinant(K)$modulus+0.5*(t(data)%*%solve(K,data))[1,1]+sum(sqrt(du^2+1e-8))/alpha
	if(runif(1)<1e-1) {par(mfrow=c(1,1));plot(p,type="l",main="Log length-scale")}
	return(val)
}

#Gradient of marg_post
grad=function(p,n,data,A,noise,prior,alpha)
{
	P=cbind(matrix(0,n,n/2),diag(n),matrix(0,n,n/2))
	P=as(P,"dgCMatrix")
	m=2*n
	u=p[1:n]
	u=c(rep(u[1],n/2),u,rep(u[n],n/2))
	l=exp(u)
	a=-l^2/sqrt(4*l)
	b=(1+2*l^2)/sqrt(4*l)
	L=bandSparse(m,m,c(-m+1,-1,0,1,m-1),list(a[m],a[2:m],b,a[1:(m-1)],a[1]))
	LTL=crossprod(L)
	K=crossprod(solve(t(L),t(P)%*%t(A)))+noise^2*Diagonal(n)
	LTLinvAT=solve(L,solve(t(L),t(P)%*%t(A)))
	Kinvg=as.numeric(solve(K,data))
	KinvALTLinv=solve(K,t(LTLinvAT))
	LTLinvATKinvy=as.numeric(LTLinvAT%*%Kinvg)
	ttt2v=as.vector(KinvALTLinv)
	tttv=as.vector(t(LTLinvAT))
	dshape=dlogdetK=rep(0,n)
	las=1
	for(j in (n/2+1):(n/2+n))
	{
		indsj=c(j-1,j,j+1)
		inds2=c((n*(j-1)+1-n):(n*(j-1)+1+2*n-1))
		da=-0.75*exp(u[j])^1.5
		db=-0.25*exp(u[j])^(-0.5)+1.5*exp(u[j])^1.5	
		dL=matrix(c(0,0,0,da,db,da,0,0,0),3,3,byrow=TRUE)	
		dLT=t(dL)
		L=matrix(c(b[indsj[1]],a[indsj[1]],0,a[indsj[2]],b[indsj[2]],a[indsj[2]],0,a[indsj[3]],b[indsj[3]]),3,3,byrow=TRUE)
		dLTL=dLT%*%L+t(L)%*%dL		
		dshape[las]=(t(LTLinvATKinvy[indsj])%*%dLTL%*%LTLinvATKinvy[indsj])[1,1]
		dlogdetK[las]=sum(diag(-matrix(ttt2v[inds2],n,3)%*%dLTL%*%matrix(tttv[inds2],3,n,byrow=TRUE)))
		las=las+1
	}	
	inds1=(n/2+1):(n/2+n)
	inds2=inds1-1
	inds3=inds1+1
	du1=u[inds1]-u[inds2]
	du2=u[inds3]-u[inds1]
	if(prior==1) grad=0.5*dlogdetK+0.5*dshape+2*du1/(du1^2+alpha^2)-2*du2/(du2^2+alpha^2)
	if(prior==2) grad=0.5*dlogdetK+0.5*dshape+2/alpha*(du1/sqrt(du1^2+1e-8)-du2/sqrt(du2^2+1e-8))
	return(grad)
}