#install.packages("Matrix") #Install package if needed
library(Matrix)

#Create data
n=240
x=seq(-1,1,len=n/2)
f=c(exp(-abs(x)/0.15),rep(0,n/6),rep(1,n/6),rep(0,n/6))
xkern=seq(-1,1,len=n)
dx=xkern[2]-xkern[1]
xkern=xkern-dx/2
kern=dnorm(xkern,0,0.05)*dx
kern=c(kern[(n/2+1):n],kern[1:(n/2)])
A=toeplitz(kern)
Af=as.numeric(A%*%f)
noise=0.02*diff(range(Af))
y=Af+noise*rnorm(n)


source("C:/Users/Arttu/OneDrive - Oulun yliopisto/Radar deconvolution/codes_to_github/functions.R") #Change path to location of R codes
init=c(rep(log(n/5),n))			#Initial guess of parameters
lowerbounds=10e-6				
upperbounds=10*n				#Bounds of parameters to not run to numerical errors
alpha=1.5					#Regularisation parameter
prior=1 					#Prior=1 is Cauchy and prior=2 is TV



##########################################################################################
#Use built-in function optim() to maximise marg_post()
pars=optim(init,marg_post,method="L-BFGS-B",lower=lowerbounds,upper=upperbounds,
	control=list(maxit=1e7),n=n,data=y,A=A,noise=noise,prior=prior,gr=grad,alpha=alpha)


#Visualise results
u=pars$par[1:n]
m=2*n
l=exp(u)
l=c(rep(l[1],n/2),l,rep(l[n],n/2))	
a=(-l^2)/sqrt(4*l)
b=(1+2*l^2)/sqrt(4*l)
L=bandSparse(m,m,c(-m+1,-1,0,1,m-1),list(a[m],a[2:m],b,a[1:(m-1)],a[1]))
C=solve(crossprod(L))
P=cbind(matrix(0,n,n/2),diag(n),matrix(0,n,n/2))
P=as(P,"dgCMatrix")
K=A%*%P%*%C%*%t(P)%*%t(A)+noise^2*Diagonal(n)

fhat=(C%*%t(P)%*%t(A)%*%solve(K,y))[(n/2+1):(n/2+n)]
par(mfrow=c(2,1))
plot(f,type="l",ylim=c(-0.2,1.2),main="Ground truth (black) and estimate (red)")
lines(fhat,col="red",lwd=2)
plot(u,type="l",main="Log length-scale")