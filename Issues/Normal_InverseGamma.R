set.seed(1234)
#prior
mu0 = 1.8; k0 = 1
sig20 = 0.010; nu0 = 1

y = c(1.64,1.7,1.72,1.74,1.82,1.82,1.82,1.90,2.08)
n = length(y);
ybar = mean(y);
s2 = var(y);

kn = k0 + n; nun = nu0 + n
mun = (k0*mu0+n*ybar)/kn
sig2n = (nu0*sig20 + (n-1)*s2 + k0*n*(ybar-mu0)^2/kn)/nun

rbind(kn,mun,sqrt(sig2n),nun)

#Monte Carlo for marginal posterior samples of sig2 and theta
S = 10000
sig2.postsample = 1/rgamma(S,nun/2,sig2n*nun/2)
theta.postsample = rnorm(S,mun,sqrt(sig2.postsample/kn))

y.postsample = rnorm(S,theta.postsample,sqrt(sig2.postsample))

#Monte Carlo again using posteriors resulting from Jeffrey's prior

sig2J = 1/rgamma(S,(n-1)/2,(sum((y-ybar)^2))/2)
thetaJ = rnorm(S,ybar,sqrt(sig2J/n))

y.J = rnorm(S,thetaJ,sqrt(sig2J))

points(density(y.J),col="green",pch=20)

H.postsample = hist(y.postsample,probability = TRUE,ylim = c(0,4))
H.J = hist(y.J,probability=TRUE)

plot(density(y.postsample))
lines(density(y.J),lty=2,lwd=2,col="red")

#predictive t-distribution resulting from Jeffrey's
location = ybar
scale = sqrt(s2)*sqrt(1+1/n)
yt = metRology::rt.scaled(S,n-1,location,scale)

plot(density(y.postsample))
lines(density(y.J),lty=2,lwd=2,col="red")
lines(density(yt),col="blue",lty=3,lwd=2)

legend("topright",legend=c("postsample","Jeffrey's posteriors","Jeffrey's t"),lty=c(1,2,3),col=c("black","red","blue"))
