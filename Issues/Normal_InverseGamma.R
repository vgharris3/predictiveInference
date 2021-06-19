set.seed(1234)

if(0){

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
myyt = location + scale * stats::rt(S,n-1)

plot(density(y.postsample))
#lines(density(y.J),lty=2,lwd=2,col="red")
lines(density(yt),col="blue",lty=3,lwd=2)
points(density(myyt),pch=5,col="green")

legend("topright",legend=c("postsample","Jeffrey's posteriors","Jeffrey's t","my Jeffrey's t"),lty=c(1,2,3,NA),col=c("black","red","blue","green"),pch=c(NA,NA,NA,5))

}

leg = rep(0,length(mu0))

for(i in 1:length(mu0)){
  c(i,mu0[i])
  leg[i] = paste("mu0 =",mu0[i])
}

###



<<fig=TRUE, echo=FALSE>>=
  if(0){
    source('C:/Users/gabe/Documents/tmp/predictiveInference/R/rpredNormIG.R')
    source('C:/Users/gabe/Documents/tmp/predictiveInference/R/dpredNormIG.R')
    source('C:/Users/gabe/Documents/tmp/predictiveInference/R/ppredNormIG.R')
    set.seed(1234)

    if(0){

      #WITH PRIOR KNOWLEDGE--BETTER PREDICTION
      #prior
      mu0 = 1.9; k0 = 1
      mu0 = 1.9;  #PLAY AROUND WITH THIS.
      mu0 = seq(1.5,2.5,0.1)
      sig20 = 0.010; nu0 = 1
      #sig20 = 10;

      y = c(1.64,1.7,1.72,1.74,1.82,1.82,1.82,1.90,2.08)
      n = length(y);
      ybar = mean(y);
      s2 = var(y);

      kn = k0 + n; nun = nu0 + n
      mun = (k0*mu0+n*ybar)/kn
      sig2n = (nu0*sig20 + (n-1)*s2 + k0*n*(ybar-mu0)^2/kn)/nun

      ####
      ####RANDOM SAMPLE:  rpredNormIG()
      ####

      S = 10000

      #Monte Carlo for marginal posterior samples of sig2 and theta

      set.seed(1234)
      random_sample_MC = matrix(0,S,length(mu0))
      densities_x = densities_y = matrix(0,512,length(mu0))
      for(i in 1:length(mu0)){
        random_sample_MC[,i] = rpredNormIG(S,y,mu0[i],k0,sig20,nu0)
        densities_x[,i] = density(random_sample_MC[,i])$x
        densities_y[,i] = density(random_sample_MC[,i])$y

      }
      #random_sample_MC = rpredNormIG(S,y,mu0,k0,sig20,nu0)

    }

    #Use Jeffrey's prior

    random_sample_J = rpredNormIG(S,y,Jeffreys=TRUE)

    #Combine with the mu0=1.9 plot on one plot for comparison

    plotTitle = "Normal-Inverse Gamma Density"
    matplot(densities_x[,pos1.9],densities_y[,pos1.9],
            type="l",lty=1,lwd=2,col="blue",
            xlim = c(1,3),
            xlab = "y", ylab = "Prediction Density", main=plotTitle)
    lines(density(random_sample_J),lty=2,lwd=2,col="red")

    legend("topright", inset = 0.05, legend = c(paste("Using mu0 =",mu0[pos1.9]),"Using Jeffreys prior"),
           lty=c(1,2), lwd=2,
           cex = 0.8, col=c("blue","red"), box.col="white",
    )

  }

if(0){

  par(mfrow=c(1,1))
  #H_MC = hist(random_sample_MC,probability = TRUE,xlim = c(1,2.5),ylim = c(0,3),breaks=50) plot=FALSE)
  #H_J = hist(random_sample_J,probability=TRUE,xlim = c(1,2.5),ylim=c(0,3),breaks=50) plot=FALSE)
  H_MC = hist(random_sample_MC,probability = TRUE,ylim = c(0,3),breaks=50, plot=FALSE)
  H_J = hist(random_sample_J,probability=TRUE,ylim=c(0,3),breaks=50, plot=FALSE)

  #COMBINE BOTH MIDS THEN TAKE SORT(UNIQUE(MIDS))
  mids = sort(unique(c(H_MC$mids,H_J$mids)))
  #if(length(H_MC$mids) > length(H_J$mids)){
  #  mids = H_MC$mids
  #}else{
  #  mids = H_J$mids
  #}

  counts_MC = numeric(length(mids))
  counts_J = numeric(length(mids))

  for(hmcInd in 1:length(H_MC$mids)){
    counts_MC[which(mids==H_MC$mids[hmcInd])] = H_MC$counts[hmcInd];
  }

  for(hjInd in 1:length(H_J$mids)){
    counts_J[which(mids==H_J$mids[hjInd])] = H_J$counts[hjInd];
  }

  pct_MC = counts_MC/sum(counts_MC)
  pct_J = counts_J/sum(counts_J)

  plotdata = rbind(pct_MC,pct_J)
  plotdata = rbind(counts_MC,counts_J)

  barplot(plotdata,col=c("blue","red"),beside=TRUE,names.arg = mids, main="Counts")
  legend("topright",legend=c("MC with previous \nknowledge","t-distribution using \nJeffreys prior"),fill=c("blue","red"),y.intersp=2,box.col="white")

  #MAYBE PLOT AS LINES AND SEE HOW DIFFERENT MU0 AFFECT THE CURVE

  if(1) {
    xn = 1000

    x = sort(runif(xn,min = min(random_sample),max = max(random_sample)))
    tictoc::tic()
    xd = dpredNormIG(x,y,mu0,k0,sig20,nu0)
    tictoc::toc()

    plot(x,xd,type="l")

    #Jeffrey's
    #WITHOUT PRIOR KNOWLEDGE, JUST USING OBSERVATIONS.  STILL PRETTY GOOD PREDICTIONS

    xd_J = dpredNormIG(x,y,Jeffreys=TRUE)
    lines(x,xd_J,lty=2,col="red")


    #s2 = var(y)
    #ybar = mean(y)
    #n = length(y)
    #location = ybar
    #scale = sqrt(s2)*sqrt(1+1/n)
    #xt = seq(1,2.5,len=100)
    #yt = (1/scale) * dt((xt - location)/scale,df = n-1)
    #points(xt,yt)
    #lines(xt,yt)

    #legend("topright",legend=c("dpredNormIG()","Jeffrey's"),lty=1,pch=c(NA,1),col=c("red","black"))

    #ppredNormIG()

    par(mfrow=c(1,1))

    xp = ppredNormIG(x,y,mu0,k0,sig20,nu0,Jeffreys=FALSE)
    plot(x,xp,type="l")

    xp_J = ppredNormIG(x,y,Jeffreys=TRUE)
    points(x,xp_J,col="red")

    #xp_J2 = stats::pt((x-location)/scale,df=n-1)
    #points(x,xp_J,pch=20,col="blue")
    #x = seq(0,3,0.1)
    #tictoc::tic()
    #xp = ppredNormIG(x,y,mu0,k0,sig20,nu0)
    #tictoc::toc()
    #plot(x,xp)
    #lines(ecdf(random_sample))
  }


}


@

