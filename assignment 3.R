small = scan("c:/Users/Saistout/Desktop/641/SmallBrainSize.csv")
large = scan("c:/Users/Saistout/Desktop/641/LargeBrainSize.csv")

# III
 
plot(density(small, kernel="g", bw=3),type="l",
xlab="Brainsize",ylab="PDF",
main="Small Litter Data, Gaussian, bw=3",cex=.5)

estim=rep(0,51)
for (i in 1:51)
	estim[i] = exp(-.5*((3-small[i])/3)^2)/153/sqrt(2*pi)
sum(estim)

estim=rep(0,51)
for (i in 1:51)
	estim[i] = exp(-.5*((16-small[i])/3)^2)/153/sqrt(2*pi)
sum(estim)

hist(small,breaks=5, plot=TRUE, prob=T, xlim=c(0,20),
main="Small Litter Size",
xlab="Brainsize", cex=.75)

#4
#pdf
par(mfrow=c(2,3))
plot(density(small, kernel="g", bw=3),type="l",
xlab="Brainsize",ylab="PDF",
main="Small Litter Data, Gaussian, bw=3",cex=.5)
#cdf
deny1 = density(small,kernel='g',bw=3,n=5000,from=0)
xy1 = deny1$x
pdfy1 = deny1$y
cdfy1 = c(rep(0,5000))
areay1 = c(rep(0,5000))
for (i in 1:5000)
{
areay1[i] = abs((pdfy1[i+1]+pdfy1[i])/2)*(xy1[i+1]-xy1[i])
cdfy1[i] = sum(areay1)
}
plot(xy1,cdfy1,type="l",
xlab="Brainsize",ylab="F(x)",
main="Smoothed Sample cdf for Small Litter Data",cex=.95)
#Quantile
probs = seq(0,1,.01)
Q1 = quantile(small, probs)
plot(probs,Q1,type="l",ylab="Q(u) (size)",xlab="u",ylim=c(0,20),
xlim=c(0,1),lab=c(10,11,7))
title("Empirical Quantile of Small Litter Size",cex=.75)


#pdf
plot(density(large, kernel="g", bw=3),type="l",
xlab="Brainsize",ylab="PDF",
main="Large Litter Data, Gaussian, bw=3",cex=.5)
#cdf
deny1 = density(large,kernel='g',bw=3,n=5000,from=0)
xy1 = deny1$x
pdfy1 = deny1$y
cdfy1 = c(rep(0,5000))
areay1 = c(rep(0,5000))
for (i in 1:5000)
{
areay1[i] = abs((pdfy1[i+1]+pdfy1[i])/2)*(xy1[i+1]-xy1[i])
cdfy1[i] = sum(areay1)
}
plot(xy1,cdfy1,type="l",
xlab="Brainsize",ylab="F(x)",
main="Smoothed Sample cdf for Large Litter Data",cex=.95)
#Quantile
probs = seq(0,1,.01)
Q1 = quantile(large, probs)
plot(probs,Q1,type="l",ylab="Q(u) (size)",xlab="u",ylim=c(0,40),
xlim=c(0,1),lab=c(10,11,7))
title("Empirical Quantile of Large Litter Size",cex=.75)
