small = scan("c:/Users/Saistout/Desktop/641/SmallBrainSize.csv")
large = scan("c:/Users/Saistout/Desktop/641/LargeBrainSize.csv")

#1
#a.
mean(large)
mean(large,.1)
#b.
library(MASS)
fitdistr(large,"weibull")
gamma=1.1326673
alpha=11.4982087
#c.
1-pweibull(30,gamma,scale=alpha)
#d.
(mu_hat=alpha*gamma((1+1/gamma)))
(sig_hat=sqrt(alpha^2*(gamma(1+2/gamma)-gamma(1+1/gamma)^2)))
mean(large)
sd(large)
#e.
qweibull(.5,gamma,scale=alpha)
median(large)

qweibull(.75,gamma,scale=alpha)-qweibull(.25,gamma,scale=alpha)
IQR(large)

#2
mean(large)
sd(large)
mean(small)
sd(small)

#3
median(large)
mad(large)
median(small)
mad(small)


#2
library(survival)
T=c(180,632,2240,195,76,70,13,1990,18,700,210,1296,23,
8,852,52,220,63,8,1976,1296,1460,63,1328,365)
ST=c(rep(1,6),0,0,rep(1,5),0,rep(1,5),0,0,rep(1,4))
G=c(rep(1,13),rep(2,12))
out=cbind(T,ST,G)

drug <- survfit(Surv(T, ST) ~ G)
summary(drug)
one <- c(180,632,2240,195,76,70,13,1990,18,700,210,1296,23)
two <- c(8,852,52,220,63,8,1976,1296,1460,63,1328,365)
mean(one)
mean(two)
plot(drug,ylab="Survival Function",xlab="Time to Recovery",
main="Cancer Treatment - Estimated S(t)",lty=1:2 )
legend(1500,.8,c("Group 1","Group 2"),lty=1:2,lwd=2)