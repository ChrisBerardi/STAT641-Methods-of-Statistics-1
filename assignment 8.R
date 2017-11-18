#P8
site1=c(0,0,0,4,2,2,5,4,2,1,0,14,1,30,0,3,28,2,21,8,82,12,10,2,0)
site2=c(0,1,1,2,2,1,2,29,2,2,0,13,0,19,1,3,26,30,5,4,94,1,9,3,0)
site3=c(0,0,0,2,3,0,0,4,0,5,4,22,0,64,4,4,43,3,16,19,95,6,22,0,0)

#A none are normal
#site1
shapiro.test(site1)
y = sort(site1)
n = length(site1)
i = seq(1,n,1)
u = (i-.375)/(n+.25)
q = qnorm(u)
qqnorm(site1,main="Normal Prob Plots of the hermit crab site 1", 
xlab="normal quantiles",ylab="Number of Hermit Crabs",cex=.65)
qqline(site1)

#site2
shapiro.test(site2)
y = sort(site2)
n = length(site2)
i = seq(1,n,1)
u = (i-.375)/(n+.25)
q = qnorm(u)
qqnorm(site2,main="Normal Prob Plots of the hermit crab site 1", 
xlab="normal quantiles",ylab="Number of Hermit Crabs",cex=.65)
qqline(site2)

#site3
shapiro.test(site3)
y = sort(site3)
n = length(site3)
i = seq(1,n,1)
u = (i-.375)/(n+.25)
q = qnorm(u)
qqnorm(site3,main="Normal Prob Plots of the hermit crab site 1", 
xlab="normal quantiles",ylab="Number of Hermit Crabs",cex=.65)
qqline(site3)

#B use BFL test since data is non-normal, not the same sigma
library(car)
site=c(site1,site2,site3)
n=length(site1)
grp=c(rep(1,n),rep(2,n),rep(3,n))
grp=factor(grp)

leveneTest(site,grp)

#C non-normal, so no von Neumann, none are correlated
#site1
means = mean(site1)
diff = site1-means

n.neg = rep(0,n)
n.pos = rep(0,n)
n.neg = length(diff[diff<0])
n.pos = length(diff[diff>0])

numb.runs = 0
for (j in 2:n) {
if (sign(diff[j]) != sign(diff[j-1])) {numb.runs = numb.runs + 1}
}
runs.result = as.data.frame(cbind(numb.runs, n.pos, n.neg))
names(runs.result) = c("No. runs", "N+", "N-")
runs.result

mu = 1 + (2*n.neg*n.pos)/(n.neg + n.pos)
sig2 = (2*n.neg*n.pos*(2*n.neg*n.pos - n.neg - n.pos))/((n.neg + n.pos)^2*(n.neg + n.pos-1))
z = (abs(numb.runs-mu)-.5)/sqrt(sig2)
pvalue = 2*(1-pnorm(abs(z)))
results = cbind(mu,sig2,z,pvalue)
results

#site2
means = mean(site2)
diff = site2-means

n.neg = rep(0,n)
n.pos = rep(0,n)
n.neg = length(diff[diff<0])
n.pos = length(diff[diff>0])

numb.runs = 0
for (j in 2:n) {
if (sign(diff[j]) != sign(diff[j-1])) {numb.runs = numb.runs + 1}
}
runs.result = as.data.frame(cbind(numb.runs, n.pos, n.neg))
names(runs.result) = c("No. runs", "N+", "N-")
runs.result

mu = 1 + (2*n.neg*n.pos)/(n.neg + n.pos)
sig2 = (2*n.neg*n.pos*(2*n.neg*n.pos - n.neg - n.pos))/((n.neg + n.pos)^2*(n.neg + n.pos-1))
z = (abs(numb.runs-mu)-.5)/sqrt(sig2)
pvalue = 2*(1-pnorm(abs(z)))
results = cbind(mu,sig2,z,pvalue)
results

#site3
means = mean(site3)
diff = site3-means

n.neg = rep(0,n)
n.pos = rep(0,n)
n.neg = length(diff[diff<0])
n.pos = length(diff[diff>0])

numb.runs = 0
for (j in 2:n) {
if (sign(diff[j]) != sign(diff[j-1])) {numb.runs = numb.runs + 1}
}
runs.result = as.data.frame(cbind(numb.runs, n.pos, n.neg))
names(runs.result) = c("No. runs", "N+", "N-")
runs.result

mu = 1 + (2*n.neg*n.pos)/(n.neg + n.pos)
sig2 = (2*n.neg*n.pos*(2*n.neg*n.pos - n.neg - n.pos))/((n.neg + n.pos)^2*(n.neg + n.pos-1))
z = (abs(numb.runs-mu)-.5)/sqrt(sig2)
pvalue = 2*(1-pnorm(abs(z)))
results = cbind(mu,sig2,z,pvalue)
results

#D since non-normal, different sigmas, not correlated and testing mu
#Ranked sum for higher t>2
kruskal.test(list(site1,site2,site3))

#9
#A
white = c(19,141)
black = c(17,52)
combine = rbind(white,black)
chisq.test(combine)

#B
white = c(30,184)
black = c(6,106)
combine = rbind(white,black)
chisq.test(combine)

#C
white = c(19,132,0,9)
black = c(11,52,6,97)
combine = rbind(white,black)
fisher.test(combine)

#D


x <- c(2.9, 3.0, 2.5, 2.6, 3.2) # normal subjects
y <- c(3.8, 2.7, 4.0, 2.4)      # with obstructive airway disease
z <- c(2.8, 3.4, 3.7, 2.2, 2.0) # with asbestosis
kruskal.test(list(x, y, z))