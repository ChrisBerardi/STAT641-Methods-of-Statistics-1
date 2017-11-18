#II

observed = c(121,110,38,7,4)
(lamda_hat = (1*110+2*38+3*7+4*3+5*1)/280)
i = c(0,1,2,3,4,5)
(p_hat = c(121/280,110/280,38/280,7/280,3/280,1/280))
(model = dpois(i,lamda_hat))
model[6]=1-ppois(5,lamda_hat)
model[5]=model[5]+model[6]
expected=model*280
chi_test=rep(0,4)
for (j in 1:5)
	chi_test[j]=(observed[j]-expected[j])^2/observed[j]
test = sum(chi_test)
1-pchisq(test,3)

#III
M= 500
o = 50
d= numeric(M)
for(i in 1:o)
{ d[i]= rnorm(1,10,1)
}
for(i in 51:M)
{ d[i]= rnorm(1,30,3)
}
y = sort(d)
n = length(d)
i = seq(1,n,1)
u = (i-.375)/(n+.25)
q = qnorm(u)
qqnorm(d,main="Normal Prob Plots of the CPUE data log(catch)", 
xlab="normal quantiles",ylab=expression(log(catch),cex=.65))
qqline(d)

#IV
x=c(10,12,13,16,37,42,43,45,55,63,66,82,99,100,100,101,107,117,
122,135,138,140,142,149,150,151,154,165,170,183,194,214,218,219,224,229,232,
247,268,268,298,299,325,332,379,400,434,464,499,537)
n = length(x)
i = seq(1,n,1)
y = -log(x)
y = sort(y)
# Anderson-Darling: For Weibull Model
library(MASS)
mle <- fitdistr(x,"weibull")
shape = mle$estimate[1]
scale = mle$estimate[2]
a = -log(scale)
b = 1/shape
z = exp(-exp(-(y-a)/b))
A1i = (2*i-1)*log(z)
A2i = (2*n+1-2*i)*log(1-z)
s1 = sum(A1i)
s2 = sum(A2i)
AD = -n-(1/n)*(s1+s2)
ADM = AD*(1+.2/sqrt(n))
AD
ADM
n
n = length(y)
weib= -y
weib= sort(weib)
i= 1:n
ui= (i-.5)/n
QW= log(-log(1-ui))
plot(QW,weib,abline(lm(weib~QW)),
main="Weibull Reference Plot",cex=.75,lab=c(7,11,7),
xlab="Q=ln(-ln(1-ui))",
ylab="y=ln(W(i))")
legend(-4.7,5.0,"y=-5.304193+.7214465Q")
legend(-4.7,4.5,"AD=.2171, p-value>.25")

#5
#1
x=c(.6,.7,1.1,1.3,1.8,2,2.3,2.7,2.9,3.1,3.9,4.3,4.4,4.9,5.2,5.4,6.1,6.8,
7.1,8,9.4,10.3,12.9,15.9,16,22,22.5,23,23.1,23.9,26.5,26.7,28.4,28.5,32.2,40.2,
42.5,47.2,48.3,55.8,57,57.2,64.9,67.6,71.3,79.5,114.5,128.6,239.5)
logx=log(x)

shapiro.test(logx)
y = sort(logx)
n = length(logx)
i = seq(1,n,1)
u = (i-.375)/(n+.25)
q = qnorm(u)
r = cor.test(q,y)

qqnorm(logx,main="Normal Prob Plots of the CPUE data log(catch)", 
xlab="normal quantiles",ylab=expression(log(catch),cex=.65))
qqline(logx)


#2
y=x
n = length(y)
yt0 = log(y)
s = sum(yt0)
varyt0 = var(yt0)
Lt0 = -1*s - .5*n*(log(2*pi*varyt0)+1)
th = 0
Lt = 0
t = -3.01
i = 0
while(t < 3)
{t = t+.001
i = i+1
th[i] = t
yt = (y^t -1)/t
varyt = var(yt)
Lt[i] = (t-1)*s - .5*n*(log(2*pi*varyt)+1)
if(abs(th[i])<1.0e-10)Lt[i]<-Lt0
if(abs(th[i])<1.0e-10)th[i]<-0
}
# The following outputs the values of the likelihood and theta and yields
# the value of theta where likelihood is a maximum
out = cbind(th,Lt)
Ltmax= max(Lt)
imax= which(Lt==max(Lt))
thmax= th[imax]
plot(th,Lt,lab=c(30,50,7),main="Box-Cox Transformations",
xlab=expression(theta),
ylab=expression(Lt(theta)))
cic = Ltmax-.5*qchisq(.95,1)
del= .01
iLtci = which(abs(Lt-cic)<=del)
iLtciL= min(iLtci)
iLtciU= max(iLtci)
thLci= th[iLtciL]
thUci= th[iLtciU]
abline(h=cic)
abline(v=thLci)
abline(v=thUci)
abline(v=thmax)

y2=x^.277
shapiro.test(y2)
qqnorm(y2,main="Normal Prob Plots of the CPUE data (catch)^.277", 
xlab="normal quantiles",ylab=expression(catch^.277,cex=.65))
qqline(y2)


#3
n = length(logx)
M= 10000
d= numeric(M)
for(i in 1:M)
{ d[i]= mean(sample(logx,replace=TRUE))
}
(bootmean= mean(d))
sd(d)/sqrt(M)

(sd(logx)/sqrt(50))

#4
(median(d))
(bootstd= sd(d))
(mad(d))

#5
.5/100*pnorm(3,3,1.5)

#VI
#2
M= 10000
d= numeric(M)
for(i in 1:M)
{ d[i]= mean(rexp(25, rate=.2))
}
y = sort(d)
n = length(d)
i = seq(1,n,1)
u = (i-.375)/(n+.25)
q = qnorm(u)
r = cor.test(q,y)

qqnorm(d,main="Normal Prob Plots of the 10000 Exponential sample means, lamda=.2)", 
xlab="normal quantiles",ylab="exponential sample means")
qqline(d)

#3
#a
pgamma(5.2,25,rate=5)-pgamma(4.8,25,rate=5)
#b
pnorm(.2)-pnorm(-.2)
#c
j=0
for(i in 1:M)
{ if (d[i] >= 4.8 & d[i] <= 5.2)
{ j=j+1}
}
(j/10000)
