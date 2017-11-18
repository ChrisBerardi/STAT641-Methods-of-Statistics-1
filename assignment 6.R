#Problem I:

add = c(151, 450, 124, 235, 357, 110, 302, 671, 118, 115,
275, 275, 2550, 243, 201, 199, 130, 119, 92, 91, 92, 98, 1650, 
1200, 1180, 900, 700, 460, 340, 330)

noadd = c(246, 268, 275, 348, 305, 311, 206, 279, 426,
269, 257, 299, 337, 329, 319, 312, 327, 342, 351, 205, 151,
426, 154, 353, 396, 441, 254, 263, 278, 281)

#With additive
#1
y=add
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
shapiro.test(add^-.36)

#flip since g is a decreasing function
(mean(add^-.36)+2.22*sd(add^-.36))^(1/-.36)

y = sort(add^-.36)
n = length(add^-.36)
i = seq(1,n,1)
u = (i-.375)/(n+.25)
q = qnorm(u)
qqnorm(add^-.36,main="Normal Prob Plots of the Data with additive", 
xlab="normal quantiles",ylab="Additive",cex=.65)
qqline(add^-.36)


#2
n= length(add)
thest = mean(add)
B = 9999
thestS = numeric(B)
thestS = rep(0,times =B)
for (i in 1:B)
thestS[i] = mean(sample(add,replace=T))
RS= sort(thestS-thest)
LRS = RS[250]
URS = RS[9750]
thL = thest-URS
thU = thest-LRS

#3
#using chart, r=10
y=sort(add)
(lower = y[10])
(upper= y[30-10+1])

#No additive
#normal looks good
shapiro.test(noadd)
y = sort(noadd)
n = length(noadd)
i = seq(1,n,1)
u = (i-.375)/(n+.25)
q = qnorm(u)
qqnorm(noadd,main="Normal Prob Plots of the Data without additive", 
xlab="normal quantiles",ylab="No Additive",cex=.65)
qqline(noadd)

#1
y=sort(noadd)
me=mean(noadd)
st=sd(noadd)
mean(noadd)-2.22*sd(noadd)

#2
#using noadd is normal and the standrd pivot for mean
(lower=me-dt(.025,29)*st/sqrt(n))
(upper=me+dt(.025,29)*st/sqrt(n))
#3
#using chart, r=10
(lower = y[10])
(upper= y[30-10+1])

#Problem II:
#1 Use sensoring stuff
x = c(2.526,  2.546,  2.628,  2.669,  2.869,  2.710,  2.731,  2.751,  2.771,  
     2.772,  2.782,  2.789,  2.793, 2.834,  2.844, 2.854,  2.875,  
     2.876,  2.895,  2.916,  2.919,  2.957,  2.977,  2.988,  
     3,  3,  3,  3)
 library(MASS)
 library(survival)
stcens = c(rep(1,24),rep(0,4))
Surv(x, stcens)
cords.surv <- survfit(Surv(x, stcens) ~ 0,conf.type="log-log")
summary(cords.surv)
print(cords.surv,print.rmean=TRUE)
plot(cords.surv,conf.int=F,log=FALSE,ylim=c(0,1),
main="Kaplan-Meier Estimator of Survival Function",xlab="Strength of Cord", 
ylab="Survival Function")

#2 Assuming weibull

stcens = c(rep(1,24),rep(0,4))
Surv(x, stcens)
cords.surv <- survreg(Surv(x, stcens) ~ 1,dist='weibull')
summary(cords.surv)

gamma= 1/.0439
alpha= exp(1.07)
cvD= alpha*gamma(1+1/gamma)
V=
n = length(x)
B = 9999
W = matrix(0,B,n)
cv = numeric(B)
cv = rep(0,B)
a = numeric(B)
a = rep(0,B)
b = numeric(B)
b = rep(0,B)
mleest = matrix(0,B,2)
{
for (i in 1:B)
W[i,] = rweibull(n,alpha,gamma)
}
{
for (i in 1:B)
mleest[i,] = coef(fitdistr(W[i,],"weibull"))
}
a = mleest[,1]
b = mleest[,2]
cv = a*gamma(1+1/b)
R = cv-cvD
R = sort(R)
L = R[250]
U = R[9750]
ci = c(cvD-U, cvD-L)
ci



#Problem IV:

w = c(.18, 3.1,   4.2,    6.0,    7.5,    8.2,    8.5,   10.3,   10.6,   24.2,

 29.6,  31.7,  41.9,   44.1,   49.5,   50.1,   59.7,   61.7,   64.4,   69.7
,
 70.0,  77.8,  80.5,   82.3,   83.5,   84.2,   87.1,   87.3,   93.2,  103.4,
104.6, 105.5, 108.8,  112.6,  116.8,  118.0,  122.3,  123.5,  124.4,  125.4, 

129.5, 130.4, 131.6,  132.8,  133.8,  137.0,  140.2,  140.9,  148.5,  149.2
, 
152.2, 152.9, 157.7,  160.0,  163.6,  166.9,  170.5,  174.9,  177.7,  179.2
, 
183.6, 183.8, 194.3,  195.1,  195.3,  202.6,  220.0,  221.3,  227.2,  251.0,
266.5, 267.9, 269.2,  270.4,  272.5,  285.9,  292.6,  295.1,  301.1,  304.3,

316.8, 329.8, 334.1,  346.2,  351.2,  353.3,  369.3,  372.3,  381.3,  393.5,

451.3, 461.5, 574.2,  656.3,  663.0,  669.8,  739.7,  759.6,  894.7,  974.9)

#A
x=w
n = length(x)
mn= mean(x)
thhat = exp(-300/mn)
S2 = var(x)
Vest= ((300/(mn**2))**2)*(exp(-600/mn))*(S2/n)
R = 9999
z = numeric(R)
z = rep(0,times =R)
for (i in 1:R)
{t= sample(x,replace=T)
V= ((300/(mean(t)**2))**2)*(exp(-600/mean(t)))*(var(t)/n)
z[i] = (exp(-300/mean(t))-thhat)/sqrt(V)}
z = sort(z)
L = z[25]
U = z[9975]
(thL = thhat-U*sqrt(Vest))
(thU = thhat-L*sqrt(Vest))

l = .22-qnorm(.995)*sqrt(.22*(1-.22))/sqrt(n)
u= .22+qnorm(.995)*sqrt(.22*(1-.22))/sqrt(n)

#b
#not normal
n  =  length(w)
lam  =  mean(w)
w  =  sort(w)
z  =  1-exp(-w/lam)   #computes F0(X(i))
i  =  seq(1,n,1)
a1i  =  (2*i-1)*log(z)
a2i  =  (2*n+1-2*i)*log(1-z)
s1  =  sum(a1i)
s2  =  sum(a2i)
AD  =  -n-(1/n)*(s1+s2)
MAD  =  AD^2*(1+.6/n)
MAD
j  =  (i-.5)/n
quexp  =  -log(1-j)
plot(quexp,w,xlab="Exponential Quantiles",ylab="Sample Quantiles")
title("Exponential Probability Plot")
abline(lm(w~quexp))

b_hat=mean(w)
(-b_hat*((2*n)/qchisq(.99,2*n))*log(.95))

#C
(l_bound=2*n*mean(w)/qchisq(.975,2*n))
(u_bound=2*n*mean(w)/qchisq(.025,2*n))

Problem V:

y = c(19.7, 21.6, 21.9, 23.5, 24.2, 24.4, 24.9, 25.1,
 26.4, 26.9, 27.6, 27.7, 
27.9, 28.4, 29.8, 30.7,
 31.1, 31.1, 31.7, 31.8, 32.6, 34.0, 34.8, 34.9,
 35.1, 
36.6, 37.0, 37.7, 38.7, 38.7, 39.0, 39.6,
 40.0, 41.4, 41.4, 41.8, 42.2, 43.5, 
44.5, 45.0,
 45.5, 45.9, 46.3, 46.7, 46.7, 47.0, 47.0, 47.4,
 47.6, 48.6, 48.8, 
57.9, 58.3, 67.9, 84.2, 97.3)

y = sort(y)
n = length(y)
i = 1:n
ui = (i-.5)/n
fitdistr(y,'gamma')
QG = qgamma(ui,9.34843040,.23817682 )
plot(QG,y,abline(c(0,1)))
#text(7,25.0,expression(hat(Q)(u) == .14 + .9445*Q(u[i])))
abline(h=0)


#A
n= length(y)
m= sum(y<50)
thest = m/n
B = 9999
thestS = numeric(B)
thestS = rep(0,times =B)
for (i in 1:B)
thestS[i] = sum(sample(y,replace=T)<50)/56
RS = sort(thestS-thest)
LRS = RS[250]
URS = RS[9750]
thL = thest-URS
thU = thest-LRS

l = .95-qnorm(.975)*sqrt(.95*(1-.95))/sqrt(n)
u= .95+qnorm(.975)*sqrt(.95*(1-.95))/sqrt(n)

#B
n= length(y)
thest = median(y)
B = 9999
thestS = numeric(B)
thestS = rep(0,times =B)
for (i in 1:B)
thestS[i] = median(sample(y,replace=T))
RS= sort(thestS-thest)
LRS = RS[250]
URS = RS[9750]
thL = thest-URS
thU = thest-LRS
