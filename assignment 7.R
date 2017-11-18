#I.
mu=10700
sigma=1000
mu_0=10500
n=10
alpha = pnorm((mu-mu_0)*sqrt(n)/sigma)
(1-alpha)
(type_2=pnorm(qnorm(1-alpha)+sqrt(n)*(mu_0-mu)/sigma))

#II
(Y = mu_0+qnorm(1-.01)*sigma/sqrt(n))


mu=c(10600,10800,11000,11500)
(type_2=pnorm(qnorm(1-alpha)+sqrt(n)*(mu_0-mu)/sigma))

#III can assume normal data, sigma unknown, test if as H1 mu<m_0
mu_0=10
n=15
Y=8.7
S=2
alpha=.01
#1
(T=sqrt(n)*(Y-mu_0)/S)
(qt(alpha,n-1))
(pt(sqrt(n)*(Y-mu_0)/S,n-1))

#2 able to detect difference, type 2 error
(1-pt(-qt(1-alpha,n-1),n-1,sqrt(n)*(8.5-mu_0)/sigma))

#3 use chart
alpha=.05
beta=0.8
delta=1
sigma=2
27

#IV data is extremely normal
x=c(125,123,117,123,115,112,128,118,124,111,
116,109,125,120,113,123,112,118,121,118,122,
115,105,118,131)
shapiro.test(x)
y = sort(x)
n = length(x)
i = seq(1,n,1)
u = (i-.375)/(n+.25)
q = qnorm(u)
qqnorm(x,main="Normal Prob Plots of the blood sugar data", 
xlab="normal quantiles",ylab="Blood Sugar",cex=.65)
qqline(x)

#1 
alpha=.1
sigma_0=10
(TS=(n-1)*sd(x)^2/(sigma_0)^2)
(qchisq(alpha,n-1))

#2
sigs=c(5,6,7,8,9,10)
(type2=(1-pchisq(sigma_0^2/sigs^2*qchisq(alpha,n-1),n-1)))

#3 From chart of handout 11
(upper= sqrt(n-1)*sd(x)/sqrt(qchisq(.1/2,n-1)))

#V 
#1 use sign test
alpha=.05
mu_tilda=120
Xi=x-mu_tilda
S_plus=0
for (i in 1:n){
	if (Xi[i] >0) S_plus= S_plus+1 }
qbinom(alpha,n,.5)
pbinom(S_plus,n,.5)

#2 use Wilcox signed rank test
y=rep(median(x),length(x))
wilcox.test(x,y,alternative="l",paired=T,conf.level=.95)
qsignrank(.05,length(x),TRUE)
psignrank(qsignrank(.05,length(x),TRUE),length(x),TRUE)

#3 Use table, with n=25, alpha=.05, k=8
k=8
y=sort(x)
y[length(y)-8]

#VI
n=50
y=46
p_hat=46/50
p_0=.8
alpha=.05
#1 Use A-C test since n>40
y_tilda=y+(qnorm(alpha/2)^2)/2
n_tilda=n+qnorm(alpha/2)^2
p_tilda=y_tilda/n_tilda
(upper=p_tilda+qnorm(1-alpha/2)*sqrt(p_tilda*(1-p_tilda)/n_tilda))
(lower=p_tilda-qnorm(1-alpha/2)*sqrt(p_tilda*(1-p_tilda)/n_tilda))

#2 use asymptotic test
(TS=(p_hat-p_0)/sqrt(p_0*(1-p_0)/n))
qnorm(1-alpha)
(1-pnorm(TS))

#3 power function
p_1=c(.75,.8,.85,.9,.95)
(power=1-pnorm(qnorm(1-alpha)*sqrt(p_0*(1-p_0)/p_1/(1-p_1))+sqrt(n)*
(p_0-p_1)/sqrt(p_1*(1-p_1))))

#4
p_1=.9
beta=.8
(num=((qnorm(1-alpha)*sqrt(p_0*(1-p_0))+qnorm(1-beta)*p_1*(1-p_1))/(p_1-p_0))^2)