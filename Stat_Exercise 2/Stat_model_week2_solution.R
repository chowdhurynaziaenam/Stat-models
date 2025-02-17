### Question 1

rats<-read.table("ratstime.txt", sep="\t", dec=".", header=TRUE)
attach(rats)

boxplot(time~poison+treat)
interaction.plot(poison,treat, time)


### a)

model<-lm(time~poison+treat, data=rats)
summary(model)

interaction.plot(poison,treat,fitted(model))

newdata<-expand.grid(poison=levels(poison),treat=levels(treat))
Xc<-model.matrix(~poison+treat, data=newdata)
betahat<-coef(model)

xAmean<-(Xc[1,]+Xc[2,]+Xc[3,])/3
xBmean<-(Xc[4,]+Xc[5,]+Xc[6,])/3
xCmean<-(Xc[7,]+Xc[8,]+Xc[9,])/3
xDmean<-(Xc[10,]+Xc[11,]+Xc[12,])/3

xAmean%*%betahat
xBmean%*%betahat
xCmean%*%betahat
xDmean%*%betahat


kAB<-xAmean-xBmean
K<-cbind(kAB)
q<-1
Wald.AB<-(t(t(K)%*%betahat)%*%solve(t(K)%*%vcov(model)%*%K)%*%t(K)%*%betahat)/q
Wald.AB
p.AB<-pf(Wald.AB, 1, summary(model)$df[2], lower.tail = FALSE)
p.AB


kAC<-xAmean-xCmean
K<-cbind(kAC)
q<-1
Wald.AC<-(t(t(K)%*%betahat)%*%solve(t(K)%*%vcov(model)%*%K)%*%t(K)%*%betahat)/q
Wald.AC
p.AC<-pf(Wald.AC, 1, summary(model)$df[2], lower.tail = FALSE)
p.AC

kAD<-xAmean-xDmean
K<-cbind(kAD)
q<-1
Wald.AD<-(t(t(K)%*%betahat)%*%solve(t(K)%*%vcov(model)%*%K)%*%t(K)%*%betahat)/q
Wald.AD
p.AD<-pf(Wald.AD, 1, summary(model)$df[2], lower.tail = FALSE)
p.AD

kBC<-xBmean-xCmean
K<-cbind(kBC)
q<-1
Wald.BC<-(t(t(K)%*%betahat)%*%solve(t(K)%*%vcov(model)%*%K)%*%t(K)%*%betahat)/q
Wald.BC
p.BC<-pf(Wald.BC, 1, summary(model)$df[2], lower.tail = FALSE)
p.BC

kBD<-xBmean-xDmean
K<-cbind(kBD)
q<-1
Wald.BD<-(t(t(K)%*%betahat)%*%solve(t(K)%*%vcov(model)%*%K)%*%t(K)%*%betahat)/q
Wald.BD
p.BD<-pf(Wald.BD, 1, summary(model)$df[2], lower.tail = FALSE)
p.BD

kCD<-xCmean-xDmean
K<-cbind(kCD)
q<-1
Wald.CD<-(t(t(K)%*%betahat)%*%solve(t(K)%*%vcov(model)%*%K)%*%t(K)%*%betahat)/q
Wald.CD
p.CD<-pf(Wald.CD, 1, summary(model)$df[2], lower.tail = FALSE)
p.CD

Wald<-c(Wald.AB,Wald.AC,Wald.AD,Wald.BC,Wald.BD,Wald.CD)
p.value<-c(p.AB,p.AC,p.AD,p.BC,p.BD,p.CD)

Wald
p.value

#

library(multcomp)
K<-cbind(kAB,kAC,kAD,kBC,kBD,kCD)
meanpairwise<-glht(model, linfct=t(K))
summary(meanpairwise)
summary(meanpairwise)$test$tstat^2
Wald


### b)

# predictive hypothesis testing

X<-model.matrix(model)
sigma2<-sigma(model)^2
kt<-t(kAB)
pred<-kt%*%betahat
T.AB<-pred/sqrt(sigma2*(2+(kt)%*%solve(t(X)%*%X)%*%t(kt)))   ### mean(Y)~N(mean(mu), sigma^2)
#k<-3
#T.AB<-pred/sqrt(sigma2/k+sigma2/k+sigma2*(kt)%*%solve(t(X)%*%X)%*%t(kt))   ### mean(Y)~N(mean(mu), sigma^2/k)
T.AB
d.AB<-2*pt(abs(T.AB),df=summary(model)$df[2], lower.tail = FALSE)
d.AB

kt<-t(kAC)
pred<-kt%*%betahat
T.AC<-pred/sqrt(sigma2*(2+(kt)%*%solve(t(X)%*%X)%*%t(kt)))
T.AC
d.AC<-2*pt(abs(T.AC),df=summary(model)$df[2], lower.tail = FALSE)
d.AC

kt<-t(kAD)
pred<-kt%*%betahat
T.AD<-pred/sqrt(sigma2*(2+(kt)%*%solve(t(X)%*%X)%*%t(kt)))
T.AD
d.AD<-2*pt(abs(T.AD),df=summary(model)$df[2], lower.tail = FALSE)
d.AD

kt<-t(kBC)
pred<-kt%*%betahat
T.BC<-pred/sqrt(sigma2*(2+(kt)%*%solve(t(X)%*%X)%*%t(kt)))
T.BC
d.BC<-2*pt(abs(T.BC),df=summary(model)$df[2], lower.tail = FALSE)
d.BC

kt<-t(kBD)
pred<-kt%*%betahat
T.BD<-pred/sqrt(sigma2*(2+(kt)%*%solve(t(X)%*%X)%*%t(kt)))
T.BD
d.BD<-2*pt(abs(T.BD),df=summary(model)$df[2], lower.tail = FALSE)
d.BD

kt<-t(kCD)
pred<-kt%*%betahat
T.CD<-pred/sqrt(sigma2*(2+(kt)%*%solve(t(X)%*%X)%*%t(kt)))
T.CD
d.CD<-2*pt(abs(T.CD),df=summary(model)$df[2], lower.tail = FALSE)
d.CD

T<-c(T.AB,T.AC,T.AD,T.BC,T.BD,T.CD)
d.value<-c(d.AB,d.AC,d.AD,d.BC,d.BD,d.CD)

T
d.value

y<-seq(0,1,0.01)
plot(y,dnorm(y, mean=as.numeric(xAmean%*%betahat), sd=sqrt(sigma2)), type="l", lwd=3)
lines(y,dnorm(y, mean=as.numeric(xBmean%*%betahat), sd=sqrt(sigma2)), col="red", lwd=3)

ya<-rnorm(100000, mean=as.numeric(xAmean%*%betahat), sd=sqrt(sigma2))
yb<-rnorm(100000, mean=as.numeric(xBmean%*%betahat), sd=sqrt(sigma2))

2*(1-sum(ya<yb)/100000)

### Simulation d-value

mu<-2.8
ya<-rnorm(100000, mean=0, sd=1)
yb<-rnorm(100000, mean=mu, sd=1)
y<-seq(-3,5,0.01)
plot(y,dnorm(y, mean=0, sd=1), type="l", lwd=3)
lines(y,dnorm(y, mean=mu, sd=1), col="red", lwd=3)

y<-c(ya,yb)
#x<-rep(0:1, each=10)
x<-rep(0:1, each=100000)
model<-lm(y~x)
anova(model, test="F")

kt<-t(c(0,1))
betahat<-coef(model)
X<-model.matrix(model)
pred<-kt%*%betahat
sigma2<-sigma(model)^2
T<-pred/sqrt(sigma2*(2+(kt)%*%solve(t(X)%*%X)%*%t(kt)))
T
d<-2*pt(abs(T),df=summary(model)$df[2], lower.tail = FALSE)
d

2*(1-sum(ya<yb)/100000)

### c)

main.normal<-glm(time~poison+treat, family=gaussian(link="identity"), data=rats)
main.normalL<-glm(time~poison+treat, family=gaussian(link="log"), data=rats)
main.normalI<-glm(time~poison+treat, family=gaussian(link="inverse"), data=rats)

par(mfrow=c(2,2))
interaction.plot(poison, treat, time)
interaction.plot(poison, treat, fitted(main.normal, type="response"))
interaction.plot(poison, treat, fitted(main.normalL, type="response"))
interaction.plot(poison, treat, fitted(main.normalI, type="response"))

# MSE values

mean(residuals(main.normal, type="response")^2)
mean(residuals(main.normalL, type="response")^2)
mean(residuals(main.normalI, type="response")^2)

AIC(main.normal)
AIC(main.normalL)
AIC(main.normalI)


### Question 2

data<-read.table("Alba.txt", sep="\t", dec=".", header=TRUE)
attach(data)

plot(Dose, DryMatter)
coplot(DryMatter~Dose| Herbicide)


### a)

model<-glm(DryMatter~Dose*Herbicide, family=Gamma(link="inverse"), data=data)
summary(model)

newdata<-data.frame(Dose=50, Herbicide="Glyphosate")
predict(model, newdata=newdata, type="response")

newdata<-expand.grid(Dose=0:700, Herbicide=levels(Herbicide))
pred<-predict(model, newdata=newdata, type="response")

plot(Dose, DryMatter)
lines(newdata$Dose[1:701],pred[1:701], lwd=2)
lines(newdata$Dose[-(1:701)],pred[-(1:701)], lwd=2, col="red")

### b)

newdata<-data.frame(Dose=50, Herbicide="Glyphosate")

pred<-predict(model, newdata=newdata, type="response")
pred

model.matrix(model)
xf<-t(t(c(1,50,1,50)))

phi<-summary(model)$dispersion
Var.Yf<-phi*(pred^2)
D.f<--(pred^2)
Var.ef<-Var.Yf+(D.f^2)*t(xf)%*%vcov(model)%*%xf

lower.yf<-pred-qnorm(0.9)*sqrt(Var.ef)
upper.yf<-pred+qnorm(0.9)*sqrt(Var.ef)

lower.yf
upper.yf

###

newdata<-data.frame(Dose=0:700, Herbicide="Glyphosate")
pred<-predict(model, newdata=newdata, type="response")

Xf<-cbind(1,newdata$Dose,1,newdata$Dose)
Xf

phi<-summary(model)$dispersion
Var.Yf<-phi*(pred^2)
D.f<--(pred^2)
Var.ef<-Var.Yf+(D.f^2)*diag(Xf%*%vcov(model)%*%t(Xf))

lower.yf<-pred-qnorm(0.9)*sqrt(Var.ef)
upper.yf<-pred+qnorm(0.9)*sqrt(Var.ef)


plot(c(0,700),c(0,6), type="n")
points(Dose, DryMatter)
lines(newdata$Dose,pred, lwd=2)
lines(newdata$Dose,lower.yf, lwd=2, col="red")
lines(newdata$Dose,upper.yf, lwd=2, col="red")


### c)

modelIG<-glm(DryMatter~Dose*Herbicide, family=inverse.gaussian(link="1/mu^2"), data=data)
summary(modelIG)

newdata<-data.frame(Dose=50, Herbicide="Glyphosate")

eta<-predict(modelIG, newdata=newdata, type="link", se.fit=TRUE)
link.lowerbound<-eta$fit-qnorm(0.975)*eta$se.fit
link.upperbound<-eta$fit+qnorm(0.975)*eta$se.fit
link.lowerbound
link.upperbound

lower.mu<-sqrt(1/link.upperbound)
upper.mu<-sqrt(1/link.lowerbound)

lower.mu
upper.mu
