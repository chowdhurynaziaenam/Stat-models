data<-read.table("NitrogenYield.txt", sep="\t", dec=".", header=TRUE)
attach(data)
plot(Nitrogen, Yield)
library(nlme)

## a 
model.poly<-lm(Yield~Nitrogen+I(Nitrogen^2), data=data)
summary(model.poly)
coef(model.poly)[3]

plot(Nitrogen, Yield)
lines(Nitrogen, fitted(model.poly), col="red")

## maximum likelihood estimate is - 0.0022237

## b 
model.expon<-glm(Yield~log(Nitrogen),family=gaussian(link="log"), data=data)
summary(model.expon)
predict(model.expon, newdata=data.frame(Nitrogen=150), type="response")

plot(Nitrogen, Yield)
lines(Nitrogen, fitted(model.poly), col="red")
lines(Nitrogen, fitted(model.expon, type="response"), col="blue")

## maximum likelihood estimate is 90.66521.

## c
model.asymp<-nls(Yield~SSasymp(Nitrogen, Asym,R0,lrc), data=data)
summary(model.asymp)

plot(Nitrogen, Yield)
lines(Nitrogen, fitted(model.poly), col="red", lwd=3)
lines(Nitrogen, fitted(model.expon, type="response"), col="blue")
lines(Nitrogen, fitted(model.asymp, type="response"), col="black", lwd=3)

## maximum likelihood estimate is 102.0634. 

## d 
model.mm<-nls(Yield~SSmicmen(Nitrogen, Vm, K), data=data)
summary(model.mm)
newdata<-data.frame(Nitrogen=150)
predict(model.mm,newdata=newdata)
## Maximum likelihood estimate is 90.57802

plot(Nitrogen, Yield)
lines(Nitrogen, fitted(model.poly), col="red")
lines(Nitrogen, fitted(model.expon, type="response"), col="blue")
lines(Nitrogen, fitted(model.asymp, type="response"), col="black")
lines(Nitrogen, fitted(model.mm, type="response"), col="green", lwd=3)

gmodel<-glm(Yield~I(1/Nitrogen), family=gaussian(link="inverse"), data=data)
lines(Nitrogen, fitted(gmodel, type="response"), col="brown", lwd=3)

AIC(model.poly)
AIC(model.expon)
AIC(model.asymp)
AIC(model.mm)

## e
model.asymp<-nls(Yield~SSasymp(Nitrogen,Asym,R0,lrc), data=data)
summary(model.asymp)

beta<-coef(model.asymp)
cov.beta<-vcov(model.asymp)
predict(model.asymp, newdata=data.frame(Nitrogen=150), type="response")
library(mvtnorm)
beta.star<-rmvnorm(1000, mean = beta, sigma = cov.beta)
Asym<-beta.star[,1]
R0<-beta.star[,2]
lrc<-beta.star[,3]

mu.star<-Asym+(R0-Asym)*exp(-exp(lrc)*newdata$Nitrogen)
conf.lowerbound<-quantile(mu.star, c(0.025))
conf.upperbound<-quantile(mu.star, c(0.975))
conf.lowerbound
conf.upperbound

sigma2<-sigma(model.asymp)[1]^2
yf.star<-rnorm(10000, mean=mu.star, sd=sqrt(sigma2))
pred.lowerbound<-quantile(yf.star, c(0.1))
pred.upperbound<-quantile(yf.star, c(0.9))

pred.lowerbound # 83.16026 
pred.upperbound # 99.13124 

newdata<-data.frame(Nitrogen=0:200)
pred.lowerbound<-numeric()
pred.upperbound<-numeric()

for(i in 1:dim(newdata)[1]){
  
  mu.star<-Asym+(R0-Asym)*exp(-exp(lrc)*newdata$Nitrogen[i])

  sigma2<-sigma(model.asymp)[1]^2
  yf.star<-rnorm(10000, mean=mu.star, sd=sqrt(sigma2))
  pred.lowerbound[i]<-quantile(yf.star, c(0.1))
  pred.upperbound[i]<-quantile(yf.star, c(0.9))
  
}

plot(Nitrogen, Yield)
lines(Nitrogen, fitted(model.asymp, type="response"), col="black")
lines(newdata$Nitrogen, pred.lowerbound, col="red")
lines(newdata$Nitrogen, pred.upperbound, col="red")
