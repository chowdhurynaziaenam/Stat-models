## Question_1
data<-read.table("tirereliability.txt", sep="\t", dec=".", header=TRUE)
attach(data)

library(survival)

## a

plot(wedge,survival)

survival
Surv(survival, complete)

model<-coxph(Surv(survival, complete)~wedge, data=data)
summary(model)
coef(model)

# -5.66 is the estimate of parameter

## b

newdata<-data.frame(wedge=0.6)
surf<-survfit(model, newdata=newdata, conf.type="plain")
summary(surf, times=1.00)

surv.curve<-summary(surf,times=seq(0,2,0.001))
plot(surf)
lines(surv.curve$time,surv.curve$surv, col="red", lwd=2)
lines(surv.curve$time,surv.curve$lower, col="blue", lwd=2)
lines(surv.curve$time,surv.curve$upper, col="blue", lwd=2)

# value of the survival function is 0.481

## c

hzratio<-exp(0.6*coef(model))/exp(1.6*coef(model))
hzratio

newdata<-data.frame(wedge=c(0.6,1.6))
exp.fit<-predict(model, newdata=newdata, type="risk")
exp.fit

exp(0.6*coef(model))/exp(mean(wedge)*coef(model))
exp(1.6*coef(model))/exp(mean(wedge)*coef(model))

hzratio<-exp.fit[1]/exp.fit[2]
hzratio

surf<-survfit(model, newdata=newdata)
plot(surf)
surv.curve<-summary(surf,times=seq(0,2,0.001))
lines(surv.curve$time,surv.curve$surv[,1], col="red", lwd=2)
lines(surv.curve$time,surv.curve$surv[,2], col="blue", lwd=2)

# hazard ratio is 288.5765

## d

model.3<-coxph(Surv(survival, complete)~wedge*peelForce+interBelt, data=data)
summary(model.3)
model.H0<-coxph(Surv(survival, complete)~peelForce+interBelt, data=data)
summary(model.H0)
anova(model.H0, model.3)

## e

newdata<-data.frame(wedge=0.6,peelForce=0.8,interBelt=0.7)
surf<-survfit(model.3, newdata=newdata, conf.type="plain")
summary(surf, times=1.00)


## Question_2

library(eha)

## a

model.weiph<-phreg(Surv(survival, complete)~wedge, data = data, dist="weibull")
summary(model.weiph)

beta<-coef(model.weiph)[1]
hzratio<-exp(beta*0.6)/exp(beta*1.6)
hzratio

x1<-0.6
x2<-1.6

hzratio<-exp(beta*x1)/exp(beta*x2)
hzratio

# hazard ratio is 528.249

p<-exp(coef(model.weiph)[3])
lambda<-exp(coef(model.weiph)[2])
beta<-coef(model.weiph)[1]

lambda.star1<-lambda/exp((x1*beta)/p)
lambda.star2<-lambda/exp((x2*beta)/p)

times<-seq(0,2,0.001)
survival.curve1<-1-pweibull(times,shape=p, scale=lambda.star1)
survival.curve2<-1-pweibull(times,shape=p, scale=lambda.star2)

plot(times, survival.curve1, type="n")
lines(times, survival.curve1, col="red", lwd=2)
lines(times, survival.curve2, col="blue", lwd=2)

## b 

p<-exp(coef(model.weiph)[3])
lambda<-exp(coef(model.weiph)[2])
beta<-coef(model.weiph)[1]
x<-1.6
lambda.star<-lambda/exp((x*beta)/p)
mu<-lambda.star*gamma(1+(1/p))
mu

## c 

p<-exp(coef(model.weiph)[3])
lambda<-exp(coef(model.weiph)[2])
beta<-coef(model.weiph)[1]
x<-1.6
lambda.star<-lambda/exp((x*beta)/p)

t.star<-rweibull(100000, shape=p, scale=lambda.star)
lowerbound<-quantile(t.star, c(0.1))
upperbound<-quantile(t.star, c(0.9))
lowerbound # 1.523983
upperbound # 2.048084

qweibull(0.1, shape=p, scale=lambda.star)
qweibull(0.9, shape=p, scale=lambda.star)

## d 

model.w<-phreg(Surv(survival, complete)~wedge*peelForce+interBelt, data = data, dist="weibull")
summary(model.w)

p<-exp(coef(model.w)[6])
lambda<-exp(coef(model.w)[5])
beta<-coef(model.w)[1:4]
x<-t(t(c(0.6,0.8,0.7,0.6*0.8)))
lambda.star<-lambda/exp((t(x)%*%beta)/p)
mu<-lambda.star*gamma(1+(1/p))
mu

survival.curve<-1-pweibull(1,shape=p, scale=lambda.star)
survival.curve

survival.curve<-1-pweibull(times,shape=p, scale=lambda.star) # 0.801
plot(times, survival.curve, type="n") 
lines(times, survival.curve, col="red",lwd=2) 



