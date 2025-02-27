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
