##Question_1
data<-read.table("chromoabnormal.txt", sep="\t", dec=".", header=TRUE)
attach(data)

##a

y<-ca/cells

plot(ca/cells~doseamt)
plot(ca/cells~doserate)
coplot(ca/cells~doseamt| doserate)

model.main<-glm(ca~offset(log(cells))+doseamt+doserate, family=poisson(link="log"), data=data)
summary(model.main)
model<-glm(ca~offset(log(cells))+doseamt+doserate+doseamt:doserate, family=poisson(link="log"), data=data)
summary(model)

AIC(model.main)
AIC(model)

model.matrix(model)

newdata<-data.frame(cells=64070,doseamt=4, doserate=0.75)
mu.hat<-predict(model, newdata=newdata, type="response")
mu.hat # 311.3886

ratio.estimate<-mu.hat/newdata$cells
ratio.estimate

xf<-t(cbind(1,4,0.75,4*0.75))
exp(log(64070)+t(xf)%*%coef(model))


##b

newdata<-data.frame(cells=median(data$cells),doseamt=4, doserate=0.75)
pred<-predict(model, newdata=newdata, type="response")
pred

xf<-t(cbind(1,4,0.75,4*0.75))

Var.eYf<-pred*(1+pred*t(xf)%*%vcov(model)%*%xf)

lower.Yf<-pred-qnorm(0.9)*sqrt(Var.eYf)
upper.Yf<-pred+qnorm(0.9)*sqrt(Var.eYf)
lower.Yf
upper.Yf

ratio.prediction<-pred/newdata$cells
ratio.prediction # 0.004860

Var.eZf<-((1/newdata$cells)^2)*Var.eYf

lower.Zf<-ratio.prediction-qnorm(0.9)*sqrt(Var.eZf)
upper.Zf<-ratio.prediction+qnorm(0.9)*sqrt(Var.eZf)
lower.Zf
upper.Zf


model.anova<-glm(ca~offset(log(cells))+factor(doseamt)*factor(doserate), family=poisson(link="log"), data=data)
summary(model.anova)

model.matrix(model.anova)


##c

model.H0<-glm(ca~offset(log(cells))+doseamt, family=quasipoisson(link="log"), data=data)
model.H1<-glm(ca~offset(log(cells))+doseamt*doserate, family=quasipoisson(link="log"), data=data)
summary(model.H1)
anova(model.H0, model.H1, test="F")
anova(model.H0, model.H1, test="F")$F[2] # 5.363

model.H0<-glm(ca~offset(log(cells))+doseamt, family=poisson(link="log"), data=data)
model.H1<-glm(ca~offset(log(cells))+doseamt*doserate, family=poisson(link="log"), data=data)
anova(model.H0, model.H1, test="Chi")

##d

model.H1<-glm(ca~offset(log(cells))+doseamt*doserate, family=quasipoisson(link="log"), data=data)

MSE.i<-mean((ca-predict(model, newdata=data, type="response"))^2)
MSE.ii<-mean((ca-predict(model.H1, newdata=data, type="response"))^2)
library(MASS)
model.NB<-glm.nb(ca~offset(log(cells))+doseamt*doserate, data=data)
summary(model.NB)
MSE.iii<-mean((ca-predict(model.NB, newdata=data, type="response"))^2)

MSE.i
MSE.ii
MSE.iii

AIC(model)  
AIC(model.H1) 
AIC(model.NB) 

plot(fitted(model.H1, type="response"), residuals(model.H1, type="response"))
plot(fitted(model.NB, type="response"), residuals(model.NB, type="response"))

plot(fitted(model, type="response"), residuals(model, type="pearson")^2)
plot(fitted(model.H1, type="response"), residuals(model.H1, type="pearson")^2)
plot(fitted(model.NB, type="response"), residuals(model.NB, type="pearson")^2)

# poisson distribution



##Question_2
data<-read.table("applejuiceCRA7152.txt", sep="\t", dec=".", header=TRUE)
attach(data)

##a,b

model<-glm(Growth~pH+Nisin+Temperature+Brix, data=data, family=binomial(link="logit"))
summary(model)


model.full<-glm(Growth~pH*Nisin*Temperature*Brix, data=data, family=binomial(link="logit"))
summary(model.full)


training<-data[-c(7,34,56,78),]
new<-data[c(7,34,56,78),]
plot(jitter(fitted(model.full, type="response")),jitter(Growth))
model.full<-glm(Growth~pH*Nisin*Temperature*Brix, data=training, family=binomial(link="logit"))
summary(model.full)
predict(model.full, newdata=new, type="response")


step(model)
AIC(model)

model.12<-glm(Growth~pH*Nisin+Temperature+Brix, data=data, family=binomial(link="logit"))
summary(model.12)

model.P<-glm(Growth~pH+Nisin+Temperature+Brix, data=data, family=binomial(link="probit"))
summary(model.P)
model.C<-glm(Growth~pH+Nisin+Temperature+Brix, data=data, family=binomial(link="cauchit"))
summary(model.C)
model.cll<-glm(Growth~pH+Nisin+Temperature+Brix, data=data, family=binomial(link="cloglog"))
summary(model.cll)


AIC(model.12)
AIC(model.P)
AIC(model.C)
AIC(model.cll)

#Probit link 

plot(Nisin, fitted(model, type="response"))
plot(Nisin, fitted(model.C, type="response"))

plot(fitted(model, type="response"), fitted(model.C, type="response"))

mean(residuals(model,type="response")^2)
mean(residuals(model.P,type="response")^2)
mean(residuals(model.C,type="response")^2)
mean(residuals(model.cll,type="response")^2)


newdata<-data.frame(pH=4.5,Nisin=20,Temperature=30,Brix=17)
predict(model, newdata=newdata, type="response")
predict(model.P, newdata=newdata, type="response") # 0.01011

newdata<-data.frame(pH=4.5,Nisin=0:70,Temperature=30,Brix=17)
pred.L<-predict(model, newdata=newdata, type="response")
pred.P<-predict(model.P, newdata=newdata, type="response")

plot(jitter(Nisin, 0.4), jitter(Growth,0.4))
lines(newdata$Nisin, pred.L, lwd=3)
lines(newdata$Nisin, pred.P, lwd=3, col="red")


##c

newdata<-data.frame(pH=4.5,Nisin=20,Temperature=30,Brix=17)
eta<-predict(model, newdata=newdata, type="link", se.fit = TRUE) 
link.lowerbound<-eta$fit-qnorm(0.975)*eta$se.fit
link.upperbound<-eta$fit+qnorm(0.975)*eta$se.fit

newdata<-data.frame(pH=4.5,Nisin=20,Temperature=30,Brix=17)
eta<-predict(model.P, newdata=newdata, type="link", se.fit = TRUE) 
link.lowerbound<-eta$fit-qnorm(0.975)*eta$se.fit
link.upperbound<-eta$fit+qnorm(0.975)*eta$se.fit

pnorm(link.lowerbound) # 0.01412
pnorm(link.upperbound) # 0.36093


## d)

model<-glm(Growth~pH+Nisin+Temperature+Brix, data=data, family=binomial(link="logit"))
summary(model)
newdata<-data.frame(pH=4.5,Nisin=20,Temperature=30,Brix=17)

mu.f<-predict(model, newdata=newdata, type="response")
YS.pred<-100*mu.f

mu.hat<-predict(model, newdata=data, type="response")
n<-dim(data)[1]

e.b<-numeric()

for(b in 1:1000){
  
  yb<-numeric()
  for(i in 1:n){
    
    yb[i]<-sample(0:1,1,prob=c(1-mu.hat[i],mu.hat[i]))
    
  }
  
  model.b<-glm(yb~pH+Nisin+Temperature+Brix, family=binomial(link="logit"), data=data)
  newdata<-data.frame(pH=4.5,Nisin=20,Temperature=30,Brix=17)
  mu.fb<-predict(model.b, newdata=newdata, type="response")
  YS.predB<-100*mu.fb
  
  yf.b<-sample(0:1,100,prob=c(1-mu.f,mu.f), replace=TRUE)
  
  e.b[b]<-sum(yf.b)-YS.predB
  
}

var.error<-var(e.b)
var.error

z<-qnorm(c(0.1), lower.tail=FALSE)
lower.bound<-YS.pred-z*sqrt(var.error)
upper.bound<-YS.pred+z*sqrt(var.error)
lower.bound # -0.241699 
upper.bound # 22.61417 
