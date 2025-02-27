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
