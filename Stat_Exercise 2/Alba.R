data<-read.table("Alba.txt", sep="\t", dec=".", header=TRUE)
attach(data)

plot(Dose, DryMatter)
coplot(DryMatter~Dose| Herbicide)

##a)

model<-glm(DryMatter~Dose*Herbicide, family=Gamma(link="inverse"), data=data)
summary(model)

newdata<-data.frame(Dose=50, Herbicide="Glyphosate")
predict(model, newdata=newdata, type="response")

newdata<-expand.grid(Dose=0:700, Herbicide=levels(Herbicide))
pred<-predict(model, newdata=newdata, type="response")

plot(Dose, DryMatter)
lines(newdata$Dose[1:701],pred[1:701], lwd=2)
lines(newdata$Dose[-(1:701)],pred[-(1:701)], lwd=2, col="red")

##b)

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


##c)

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
