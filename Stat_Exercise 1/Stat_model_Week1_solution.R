##Question_1

library(tidyverse)
library(plotly)

data<-read.table("paper.txt", sep="\t", dec=".", header=TRUE)
attach(data)

coplot(strength~hardwood |pressure)
plot(hardwood,strength)
plot(pressure, strength)

plot_ly()%>%
  add_markers(x=~hardwood,y=~pressure, z=~strength) 

model.main<-lm(strength~hardwood+pressure)
summary(model.main)

#a
coef(model.main)
coef(model.main)[2]

#b
sigma(model.main)^2

#c
fitted(model.main)
fitted(model.main)[36]
tail(data)
predict(model.main, newdata=data.frame(hardwood=8, pressure=650))

#d
newdata<-data.frame(hardwood=8,pressure=550)
mu.hat<-predict(model.main, newdata=newdata, interval="confidence")

mu.hat[1] # maximum likelihood point estimate

#e
newdata<-data.frame(hardwood=8,pressure=550)
y.hat<-predict(model.main, newdata=newdata, interval="prediction", level=0.8)

y.hat[2] # lower bound of the 80% prediction interval

#f
model.1<-lm(strength~hardwood)
anova(model.1, model.main,test="F")
anova(model.1, model.main,test="F")$F[2] ## ans

summary(model.main)
model.2<-lm(strength~pressure)
anova(model.2, model.main,test="F")
anova(model.1,model.2, test="F")$F[2]


##Question_2

data<-read.table("makiwaraboard.txt", sep="\t", dec=".", header=TRUE)
attach(data)

plot_ly()%>%
  add_markers(x=~WoodType,y=~BoardType, z=~Deflection)

model.main<-lm(Deflection~factor(WoodType)+factor(BoardType))
summary(model.main)
interaction.plot(WoodType, BoardType, fitted(model.main))

model.12<-lm(Deflection~factor(WoodType)+factor(BoardType)+factor(WoodType):factor(BoardType))
summary(model.12)
interaction.plot(WoodType, BoardType, fitted(model.12))

##a
newdata<-data.frame(WoodType="Oak", BoardType="Tapered")
predict(model.main, newdata=newdata)

##b
summary(model.12)

betahat<-coef(model.12)

k1<-c(0,0,0,0,0,1,0,0)
k2<-c(0,0,0,0,0,0,1,0)
k3<-c(0,0,0,0,0,0,0,1)

K<-cbind(k1,k2,k3)

q<-3
Wald<-(t(t(K)%*%betahat)%*%solve(t(K)%*%vcov(model.12)%*%K)%*%t(K)%*%betahat)/q
Wald
p.value<-pf(Wald, q, summary(model.12)$df[2], lower.tail = FALSE)
p.value

anova(model.main,model.12, test="F")

##c
newdata<-data.frame(WoodType=c("Oak","Cherry"), BoardType=c("Tapered","Stacked"))
predict(model.12, newdata=newdata)

x1<-cbind(c(1,1,0,0,0,0,0,0)) ### Cherry and Stacked
x2<-cbind(c(1,0,0,1,1,0,0,1)) ## Oak and Tapered
betahat<-cbind(coef(model.12))

pred<-(t(x2)-t(x1))%*%betahat
pred
sigma2<-sigma(model.12)^2

# predictive hypothesis testing
T<-pred/sqrt(sigma2*(2+(t(x2)-t(x1))%*%solve(vcov(model.12))%*%(x2-x1)))
T
p<-2*pt(abs(T),df=328, lower.tail = FALSE)
p

