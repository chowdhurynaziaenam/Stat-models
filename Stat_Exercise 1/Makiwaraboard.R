data<-read.table("makiwaraboard.txt", sep="\t", dec=".", header=TRUE)
attach(data)

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

