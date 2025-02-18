library(tidyverse)
library(plotly)

data<-read.table("paper.txt", sep="\t", dec=".", header=TRUE)
attach(data)

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


