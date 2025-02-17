#Question_1

data<-read.table("ozone.txt", header=TRUE, sep="\t", dec=".")
attach(data)

##a
normal.identity<-glm(ozone~rad+temp+wind, family=gaussian(link="identity"), data=data)
summary(normal.identity)
normal.inverse<-glm(ozone~rad+temp+wind, family=gaussian(link="inverse"), data=data)
summary(normal.inverse)
normal.log<-glm(ozone~rad+temp+wind, family=gaussian(link="log"), data=data)
summary(normal.log)
normal.exponential<-glm(ozone~log(rad)+log(temp)+log(wind), family=gaussian(link="log"), data=data)
summary(normal.exponential)

AIC(normal.identity) ## 998.6276
AIC(normal.inverse) ## 999.0104
AIC(normal.log) ## 972.5169
AIC(normal.exponential) ## 974.9654

#the best model = model.log

##b

gamma.log <- glm(ozone~rad+temp+wind, family = Gamma(link = "log"), data=data)
summary(gamma.log)

inversegaussian.log <- glm(ozone~rad+temp+wind, family = inverse.gaussian(link = "log"), data=data)
summary(inversegaussian.log)

plot(fitted(normal.log, type="response"), residuals(normal.log, type="pearson")^2)
plot(fitted(gamma.log, type="response"), residuals(gamma.log, type="pearson")^2)
plot(fitted(inversegaussian.log, type="response"), residuals(inversegaussian.log, type="pearson")^2)

pearsonresidual.o1 <- residuals(normal.log, type = "pearson")^2
mu.1 <- fitted(normal.log, type= "response")
model.pearsonN <- lm(pearsonresidual.o1~mu.1)
summary(model.pearsonN)

pearsonresidual.o2 <- residuals(gamma.log, type = "pearson")^2
mu.2 <- fitted(gamma.log, type= "response")
model.pearsonG <- lm(pearsonresidual.o2~mu.2)
summary(model.pearsonG)

pearsonresidual.o3 <- residuals(inversegaussian.log, type = "pearson")^2
mu.3 <- fitted(inversegaussian.log, type= "response")
model.pearsonIG <- lm(pearsonresidual.o3~mu.3)
summary(model.pearsonIG)

# here gamma.log model has minimum negative effect to the response variable 
# with t value -2.023, where (-) represents negative effect and
# 2.023 represents the magnitude.

##c
gammaH0 <- glm(ozone~rad, family = Gamma(link = "log"), data = data)
testhypothesis <- anova(gammaH0, gamma.log, test = "F")
testhypothesis$F[2] ## 69.28196

##gamma.log model provides better fit than that of gammaH0 model



#Question_2

data<-read.table("weld.txt", sep="\t", dec=".", header=TRUE)
attach(data)

##a

interaction.plot(Drying, Material, Strength)

#Checking Normal Distribution
normal.1<-glm(Strength~Drying+Material, family=gaussian(link="identity"), data=data)
normal.2<-glm(Strength~Drying+Material, family=gaussian(link="log"), data=data)
normal.3<-glm(Strength~Drying+Material, family=gaussian(link="inverse"), data=data)
normal.4<-glm(Strength~Drying*Material, family=gaussian(link="inverse"), data=data)

AIC(normal.1) # 30.68635
AIC(normal.2) # 30.46339
AIC(normal.3) # 30.44052
AIC(normal.4) # 32.42567

#Checking gamma

gamma.1<-glm(Strength~Drying+Material, family=Gamma(link="identity"), data=data)
gamma.2<-glm(Strength~Drying+Material, family=Gamma(link="log"), data=data)
gamma.3<-glm(Strength~Drying+Material, family=Gamma(link="inverse"), data=data)

AIC(gamma.1) # 29.861
AIC(gamma.2) # 29.6274
AIC(gamma.3) # 29.6033

#Checking inverse gaussian
inversegaussian.1<-glm(Strength~Drying+Material, family=inverse.gaussian(link="identity"), data=data)
inversegaussian.2<-glm(Strength~Drying+Material, family=inverse.gaussian(link="log"), data=data)
inversegaussian.3<-glm(Strength~Drying+Material, family=inverse.gaussian(link="inverse"), data=data)

AIC(inversegaussian.1) # 29.4747 
AIC(inversegaussian.2) # 29.2367
AIC(inversegaussian.3) # 29.21215


inversegaussianH0<-glm(Strength~Drying, family=inverse.gaussian(link="inverse"), data=data)
anova(inversegaussianH0, inversegaussian.3, test="F")


newdata<-data.frame(Drying=c(0,1), Material=c(0,1))
newdata

pred<-predict(inversegaussian.3, newdata=newdata, type="response")
pred

x1f<-cbind(c(1,0,1)) 
x2f<-cbind(c(1,0,1)) 
Xf<-t(cbind(x1f,x2f))

k<-cbind(c(-1,1))

phi<-summary(inversegaussian.3)$dispersion
Var.Y1f<-phi*(pred[1]^3)
Var.Y2f<-phi*(pred[2]^3)

D.f<-diag(pred)   

Var.ef<-Var.Y1f+Var.Y2f+t(k)%*%D.f%*%Xf%*%vcov(inversegaussian.3)%*%t(Xf)%*%D.f%*%k


Q<-(pred[2]-pred[1])/sqrt(Var.ef)
d<-2*pnorm(abs(Q), lower.tail = FALSE)
d  # 0.1946375

##b

# Loading the multcomp library
library(multcomp)

# Creating a grid of all possible combinations of Drying and Material
newdata <- expand.grid(Drying=c(0,1), Material=c(0,1))
# Converting into a model matrix
X <- model.matrix(~Drying+Material, data=newdata)
# Creating the contrast matrix for testing all possible differences
K <- cbind(t(X[-1,]) - X[1,],
           t(X[-(1:2),]) - X[2,])

# Performing the linear hypothesis test
pairwise <- glht(inversegaussian.3, linfct = t(K))


data.frame(t(K),summary(pairwise)$test$pvalues)


summary(pairwise, Chisqtest())

##c


#Checking Normal Distribution
normal.1<-glm(Strength ~ Drying + Material + Thickness + Angle + Opening + Preheating, family=gaussian(link="identity"), data=data)
normal.2<-glm(Strength ~ Drying + Material + Thickness + Angle + Opening + Preheating, family=gaussian(link="log"), data=data)
normal.3<-glm(Strength ~ Drying + Material + Thickness + Angle + Opening + Preheating, family=gaussian(link="inverse"), data=data)
normal.4<-glm(Strength ~ Drying + Material + Thickness + Angle + Opening + Preheating, family=gaussian(link="inverse"), data=data)

AIC(normal.1) # 32.28658
AIC(normal.2) # 32.02075
AIC(normal.3) # 31.76837
AIC(normal.4) # 31.76837

#Checking gamma

gamma.1<-glm(Strength ~ Drying + Material + Thickness + Angle + Opening + Preheating, family=Gamma(link="identity"), data=data)
gamma.2<-glm(Strength ~ Drying + Material + Thickness + Angle + Opening + Preheating, family=Gamma(link="log"), data=data)
gamma.3<-glm(Strength ~ Drying + Material + Thickness + Angle + Opening + Preheating, family=Gamma(link="inverse"), data=data)

AIC(gamma.1) # 31.55634
AIC(gamma.2) # 31.307
AIC(gamma.3) # 31.07271

#Checking inverse gaussian
inversegaussian.1<-glm(Strength ~ Drying + Material + Thickness + Angle + Opening + Preheating, family=inverse.gaussian(link="identity"), data=data)
inversegaussian.2<-glm(Strength ~ Drying + Material + Thickness + Angle + Opening + Preheating, family=inverse.gaussian(link="log"), data=data)
inversegaussian.3<-glm(Strength ~ Drying + Material + Thickness + Angle + Opening + Preheating, family=inverse.gaussian(link="inverse"), data=data)

AIC(inversegaussian.1) # 31.19959
AIC(inversegaussian.2) # 30.95888
AIC(inversegaussian.3) # 30.73173


gammaH0<-glm(Strength~Drying, family=Gamma(link="inverse"), data=data)
anova(gammaH0, gamma.3, test="F")


inversegaussianH0<-glm(Strength~Drying, family=inverse.gaussian(link="inverse"), data=data)
anova(inversegaussianH0, inversegaussian.3, test="F")

#From AIC and anova test it seems that we can take either gamma distribution with 
#inverse as a link function or inverse-gaussian distribution
#inverse as a link function as their AIC value is lower and F-value is also similar

#Taking Inverse-gaussian as a model


#Calculating the fitted value

fitted_value <- predict(inversegaussian.3, newdata = data[1, ])

fitted_value # 0.02286562

 
