library(Epi)
library(survival)
library(flexsurv)
library(GlobalDeviance)
library(survminer)
data(thoro)


# # #####################################################
# # clean up the data
# # #####################################################
outid = which(thoro$cause == 15 | thoro$cause ==16)
thoro = thoro[-outid,]

id = thoro$id
Age = as.numeric(thoro$exitdat-thoro$birthdat)/365
Sex = thoro$sex
dose= (thoro$volume==0)*0 + (thoro$volume>0 & thoro$volume <= 20)*1 + (thoro$volume>20)*2
exitstat = (thoro$exitstat ==2 | thoro$exitstat ==3) * 0 + (thoro$exitstat ==1) *1

newdat = data.frame("id"=id, "Age" = Age, "Sex" = Sex, "Dose" = dose, "censor" = exitstat, "cause" = thoro$cause)
##1:male, dead
##2:female, cancer
##0:censored

id40 = which(newdat$Age >= 40)  ####cut off 211
newdat = newdat[id40,]


head(newdat)


# # #####################################################
# # fit the model
# # #####################################################
fitcox = coxph(Surv(Age,censor)~as.factor(Sex) + as.factor(Dose) + as.factor(Sex)*as.factor(Dose), data = newdat, method = "efron")
plot(survfit(Surv(Age,censor)~as.factor(Sex) + as.factor(Dose), data = newdat), col = 1:6, xlim = c(40, 120))
ggsurvplot(survfit(Surv(Age,censor)~as.factor(Sex) + as.factor(Dose), data = newdat), xlim = c(40, 120),
           legend.labs = c(paste("Male", 0:2), paste("Female", 0:2)))

fit = flexsurvreg(Surv(Age,censor)~as.factor(Sex) + as.factor(Dose) + as.factor(Sex)*as.factor(Dose), data = newdat, dist = "weibull")
fit2 = flexsurvreg(Surv(Age,censor)~as.factor(Sex) + as.factor(Dose) + as.factor(Sex)*as.factor(Dose), data = newdat, dist = "gamma")
fit$AIC
fit2$AIC

coef = coef(fit)

######male && dose0
m0 = integrate(l<-function(x) {1-pweibull(x, shape = exp(coef[1]), scale = exp(coef[2]))}, 0, Inf)$value

######male && dose1
m1 = integrate(l<-function(x) {1-pweibull(x, shape = exp(coef[1]), scale = exp(coef[2]+coef[4]))}, 0, Inf)$value

######male && dose2
m2 = integrate(l<-function(x) {1-pweibull(x, shape = exp(coef[1]), scale = exp(coef[2]+coef[5]))}, 0, Inf)$value

######female && dose0
f0 = integrate(l<-function(x) {1-pweibull(x, shape = exp(coef[1]), scale = exp(coef[2]+coef[3]))}, 0, Inf)$value

######female && dose1
f1 = integrate(l<-function(x) {1-pweibull(x, shape = exp(coef[1]), scale = exp(coef[2]+coef[3]+ coef[4]+coef[6]))}, 0, Inf)$value

######male && dose2
f2 = integrate(l<-function(x) {1-pweibull(x, shape = exp(coef[1]), scale = exp(coef[2]+coef[3]+coef[5]+coef[7]))}, 0, Inf)$value



cov = fit$cov[-1,-1]
sdm0 = sqrt(c(1,0,0,0,0,0) %*% cov %*% c(1,0,0,0,0,0))
sdm1 = sqrt(c(1,0,1,0,0,0) %*% cov %*% c(1,0,1,0,0,0))
sdm2 = sqrt(c(1,0,0,1,0,0) %*% cov %*% c(1,0,0,1,0,0))
sdf0 = sqrt(c(1,1,0,0,0,0) %*% cov %*% c(1,1,0,0,0,0))
sdf1 = sqrt(c(1,1,1,0,1,0) %*% cov %*% c(1,1,1,0,1,0))
sdf2 = sqrt(c(1,1,0,1,0,1) %*% cov %*% c(1,1,0,1,0,1))

m0l = integrate(l<-function(x) {1-pweibull(x, shape = exp(coef[1]), scale = exp(coef[2] - 1.96*sdm0))}, 0, Inf)$value
m0u = integrate(l<-function(x) {1-pweibull(x, shape = exp(coef[1]), scale = exp(coef[2] + 1.96*sdm0))}, 0, Inf)$value

m1l = integrate(l<-function(x) {1-pweibull(x, shape = exp(coef[1]), scale = exp(coef[2] + coef[4] - 1.96*sdm1))}, 0, Inf)$value
m1u = integrate(l<-function(x) {1-pweibull(x, shape = exp(coef[1]), scale = exp(coef[2] + coef[4] + 1.96*sdm1))}, 0, Inf)$value

m2l = integrate(l<-function(x) {1-pweibull(x, shape = exp(coef[1]), scale = exp(coef[2] + coef[5] - 1.96*sdm2))}, 0, Inf)$value
m2u = integrate(l<-function(x) {1-pweibull(x, shape = exp(coef[1]), scale = exp(coef[2] + coef[5] + 1.96*sdm2))}, 0, Inf)$value

f0l = integrate(l<-function(x) {1-pweibull(x, shape = exp(coef[1]), scale = exp(coef[2] + coef[3] - 1.96*sdf0))}, 0, Inf)$value
f0u = integrate(l<-function(x) {1-pweibull(x, shape = exp(coef[1]), scale = exp(coef[2] + coef[3] + 1.96*sdf0))}, 0, Inf)$value

f1l = integrate(l<-function(x) {1-pweibull(x, shape = exp(coef[1]), scale = exp(coef[2] + coef[3] + coef[4] - 1.96*sdf1))}, 0, Inf)$value
f1u = integrate(l<-function(x) {1-pweibull(x, shape = exp(coef[1]), scale = exp(coef[2] + coef[3] + coef[4] + 1.96*sdf1))}, 0, Inf)$value

f2l = integrate(l<-function(x) {1-pweibull(x, shape = exp(coef[1]), scale = exp(coef[2] + coef[3] + coef[5] - 1.96*sdf2))}, 0, Inf)$value
f2u = integrate(l<-function(x) {1-pweibull(x, shape = exp(coef[1]), scale = exp(coef[2] + coef[3] + coef[5] + 1.96*sdf2))}, 0, Inf)$value

Male0 = c(m0, exp(sdm0), m0l, m0u)
Male1 = c(m1, exp(sdm1), m1l, m1u)
Male2 = c(m2, exp(sdm2), m2l, m2u)
Female0 = c(f0, exp(sdf0), f0l, f0u)
Female1 = c(f1, exp(sdf1), f1l, f1u)
Female2 = c(f2, exp(sdf2), f2l, f2u)

T <- rbind(Male0, Male1, Male2, Female0, Female1, Female2)
colnames(T) <- c("Mean Survival Time", "Std Error", "Lower 95% CI", "Upper 95% CI")
rownames(T) <- c(paste("Male, Dose ", 0:2), paste("Female, Dose ", 0:2))

library(xtable)
xtable(T)



# # #####################################################
# # problem 2
# # #####################################################

head(newdat)

cancer = rep(0, dim(newdat)[1])
cancer[which(newdat$cause == 2)] = 1

###########treated vs not treated
contro_id = which(newdat$Dose == 0)

contro = newdat[contro_id,]
canp1 = length(which(contro$cause == 2))/dim(contro)[1]



trt = newdat[-contro_id,]
canp2 = length(which(trt$cause == 2))/dim(trt)[1]



prop.test(c(length(which(contro$cause == 2)),length(which(trt$cause == 2))),c(dim(contro)[1],dim(trt)[1]) )


grp = rep(0, dim(newdat)[1])
grp[-contro_id] = 1
fit2 = glm(cancer~1+grp, family = "binomial")
table(predict(fit2, newdata = data.frame(rep(1, dim(newdat)[1]), grp), type="response"))





















