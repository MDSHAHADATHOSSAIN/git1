
library(boot)
library(MASS)
library(rcompanion)
library(ggplot2)
library(gridExtra)
library(FSA)
setwd("D:/Study/Summer Semester 2020/soil chemistry/r program")
setwd("D:/Study/Summer Semester 2020/soil chemistry/r program")
individual <- read.csv("Individual treatment.csv")
combined <- read.csv("Combined treatment Cmic Nmic.csv")

individual
combined

##checking the normality of data with graph
par(mfrow = c(2,1))
plot(density(individual$Cmic..mg.kg.))
plot(density(combined$Cmic..mg.kg.))
par(mfrow = c(2,1))
plot(density(individual$Cmic.Nmic..ratio.))
plot(density(combined$Cmic.Nmic..ratio.))

##To see the difference between median and mean
summary(individual$Cmic..mg.kg.)
summary(combined$Cmic..mg.kg.)
summary(individual$Cmic.Nmic..ratio.)
summary(combined$Cmic.Nmic..ratio.)
#significance difference between median and mean value which indicate that density of Cmic and Cmic/Nmic have long tail

#Transformation for Cmic (individual and combined)
#squareroot transformation
t <- individual$Cmic..mg.kg.
plotNormalHistogram(t)
t_sqrt <- sqrt(t)
summary(t_sqrt)
plotNormalHistogram(t_sqrt)
shapiro.test(t_sqrt)

#log transformation
t_log <- log10(t+1)
plotNormalHistogram(t_log)
summary(t_log)
shapiro.test(t_log)
# median and mean value are now closure

#Boxcox transformation
c = 1
t2 <- t+c #add  because the minimum value is zero
b <- boxcox(t2 ~ 1, lambda = seq(-6,6,0.3))
c <- data.frame(b$x, b$y)#this is to create the data frame from the fitting
c2 <- c[with(c,order(-c$b.y)),]#this is to sort the dataframe of c from the lowest value
#then extract the lambda
lambda <-c2[1,"b.x"]
#transform the original data
t.box <- (t^lambda - 1)/lambda

plotNormalHistogram(t.box)

qqnorm(t.box)
qqline(t.box,col="red", main = "Q-Q plot for individual treatment (Cmic)", cex.main = 1)
shapiro.test(t.box)

###Tukey transformation
t.tuk <- transformTukey(t)
plotNormalHistogram(t.tuk)
shapiro.test(t.tuk)

### we will take the exact value
qqnorm(individual$Cmic..mg.kg., main = "Q-Q plot for individual treatment (Cmic)", cex.main = 1)
qqline((individual$Cmic..mg.kg.), col= "blue")
par(mfrow = c(1,2))
boxplot(individual$Cmic..mg.kg. ~ individual$Treatments, xlab = "Treatments", ylab = "Cmic (mg/kg)", main = "Individual treatment", cex.main = 1)
boxplot(u.box ~ combined$Combined.Treatments, xlab = "Treatments", ylab = "Cmic (mg/kg)", main = "Combined treatment", cex.main = 1)


#############
#Transformation
#squareroot transformation (p value = 0.05977)
u <- combined$Cmic..mg.kg.
plotNormalHistogram(u)
u_sqrt <- sqrt(u)
summary(u_sqrt)
plotNormalHistogram(u_sqrt)
qqnorm(u_sqrt)
qqline(u_sqrt,col="red", main = "Q-Q plot for individual treatment (Cmic)", cex.main = 1)
shapiro.test(u_sqrt)

#log transformation (p value = 0.1122)
u_log <- log10(u+1)
plotNormalHistogram(u_log)
summary(u_log)
qqnorm(u_log)
qqline(u_log,col="red", main = "Q-Q plot for individual treatment (Cmic)", cex.main = 1)
shapiro.test(u_log)
# median and mean value are now closure

#Boxcox transformation (p value = 0.1168)
d = 1
u2 <- u+d #add  because the minimum value is zero
e <- boxcox(u2 ~ 1, lambda = seq(-6,6,0.3))
f <- data.frame(e$x, e$y)#this is to create the data frame from the fitting
f2 <- f[with(f,order(-f$e.y)),]#this is to sort the dataframe of c from the lowest value
#then extract the lambda
lambda <-f2[1,"e.x"]
#transform the original data
u.box <- (u^lambda - 1)/lambda ##cmic combined

plotNormalHistogram(u.box)

qqnorm(u.box)
qqline(u.box,col="red", main = "Q-Q plot for individual treatment (Cmic)", cex.main = 1)
shapiro.test(u.box)

#### we will use boxcox transformation
qqnorm(u.box,main = "Q-Q plot for combined treatment (Cmic)", cex.main = 1)
qqline((u.box), col= "blue")
boxplot(u.box ~ combined$Combined.Treatments, xlab = "Treatments", ylab = "Cmic (mg/kg)", main = "Combined treatment", cex.main = 1)
#For Boxcox transformationb p value is not reject null hypothesis



#################
#Transformation for Cmic/Nmic (individual and combined)
#squareroot transformation
v <- individual$Cmic.Nmic..ratio.
plotNormalHistogram(v)
v_sqrt <- sqrt(v)
summary(v_sqrt)
plotNormalHistogram(v_sqrt)
qqnorm(v_sqrt)
qqline(v_sqrt,col="red", main = "Q-Q plot for individual treatment (Cmic)", cex.main = 1)
shapiro.test(v_sqrt)

#log transformation 
v_log <- log10(v+1)
plotNormalHistogram(v_log)
summary(v_log)
qqnorm(v_log)
qqline(v_log,col="red", main = "Q-Q plot for individual treatment (Cmic)", cex.main = 1)
shapiro.test(v_log)
# median and mean value are now closure

#Boxcox transformation (p=0.09238)
a = 1
v2 <- v+a #add  because the minimum value is zero
g <- boxcox(v2 ~ 1, lambda = seq(-6,6,0.1))
h <- data.frame(g$x, g$y)#this is to create the data frame from the fitting
h2 <- h[with(h,order(-h$g.y)),]#this is to sort the dataframe of c from the lowest value
#then extract the lambda
lambda <-h2[1,"g.x"]
#transform the original data
v.box <- (v^lambda - 1)/lambda

plotNormalHistogram(v.box)

qqnorm(v.box)
qqline(v.box,col="red", main = "Q-Q plot for individual treatment (Cmic)", cex.main = 1)
shapiro.test(v.box)

#For Boxcox transformationb p value is not reject null hypothesis

####we will use boxcox transformation
qqnorm(v.box,main = "Q-Q plot for individual treatment (Cmic/Nmic)", cex.main = 1)
qqline((v.box), col= "blue")
####
par(mfrow = c(1,2))
boxplot(v.box ~ individual$Treatments, xlab = "Treatments", ylab = "Cmic/Nmic", main = "Individual treatment", cex.main = 1)
boxplot(w.box ~ combined$Combined.Treatments, xlab = "Treatments", ylab = "Cmic/Nmic", main = "Combined treatment", cex.main = 1)
par(mfrow=c(1,1))
#####
#squareroot transformation
w <- combined$Cmic.Nmic..ratio.
plotNormalHistogram(w)
w_sqrt <- sqrt(w)
summary(w_sqrt)
plotNormalHistogram(w_sqrt)
qqnorm(w_sqrt)
qqline(w_sqrt,col="red", main = "Q-Q plot for individual treatment (Cmic)", cex.main = 1)
shapiro.test(w_sqrt)

#log transformation 
w_log <- log10(w+1)
plotNormalHistogram(w_log)
summary(w_log)
qqnorm(w_log)
qqline(w_log,col="red", main = "Q-Q plot for individual treatment (Cmic)", cex.main = 1)
shapiro.test(w_log)
# median and mean value are now closure

#Boxcox transformation (0.6342)
z = 1
w2 <- w+z #add  because the minimum value is zero
i <- boxcox(w2 ~ 1, lambda = seq(-6,6,0.1))
j <- data.frame(i$x, i$y)#this is to create the data frame from the fitting
j2 <- j[with(j,order(-j$i.y)),]#this is to sort the dataframe of c from the lowest value
#then extract the lambda
lambda <-j2[1,"i.x"]
#transform the original data
w.box <- (w^lambda - 1)/lambda

plotNormalHistogram(w.box)

qqnorm(w.box)
qqline(w.box,col="red", main = "Q-Q plot for individual treatment (Cmic)", cex.main = 1)
shapiro.test(w.box)

#For Boxcox transformationb p value is not reject null hypothesis

#### we will use boxcox transformation
qqnorm(w.box,main = "Q-Q plot for combined treatment (Cmic/Nmic)", cex.main = 1)
qqline((w.box), col= "blue")
boxplot(w.box ~ combined$Combined.Treatments, xlab = "Treatments", ylab = "Cmic/Nmic", main = "Combined treatment", cex.main = 1)


#############
###
combined <- cbind(combined, u.box, w.box)
individual <- cbind(individual, v.box)
  
A <- ggplot(individual, aes(x= Treatments, y= Cmic..mg.kg., fill = Treatments)) + 
  geom_boxplot() +
  labs(x = "Treatments", y = "Cmic..mg.kg.", title = "individual") +
  theme_classic() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
B <- ggplot(combined, aes(x= Combined.Treatments, y= u.box, fill = Combined.Treatments)) + 
  geom_boxplot() +
  labs(x = "Treatments", y = "Cmic..mg.kg.(box-cox transformation)", title = "Combined") +
  theme_classic() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

grid.arrange(arrangeGrob(A, B, ncol = 2), nrow = 2, heights = c(10, 1))

par(mfrow = c(1,2))
C <- ggplot(individual, aes(x= Treatments, y= Nmic..µg..g.soil..1., fill = Treatments)) + 
  geom_boxplot() +
  labs(x = "Treatments", y = "Nmic..µg..g.soil..1.", title = "individual") +
  theme_classic() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
D <- ggplot(combined, aes(x= Combined.Treatments, y= Nmic..µg..g.soil..1., fill = Combined.Treatments)) + 
  geom_boxplot() +
  labs(x = "Treatments", y = "Nmic..µg..g.soil..1.", title = "Combined") +
  theme_classic() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

grid.arrange(arrangeGrob(C, D, ncol = 2), nrow = 2, heights = c(10, 1))


E <- ggplot(individual, aes(x= Treatments, y= v.box, fill = Treatments)) + 
  geom_boxplot() +
  labs(x = "Treatments", y = "Cmic.Nmic..ratio.(box-cox transformation)", title = "individual") +
  theme_classic() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
F <- ggplot(combined, aes(x= Combined.Treatments, y= w.box, fill = Combined.Treatments)) + 
  geom_boxplot() +
  labs(x = "Treatments", y = "Cmic.Nmic..ratio.(box-cox transformation)", title = "Combined") +
  theme_classic() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

grid.arrange(arrangeGrob(E, F, ncol = 2), nrow = 2, heights = c(10, 1))


#### for ANOVA
###Cmic Nmic and Nmic/Cmic ratio individual
bartlett.test(individual$Cmic..mg.kg.,individual$Treatments)
model11 <- kruskal.test(individual$Cmic..mg.kg.,individual$Treatments)
DT = dunnTest(individual$Cmic..mg.kg.,individual$Treatments,
              method="bh")      # Adjusts p-values for multiple comparisons;
# See ?dunnTest for options

DT


model2 <- aov(individual$Nmic..µg..g.soil..1. ~ individual$Treatments)
summary(model2)
TukeyHSD(model2)
t.test(individual$Nmic..µg..g.soil..1.)

model3 <- aov(v.box ~ individual$Treatments)
summary(model3)
TukeyHSD(model3)
t.test(individual$Cmic.Nmic..ratio.)


###combined
model4 <- aov(u.box ~ combined$Combined.Treatments)
summary(model4)
TukeyHSD(model4)
t.test(combined$Cmic..mg.kg.)

model5 <- aov(combined$Nmic..µg..g.soil..1. ~ combined$Combined.Treatments)
summary(model5)
TukeyHSD(model5)
t.test(combined$Nmic..µg..g.soil..1.)

model6 <- aov(w.box ~ combined$Combined.Treatments)
summary(model6)
TukeyHSD(model6)
t.test(combined$Cmic.Nmic..ratio.)


####Q-Q plot
#Cmic
par(mfrow = c(1,2))
qqnorm(t, main = "Q-Q plot for individual treatment (Cmic)", cex.main = 1)
qqline((t), col= "blue")
qqnorm(u.box, main = "Q-Q plot for combined treatment (Cmic)", cex.main = 1)
qqline((u.box), col= "blue")


#Nmic
par(mfrow = c(1,2))
qqnorm(individual$Nmic..µg..g.soil..1., main = "Q-Q plot for individual treatment (Nmic)", cex.main = 1)
qqline((individual$Nmic..µg..g.soil..1.), col= "blue")
qqnorm(combined$Nmic..µg..g.soil..1., main = "Q-Q plot for combined treatment (Nmic)", cex.main = 1)
qqline((combined$Nmic..µg..g.soil..1.), col= "blue")

#Cmic Nmic ratio
par(mfrow = c(1,2))
qqnorm(v.box, main = "Q-Q plot for individual treatment (Cmic/Nmic)", cex.main = 1)
qqline((v.box), col= "blue")
shapiro.test(v.box)
qqnorm(w.box, main = "Q-Q plot for combined treatment (Nmic)", cex.main = 1)
qqline((w.box), col= "blue")
shapiro.test(w.box)
par(oma=c(4,0,0,0))
mtext("Dider", side=1, outer=TRUE,adj=1, cex=0.6)

####Q-Q plot combined
par(mfrow = c(1,2))
qqnorm(u.box, main = "Q-Q plot for combined treatment (Cmic)", cex.main = 1)
qqline((u.box), col= "blue")
qqnorm(combined$Nmic..µg..g.soil..1., main = "Q-Q plot for combined treatment (Nmic)", cex.main = 1)
qqline((combined$Nmic..µg..g.soil..1.), col= "blue")
par(oma=c(4,0,0,0))
mtext("Dider", side=1, outer=TRUE,adj=1, cex=0.6)
qqnorm(w.box, main = "Q-Q plot for combined treatment Cmic/Nmic", cex.main = 1)
qqline((w.box), col= "blue")
par(oma=c(4,0,0,0))
mtext("Dider", side=1, outer=TRUE,adj=1, cex=0.6)
##### Boxplot
#Cmic
par(mfrow = c(1,2))
boxplot(individual$Cmic..mg.kg. ~ individual$Treatments, xlab = "Treatments", ylab = "Cmic (mg/kg)", main = "Individual treatment", cex.main = 1)
boxplot(u.box ~ combined$Combined.Treatments, xlab = "Treatments", ylab = "Cmic (mg/kg)", main = "Combined treatment (Box-cox transformation)", cex.main = 1)
par(oma=c(4,0,0,0))
mtext("Dider", side=1, outer=TRUE,adj=1, cex=0.6)

#Nmic
par(mfrow = c(1,2))
boxplot(individual$Nmic..µg..g.soil..1. ~ individual$Treatments, xlab = "Treatments", ylab = "Nmic (µg/(g.soil))", main = "Individual treatment", cex.main = 1)
boxplot(combined$Nmic..µg..g.soil..1. ~ combined$Combined.Treatments, xlab = "Treatments", ylab = "Nmic (µg/(g.soil))", main = "Combined treatment", cex.main = 1)
par(oma=c(4,0,0,0))
mtext("Dider", side=1, outer=TRUE,adj=1, cex=0.6)

#Cmic and Nmic ratio
par(mfrow = c(1,2))
boxplot(v.box ~ individual$Treatments, xlab = "Treatments", ylab = "Cmic/Nmic", main = "Individual treatment (Box-cox transformation)", cex.main = 1)
boxplot(w.box ~ combined$Combined.Treatments, xlab = "Treatments", ylab = "Cmic/Nmic", main = "Combined treatment (Box-cox transformation)", cex.main = 1)
par(oma=c(4,0,0,0))
mtext("Dider", side=1, outer=TRUE,adj=1, cex=0.6)


####Boxplot combined
par(mfrow = c(1,2))
boxplot(u.box ~ combined$Combined.Treatments, xlab = "Treatments", ylab = "Cmic (mg/kg)", main = "Combined treatment (Box-cox transformation)", cex.main = 1)
boxplot(combined$Nmic..µg..g.soil..1. ~ combined$Combined.Treatments, xlab = "Treatments", ylab = "Nmic (µg/(g.soil))", main = "Combined treatment", cex.main = 1)
par(oma=c(4,0,0,0))
mtext("Dider", side=1, outer=TRUE,adj=1, cex=0.6)

boxplot(w.box ~ combined$Combined.Treatments, xlab = "Treatments", ylab = "Cmic/Nmic", main = "Combined treatment (Box-cox transformation)", cex.main = 1)
par(oma=c(4,0,0,0))
mtext("Dider", side=1, outer=TRUE,adj=1, cex=0.6)
####
setwd("D:/Study/Summer Semester 2020/soil chemistry/r program")

mydat <- read.csv("barplot.csv")
head(mydat)
summary(mydat)
str(mydat)
barplot(mydat, ylim = c(0,800), names = mydat$Combined.Treatments, xlab = "Combined Treatments", ylab = "Cmic", main = "Cmic value for individual treatment", col = c(2,5))

library(agricolae)
par(mfrow = c(1,2), cex = 0.7)
require(agricolae)
treatments <- individual$Cmic..mg.kg. 
bar.treatment <- bar.err(x = kruskal(individual$Cmic..mg.kg.),
                         variation = "SE",
                         horiz = FALSE,
                         ylim = c(0,800),
                         names = individual$Treatments,
                         main = "Cmic for individual",
                         ylab = "Cmic (mg/kg)",
                         space = 0.5,
                         col = c(2,7))



##### VIF analysis
library(faraway)
###Cmic
setwd("D:/Study/Summer Semester 2020/soil chemistry/r program")
biomassvariables <- read.csv("Biomass with variables.csv")
head(biomassvariables)
mydata <- data.frame(biomassvariables[,-1])
head(mydata)
####Checking colinearity
round(cor(mydata),3)

####Check coefficient
Cmicmodel <- lm(Cmic..mg.kg.~., mydata)
Cmicmodel
summary(Cmicmodel)
vif(Cmicmodel)

###Nmic
setwd("D:/Study/Summer Semester 2020/soil chemistry/r program")
biomassvariables <- read.csv("Biomass with variables.csv")
head(biomassvariables)
mydataN <- data.frame(biomassvariables[,-1])
head(mydataN)
####Checking colinearity
round(cor(mydataN),3)

####Check coefficient
Nmicmodel <- lm(Nmic..µg..g.soil..1.~., mydata)
Nmicmodel
summary(Nmicmodel)
vif(Nmicmodel)
vif <- vif(Nmicmodel)

####Cmic Nmic ratio
setwd("D:/Study/Summer Semester 2020/soil chemistry/r program")
biomassvariables <- read.csv("Biomass with variables.csv")
head(biomassvariables)
mydataCN <- data.frame(biomassvariables[,-1])
head(mydataCN)
####Checking colinearity
round(cor(mydataCN),3)

####Check coefficient
CNmicmodel <- lm(Cmic.Nmic..ratio.~., mydata)
CNmicmodel
summary(CNmicmodel)
vif(CNmicmodel)
vif <- vif(CNmicmodel)



####
biomass <- read.csv("Biomass1.csv")
plot(biomass)


#calculate the correlation matrix
mycorr = cor(biomass)
mycorr
plot(cor(biomass))
plotNormalHistogram(mycorr)

# Eigensystem Analysis
eigen(cor(biomass))$values

#Ratio of maximum to minimum eigenvalues
max(eigen(cor(biomass))$values)/min(eigen(cor(biomass))$values)
kappa(cor(biomass), exact = TRUE)


####boxplot
setwd("D:/Study/Summer Semester 2020/soil chemistry/r program")
biomassvariables <- read.csv("Biomass with variables.csv")
head(biomassvariables)
par(mfrow = c(1,2))
boxplot(biomassvariables$pH.H2O. ~ biomassvariables$Combined.Treatments, xlab = "Treatments", ylab = "pH(H2O)", main = "pH for combined treatment", cex.main = 1)
boxplot(biomassvariables$WC.g.g. ~ biomassvariables$Combined.Treatments, xlab = "Treatments", ylab = "WC(g/g)", main = "Water content for combined treatment", cex.main = 1)
par(oma=c(4,0,0,0))
mtext("Shahadat", side=1, outer=TRUE,adj=1, cex=0.6)

par(mfrow = c(1,2))
boxplot(biomassvariables$Nitrate..mg..gsoil.. ~ biomassvariables$Combined.Treatments, xlab = "Treatments", ylab = "Nitrate ion (mg/(gsoil))", main = "Nitrate ion  for combined treatment", cex.main = 1)
boxplot(biomassvariables$NH4.mg..gsoil.. ~ biomassvariables$Combined.Treatments, xlab = "Treatments", ylab = "Ammonium ion (mg/(gsoil))", main = "NH4 ion for combined treatment", cex.main = 1)
par(oma=c(4,0,0,0))
mtext("Shahadat", side=1, outer=TRUE,adj=1, cex=0.6)

plot(density(biomassvariables$pH.H2O.))
shapiro.test(biomassvariables$pH.H2O.)
shapiro.test(biomassvariables$WC.g.g.)
shapiro.test(biomassvariables$Nitrate..mg..gsoil..)
shapiro.test(biomassvariables$NH4.mg..gsoil..)

trans