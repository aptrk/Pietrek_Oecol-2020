#Author: Alejandro Pietrek
#Looking at density dependence and spread of Pheidole megacephala
#Paper accepted in Oecologia
#Last updated 27/1/2021

#Try to look at density dependence


library(dplyr)
library(tidyverse)
library(lme4)
library(nlme)
library(sjPlot)
library(effects)
library(truncnorm)

#Import data_bha in the txt version, convert a few variables

data_bha$transect <- as.factor(data_bha$transect)
data_bha$time_point <- as.factor(data_bha$time_point)


#Sumarize three baits per point, and remove three dates with rain: 12/2/2018, 13/7/2017,5/7/2017
#Also import data_bha from the text file


data_bha <- dplyr::select (data_bha, site:bha) %>%
  filter(date != "12/2/2018" & date != "13/7/2018" & date != "5/7/2018") %>%
  group_by(site, transect, survey, time_point, dist) %>%
    summarize(bha = sum(bha,  na.rm = TRUE)) 

#Get means and sd

by_time_point <- group_by(data_bha,transect,time_point,dist) %>%
  summarize(
    mean_bha = mean(bha, na.rm = TRUE),
    sd_bha = sd (bha, na.rm = TRUE))

# Plot how densities vary over time averaging all time_points

by_time_point2 <- group_by(by_time_point, dist, transect)%>%
  summarize(
    mean_bha = mean(mean_bha,na.rm = TRUE),
    sd_bha = sd (mean_bha, na.rm = TRUE)) %>%
  filter(transect != "9" & transect != "11")

by_time_point2$transect<- fct_recode(by_time_point2$transect, "5" = "7", "6" = "8") 



ggplot(aes(x = dist, y = mean_bha, group = transect), data = by_time_point2) +  
  geom_point(aes(shape = transect)) + geom_line() +
  theme_bw() +  # use a white background 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylab ("Mean ant density") +
  xlab ("Distance on the transect (m)") +
  scale_shape_discrete(name = "Transect") +
  theme(axis.title=element_text(size=14))



#Now plot all data in a single plot and run a non-linear mixed regression


all_data <- by_time_point %>%
           dplyr::select (time_point,dist,mean_bha) %>%
           spread (value = mean_bha, key = dist)


#Estimate lamdas

all_data$transect <- as.factor(all_data$transect)
all_data$time_point <- as.factor(all_data$time_point)

lamda <- all_data %>%group_by(transect) %>% 
  transmute_if(is.double, list(lamda = ~log((.)/lag(.))))


#Plot lamdas, find out if there is density dependence


Nt <- as.matrix(all_data[-nrow(all_data), 3:32] )# There is no estimate of lambda for the last row
lambda <- as.matrix(lamda[-1,2:31]) # There is no estimate of lamda for the first row

matplot(Nt, lambda)



# Arrange data to analyze in ggplot, comprehensive, and looks nicer!


lambda <- as.data.frame(lambda) %>%
  gather('0_lamda':'900_lamda', key= "Distance", value = "Lambda" )

Nt <- all_data[-nrow(all_data),] %>%
  gather('0':'900', key = "Distance", value = Nt)
  
data <-  bind_cols(Nt,lambda) 


#We need to tide up data and square lamdas of first survey to make them comparable (2 vs 4 months)
data <- data[,-5] %>% 
  filter(is.finite(Lambda))%>%  
  mutate(lambda = ifelse(time_point == "1", 2*Lambda, Lambda)) %>%
  select(-Lambda)%>%
  rename(Lambda = lambda)
  

#Look at density dependence

ggplot(data = data) +
  geom_point(mapping = aes(x = Nt, y = Lambda, color = transect )) +
  ylab("log(lambda)") +
  xlab(expression("Worker abundance"))

# DD no transects
tiff("Fig_2.tiff", width = 7, height = 5, units = 'in', res = 300)
ggplot(data = data) +
  geom_point(mapping = aes(x = Nt, y = Lambda, cex.axis = 1.5 )) +
  ylab("Observed local population growth") +
  xlab(expression("Worker abundance"))+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title=element_text(size=16))
dev.off()
  




#Look at temporal variability in pop. growth rates

ggplot(data = data) +
  geom_point(mapping = aes(x = time_point, y = Lambda, color = transect )) 

ggplot(data = data) +
  geom_boxplot(mapping = aes(x = time_point, y = Lambda))
  
# Check at temporal variability looking for a trend
#Try ro run a linear model to detect temporal trends

data$time_point <- as.numeric(data$time_point)
m1 <- lm(Lambda ~ time_point, data=data)
 
# Look at density dependence
#Correlation 

cor(data$Lambda,data$Nt) #-0.43
 

# Run  first a linear mixed effect model
#Fit the Ricker equation N(t+1) =  N(t)*exp{ r(1-N(t)/K) } by doing a linear regression of
# log( N(t+1)/N(t) ) against N(t), using the R function lm ("linear model");
# in the resulting regression equation y=b[1] + b[2]*x, r=b[1] and K=-b[1]/b[2]


m2 <- lmer(Lambda ~ Nt + (1|transect/Distance), data = data)
m3 <-lmer(Lambda ~ Nt + (1|transect/Distance), data = data, REML = TRUE)
m4 <- lm(Lambda ~ Nt, data = data)



ggplot(data = data) +
  geom_point(mapping = aes(x = Nt, y = Lambda, color = transect)) +
  geom_line(aes(y= predict(gplusr), x = Nt, color = transect, lty = transect), size = 1) +
  ylab("r obs")


#The non-linear regression equivalent


ricker <- nls(Lambda ~ r*(1-Nt/K), data = data, start = list(r = 0.5, K = 400))
ricker2 <- nlme(Lambda ~ r*(1-Nt/K) ,fixed= list(r~1, K~1), random = r + K ~1|transect / Distance, data = data, start = c(r = 0.5, K = 200))
ricker3 <- nlmer(Lambda ~ r*(1-Nt/K),fixed= list(r~1, K~1), random = r + K ~1|transect / Distance, data = data, start = c(r = 0.5, K = 200))



#Fit a theta logistic model using a non-linear regression

library(minpack.lm)


theta1 <- nls(Lambda ~ r*(1- (Nt/K))^theta ,data = data, start = list(r = 0.5, K = 400, theta = 1))
theta2 <- nlsLM(Lambda ~ r *(1- (Nt/K))^theta, start = list(r = 0.5, K= 300, theta = 1), data = data, lower= c(0,10,0.1))



#Attempt to fit an Allee effect, but it obviously won't work!


allee <- nls(Lambda~ r*-B*Nt*log(Nt/(A+Nt)), data = data, start = list(r = 0.5, B = 0.5, A = 10)) # not working, but apparently the maximum likelihood estimate for A should end up being zero, which interestingly reduces the equation down to the Ricker model (B = r/K).



# Trying to fit a Gompertz  discrete model as a simple linear model


gompertz <- lm(Lambda ~ log(Nt),data = data)



#Trying to fit a Gompertz including random effects

gompertz2 <- lmer(Lambda ~  log(Nt) + (1|transect / Distance), data = data,REML = FALSE)

#Under this model K is K = Exp[-rmax/b], 78.46 in this case




##########################################################################################

#Are trees and rain impacting population growth rates?

ndata <- filter(data,transect != "9" & transect != "11")

#Load Trees_transect and raindist first

Trees_transect$Transect <- as.character(Trees_transect$Transect)
Trees_transect$Point <- as.factor(Trees_transect$Point)


raindist$Transect <- as.factor(raindist$Transect)

data_t <- left_join(ndata, Trees_transect,  by=c("transect" = "Transect", "Distance" = "Point"))


data_r <- select(raindist,c(Transect, Rainfall, Tpoint)) %>%
  right_join(data_t ,  by=c("Transect" = "transect", "Tpoint" = "time_point")) %>%
   rename(Saplings =`<1m`, Trees = `>1m` )
   


#Looking at  basic plots

plot (Lambda ~Saplings, data = data_r)
plot (Lambda ~Trees, data = data_r)

# Checking on ggplot

ggplot(data = data_r) +
  geom_point(mapping = aes(x = Trees, y = Lambda)) +
  ylab("r obs") +
  xlab( "Number of trees")


ggplot(data = data_r) +
  geom_point(mapping = aes(x = Saplings, y = Lambda)) +
  ylab("r obs") +
  xlab( "Number of saplings")

cor(data_r$Trees, data_r$Saplings, use = "complete.obs")




# Modelling population growth rates
#Fit a Ricker and a gompertz with extrinsic drivers

rplus <- lm(Lambda ~ Nt + Trees+ Saplings+ Rainfall +
           Nt*Trees + Nt*Saplings + Nt * Rainfall , data = data_r)

gplus <- lm(Lambda ~ log(Nt) + Trees+ Saplings+ Rainfall + Nt +
              log(Nt)*Trees + log(Nt)*Saplings + log(Nt) * Rainfall , data = data_r)



#Adding random effects, 44 AIC units difference

rplusr <- lmer(Lambda ~ scale(Nt) + scale(Trees)+ scale(Saplings)+scale(Rainfall) +
              scale(Nt)* scale(Trees) + scale(Nt)*scale(Saplings) + scale(Nt) * scale(Rainfall) + (1|Transect),
             data = data_r, REML= FALSE)

gplusr <- lmer(Lambda ~ scale(log(Nt)) + scale(Trees)+ scale(Saplings)+ scale(Rainfall) +
              scale(log(Nt))*scale(Trees) + scale(log(Nt))*scale(Saplings) + scale(log(Nt)) * scale(Rainfall) + (1|Transect) , data = data_r, REML= FALSE)



#Looking a bit at estimates and interaction terms


plot_model(gplusr, type = "est")
plot(effect(term="Rainfall:Nt", mod= gplus, xlab= "Tree cover", 
     ylab=("Tree cover change"), main = NULL, factor.names = FALSE, font.lab = 2)) 


#Impact of rain on  population growth rates, build dataset


data_r <- data %>%
  mutate( Rain = case_when(
    (transect == 1|transect ==2|transect == 3|transect == 4) & time_point == 1 ~ 7,
    (transect == 1|transect ==2|transect == 3|transect == 4) & time_point == 2 ~ 216,
    (transect == 1|transect ==2|transect == 3|transect == 4) & time_point == 3 ~ 335,
    (transect == 1|transect ==2|transect == 3|transect == 4) & time_point == 4 ~ 92,
    (transect == 1|transect ==2|transect == 3|transect == 4) & time_point == 5 ~ 585,
    (transect == 7|transect ==8) & time_point == 2 ~ 173.5,
    (transect == 7|transect ==8) & time_point == 3 ~ 341,
    (transect == 7|transect ==8) & time_point == 4 ~ 118,
    (transect == 7|transect ==8) & time_point == 5 ~ 656))


ggplot(data = data_r) +
  geom_point(mapping = aes(x = Rain, y = Lambda)) +
  ylab("r obs") +
  xlab( "Rain")




#########################################################################################


# Is rainfall affecting rates of spread

#Load rainfall data



raindist$Transect <- as.factor(raindist$Transect)

# Give a try to saplings + trees


raindist <- mutate(raindist,
                   all = Saplings + Trees)
 

               
Cand.models2 <- list( )

Cand.models2[[1]] <- lmer(Distance ~ Rainfall + (1|Transect), data = raindist, REML = FALSE) #run a model
Cand.models2[[2]] <- lmer(Distance ~ Rainfall + Trees+ (1|Transect), data = raindist, REML = FALSE)
Cand.models2[[3]] <- lmer(Distance ~ Trees+ (1|Transect), data = raindist, REML = FALSE)
Cand.models2[[4]] <- lmer(Distance ~ Trees+ Saplings+ (1|Transect), data = raindist, REML = FALSE)
Cand.models2[[5]] <- lmer(Distance ~ Saplings+ (1|Transect), data = raindist, REML = FALSE)
Cand.models2[[6]]<- lmer(Distance ~ Saplings+ Rainfall+ (1|Transect), data = raindist, REML = FALSE)
Cand.models2[[7]]<- lmer(Distance ~  (1|Transect), data = raindist, REML = FALSE)


Modnames <- paste("mod", 1:length(Cand.models2), sep = " ")
#Generate AIC table
aictab(cand.set = Cand.models2, modnames = Modnames, sort = TRUE)
##round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = Cand.models2, modnames = Modnames, sort = TRUE),
      digits = 4, LL = TRUE)



#With transect as a fixed effect



m13 <- lm(Distance ~ Rainfall, data = raindist) 
m14 <- lm(Distance ~  Trees, data = raindist)
m15 <- lm(Distance ~ Trees + Rainfall, data = raindist)
m16 <- lm(Distance ~ Trees * Rainfall, data = raindist)
m17 <- lm(Distance ~ 1, data = raindist) 
m18 <- lm(Distance ~  Saplings, data = raindist)

ggplot(data = raindist) +
  geom_point(mapping = aes(x = Rainfall , y = Distance , color = Site), size = 2) + #Plot spread vs rain
  geom_line(aes(y= predict(m13), x = Rainfall), size = 1) +
  ylab("Distance moved from the front (m)") + xlab("Cumulative rainfall (mm)")


##########################################################################################

#Add uncertainty to dd population models with drivers


#Get means and sd


dat<- by_time_point



listnew <-
  replicate(n = 1000, expr = dat, simplify = F) %>% # Generate a list with the dataframes
  lapply(function(dat) {
    dat$mean_bha <-
      rtruncnorm(
        nrow(dat),
        a = 0,
        b = Inf,
        mean = dat$mean_bha,
        sd = dat$sd_bha
      )
    return(dat)
  })                                  # Generate new means from a N~(bha_mean, sd_bha)                                      



listdata <- listnew %>%
  lapply(function(x)
    x %>% select (transect:mean_bha)) %>%  # Ommit SD
  lapply(function(x)
    x %>% spread (value = 'mean_bha', key = 'dist')) #change disposition



listlamda <- listdata %>%
  lapply(function(x)
    x %>% group_by(transect)) %>%
  lapply(function(x)
    x %>% transmute_if(is.double, list(
      lamda = ~ log((.) / lag(.))))) %>%
  lapply(function(x)
    x[-1,]) %>%
  lapply(function(x)
    x %>% gather(`0_lamda`:`900_lamda`, key = "Distance", value = "Lambda"))  #Estimate population growth rates
          
                                   


Ntlist <- listdata %>% lapply(function(x)
  x [-42,]) %>% 
  lapply(function(x) 
    x %>% gather('0':'900', key = "Distance", value = Nt))
                   

ldata <-  mapply(cbind, Ntlist,listlamda, SIMPLIFY=FALSE) %>%
  lapply(function (x)
    x[,-(5:6)]) %>%
  lapply(function(x)
    x %>%
      filter(is.finite(Lambda))) %>%
  lapply(function(x)
    x %>%
      mutate(lambda = ifelse(time_point == "1",abs(Lambda)*Lambda, Lambda))) %>%
  lapply(function(x)
    x %>%
      select(-Lambda)) %>% 
  lapply(function(x)
    x %>%
      rename(Lambda = lambda)) %>%
  lapply(function(x)
    x %>%
      left_join(Trees_transect,  by=c("transect" = "Transect", "Distance" = "Point"))) #change transect and point to trabsect in that file
      
# Load rain data if we haven't before    

data_r <- select(raindist,c(Transect, Rainfall, Tpoint)) %>%
  rename("time_point" = "Tpoint", "transect"= "Transect" )
  
data_r$time_point <- as.factor(data_r$time_point)
  
final <- ldata %>%
  lapply(function(x)
    x %>%
      rename(Saplings =`<1m`, Trees = `>1m` )) %>%
  lapply(function(x)
    x %>%
      right_join(data_r,by=c( "transect", "time_point"))) 
  
  


lms <- lapply(final, function(a) 
  lmer(Lambda ~ scale(log(Nt)) + scale(Trees)+ scale(Saplings)+ scale(Rainfall)+
  scale(log(Nt))*scale(Trees) + scale(log(Nt))*scale(Saplings) + scale(log(Nt)) * scale(Rainfall) + (1|transect), data=a))

coef <- t(sapply(lms, fixef))
Intercept <- quantile(coef[,1], c(.025,.975))
logNt <- quantile(coef[,2], c(.025,.975))
Trees <- quantile(coef[,3], c(.025,.975))
Saplings <- quantile(coef[,4], c(.025,.975))
Rainfall <- quantile(coef[,5], c(.025,.975))
Ntrees <- quantile(coef[,6], c(.025,.975))
Ntsap <- quantile(coef[,7], c(.025,.975))
Ntrain <- quantile(coef[,8], c(.025,.975))

#Generate a dataframe to plot regression estimates


estimate <- c(mean(coef[,2]), mean(coef[,3]),mean(coef[,4]),mean(coef[,5]),mean(coef[,6]),
              mean(coef[,7]),mean(coef[,8]))
conf.low <- c(quantile(coef[,2],.025),quantile(coef[,3],.025),quantile(coef[,4],.025),quantile(coef[,5],.025),
            quantile(coef[,6],.025),quantile(coef[,7],.025),quantile(coef[,8],.025))

conf.high <- c(quantile(coef[,2],.975),quantile(coef[,3],.975),quantile(coef[,4],.975),quantile(coef[,5],.975),
            quantile(coef[,6],.975),quantile(coef[,7],.975),quantile(coef[,8],.975))

term <- c("log(Nt)","Tree cover", "Sapling cover", "Rain", "log(Nt) * Tree cover",
          "log(Nt) * Saplings cover", "log(Nt) * Rain")

df <- data.frame(term,estimate,conf.low, conf.high)

tiff("Fig_4.tiff", width = 7, height = 5, units = 'in', res = 300)
ggplot(data=df, aes(x= term, y= estimate, ymin= conf.low, ymax= conf.high)) +
  geom_pointrange() + 
  geom_hline(yintercept= 0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  ylab("Regression coefficients (95% CI)") +
  theme_bw() +  # use a white background 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title=element_text(size=14))
dev.off()



############################################################################

#Plotting a Gompertz with confidence bands extracted from a model without drivers

Ns=seq(1,400,1)
ymin= 2.33- 0.75 *log(Ns)
ymax= 3.20 - 0.549 *log(Ns)
d1= as.data.frame(cbind(Ns,ymin, ymax))


ggplot(data = data) +
  geom_point(mapping = aes(x = Nt, y = Lambda, color = transect )) +
  geom_line(aes(y= predict(gompertz2), x = Nt, color = transect, lty = transect), size = 1) +
  ylab("r obs") +
  geom_ribbon( data= d1, aes( x= Ns, ymin = ymin, ymax = ymax , color = NULL),fill = "red", alpha = .15)

########################################################################################################

#Add uncertainty with dd models without drivers


#Get means and sd


dat<- by_time_point %>%
  filter(transect != "9" & transect != "11")

Trees_transect$Transect <- as.character(Trees_transect$Transect)
Trees_transect$Point <- as.factor(Trees_transect$Point)



data_t <- left_join(ndata, Trees_transect,  by=c("transect" = "Transect", "Distance" = "Point"))
data_r <- select(raindist,c(Transect, Rainfall, Tpoint)) %>%
  right_join(data_t ,  by=c("Transect" = "transect", "Tpoint" = "time_point")) %>%
  rename(Saplings =`<1m`, Trees = `>1m` )




listnew <-
  replicate(n = 1000, expr = dat, simplify = F) %>% # Generate a list with the dataframes
  lapply(function(dat) {
    dat$mean_bha <-
      rtruncnorm(
        nrow(dat),
        a = 0,
        b = Inf,
        mean = dat$mean_bha,
        sd = dat$sd_bha
      )
    return(dat)
  })                                  # Generate new means from a N~(bha_mean, sd_bha)                                      



listdata <- listnew %>%
  lapply(function(x)
    x %>% select (transect:mean_bha)) %>%  # Ommit SD
  lapply(function(x)
    x %>% spread (value = 'mean_bha', key = 'dist')) #change disposition



listlamda <- listdata %>%
  lapply(function(x)
    x %>% group_by(transect)) %>%
  lapply(function(x)
    x %>% transmute_if(is.double, list(
      lamda = ~ log((.) / lag(.))))) %>%
  lapply(function(x)
    x[-1,]) %>%
  lapply(function(x)
    x %>% gather(`0_lamda`:`900_lamda`, key = "Distance", value = "Lambda"))  #Estimate population growth rates




Ntlist <- listdata %>% lapply(function(x)
  x [-42,]) %>% 
  lapply(function(x) 
    x %>% gather('0':'900', key = "Distance", value = Nt))


ldata <-  mapply(cbind, Ntlist,listlamda, SIMPLIFY=FALSE) %>%
  lapply(function (x)
    x[,-(5:6)]) %>%
  lapply(function(x)
    x %>%
      filter(is.finite(Lambda)))





lms <- lapply(ldata, function(a) lmer(Lambda~log(Nt) + (1|transect/Distance), data=a))
coef <- t(sapply(lms, fixef))
quantile(coef[,1], c(.025,.975))
quantile(coef[,2], c(.025,.975))


Ns=seq(1,400,1)
ymin= 2.33- 0.75 *log(Ns)
ymax= 3.20 - 0.549 *log(Ns)
d1= as.data.frame(cbind(Ns,ymin, ymax))


ggplot(data = data) +
  geom_point(mapping = aes(x = Nt, y = Lambda, color = transect )) +
  geom_line(aes(y= predict(gompertz2), x = Nt, color = transect, lty = transect), size = 1) +
  ylab("r obs") +
  geom_ribbon( data= d1, aes( x= Ns, ymin = ymin, ymax = ymax , color = NULL),fill = "red", alpha = .15)


###################################################################################################

#Estimate correlation between surveys, asked by a reviewer

cor_bha <- pivot_wider(data_bha,names_from = survey, values_from = bha) %>%
  filter(dist !="700" & dist != "800" & dist != "900")

cor(cor_bha$`1`,cor_bha$`2`, use = "complete.obs")  














