library(readxl) #import data
library(MASS) #polr
library(ggplot2)
library(ggpubr)#ggarrange
library(marginaleffects)
library(dplyr)
speech <- read_excel("C:/Users/speech.xlsx")

#Create dummy variables
speech$Male=ifelse(speech$sex=='M',1,0)
speech$Gene_other=ifelse(is.na(speech$mutation),NA,ifelse(speech$mutation=='OTHER',1,0))
speech$Gene_PCDH19=ifelse(is.na(speech$mutation),NA,ifelse(speech$mutation=='PCDH19',1,0))
speech$Gene_SCN1A=ifelse(is.na(speech$mutation),NA,ifelse(speech$mutation=='SCN1A',1,0))


#Discard unrealistic ages
speech=subset(speech,age<100)

#Create the lag#####
speech <- speech[order(speech$id),]

lag_speech=function(i,l){
  if((i-l)<0){
    return('NA')
  }
  else{
    return(speech$Speech2[i-l])
  }
}

speech$l1=rep(NA,nrow(speech))
for(v in 2:nrow(speech)){
  speech$l1[v]=lag_speech(v,1)
  if(speech$time[v]==0){
    speech$l1[v]=NA
  }
}


###fit transition model#####
fulldata=(speech[!is.na(speech$l1)&!is.na(speech$Speech2)&!is.na(speech$baseline_age)&!is.na(speech$Male)&!is.na(speech$Gene_other),])
Markov=polr(as.factor(Speech2)~baseline_age+Male+Gene_other+Gene_PCDH19+time+l1+
              l1*baseline_age+l1*Male+l1*time,data=fulldata,
     Hess = TRUE, model = TRUE,
     method = c("probit"))

summary(Markov)

#get predicted values
Markov$fitted.values
predicted=cbind(fulldata,Markov$fitted.values)

#method with confidence intervals
nd=expand.grid(
  baseline_age=0.3,Male=c(0,1),Gene_other=0,Gene_PCDH19=0,time=seq(from=0,to=18,by=0.1),l1=c("0.Absent","1.Poor","2.Normal"))
nd$sex=ifelse(nd$Male==0,'F','M')
probs <- marginaleffects::predictions(Markov, 
                                      newdata = nd,
                                      type = "probs")

#visualise predicted values
probs$Category=paste(probs$l1,probs$sex)
probs$age=probs$baseline_age+probs$time
min(probs$conf.low)
max(probs$conf.high)
names(predicted)[14:16]=c('P_absent','P_poor','P_normal')
smooth_absent=ggplot(subset(probs,probs$group=='0.Absent'), aes(x = age, y =estimate,colour= Category))+
  geom_smooth(method = "loess", level=0.95)+theme_bw()+ylab('Predicted probability')+ ylim(0, 1)+ xlim(0, 18)+ggtitle('Predicted Probability Absent Speech')
smooth_absent <- ggplot(subset(probs, probs$group == '0.Absent'), 
                        aes(x = age, y = estimate, colour = Category)) +
  geom_smooth(method = "loess", level = 0.95, se = FALSE) + # Remove default CI from geom_smooth
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = Category), alpha = 0.2) + # Add custom CI
  theme_bw() +
  ylab('Predicted probability') + 
  ylim(0, 1) + 
  xlim(0, 15) + 
  ggtitle('Predicted Probability Absent Speech')
smooth_poor <- ggplot(subset(probs, probs$group == '1.Poor'), 
                        aes(x = age, y = estimate, colour = Category)) +
  geom_smooth(method = "loess", level = 0.95, se = FALSE) + # Remove default CI from geom_smooth
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = Category), alpha = 0.2) + # Add custom CI
  theme_bw() +
  ylab('Predicted probability') + 
  ylim(0, 1) + 
  xlim(0, 15) + 
  ggtitle('Predicted Probability Poor Speech')

ggarrange(smooth_absent,smooth_poor, common.legend = TRUE, legend="bottom")

