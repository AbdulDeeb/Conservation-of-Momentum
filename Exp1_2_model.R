
#this code, generates figures 2, 3, 6, & 7 and has the data from both Experiment 1 and 2. 

#libraries 
library(nloptr)
library(cowplot)
library(ggplot2)
library(boot)
library(ggpubr)
library(sciplot) # se function
library(gtools)
library(data.table)
library(Metrics)
rm(list=ls())

dir = paste("~/Documents/GitHub/Experiment 1_2/") #change for your own OS
setwd(dir)
#Experiment 2: velocity study data 
load('allData')


ggplot(data=allData,aes(x=as.numeric(Vb*100),y=((probeSpeed)*100), color = as.factor(Ua),
                        group  = as.factor(Ua))) +
  facet_grid(.~VolumeRatio)+   scale_color_manual(values=c('black','darkgray')) +
  stat_summary(fun.data=mean_se,geom='point') +
  stat_summary(fun.data=mean_se,geom='line') +
  stat_summary(fun.data=mean_se,geom='errorbar',width=.01) +
  geom_abline(intercept=0,slope=1,lty='longdash',col='gray') +
  xlab("Displayed Vb" ) + ylab("Probe Speed" )+ coord_equal()+
  scale_x_continuous(breaks = c(36, 61))+
  theme_minimal()



ModelData <- allData
ModelData$vbp <- ModelData$VolumeRatio * (1 + 1) / (1 + ModelData$VolumeRatio ) * ModelData$Ua  # predicted vb
#noise test 
ModelData_sub <- subset(ModelData, VolumeRatio == 1.0000000)

varset = with(ModelData_sub,aggregate(probeSpeed,
                                      by=list(subjName = subjName,
                                              Vb=Vb,
                                              Ua = Ua
                                      ), var))

colnames(varset)[4] <- 'var'
t.test(varset$var[ varset$Vb == 0.6155], varset$var[ varset$Vb == 0.3611], alternative = "greater")

ggplot(data=varset,aes(x=as.numeric(Vb*100),y=(var))) +
  scale_color_manual(values=c('black','darkgray')) +
  stat_summary(fun.data=mean_se,geom='point') +
  stat_summary(fun.data=mean_se,geom='line') +
  stat_summary(fun.data=mean_se,geom='errorbar',width=.01) +
  theme_minimal()


#GET MODEL OF VB
#ssm function set up 
Bm2 = function(guess,vbp,vb,ua,r, data){
  # Parameters
  kp = guess[1] # we are fitting the reliability of the prediction
  vbp = vbp
  vb = vb
  ua = ua 
  e = 1
  # Model Equation
  #vbp = r*(1+e)/(1+r)*ua
  
  vb_hat = sqrt((kp^2* vbp^2) + .75*vb^2)
  
  #vb_hat = kp*vbp + (1-kp)*vb
  simulation_error = sqrt(mean((vb_hat-data)^2,na.rm=T)) #RMSE
  
  return (simulation_error)
}
# Setup for the optimizer

vec = NULL
guess = c(runif(1, min = 0.0, max = 10)) # initial guess for the fit parameter its a weight between 0 and 1 now. 
LB =  c(0.0)
UB =  c(10)
opt_print=0
options = list('algorithm'='NLOPT_LN_COBYLA',#'NLOPT_LN_BOBYQA',#'NLOPT_LN_NELDERMEAD',#'NLOPT_LN_SBPLX','NLOPT_LN_COBYLA'
               'print_level'=opt_print,
               'maxtime'=240) # better not take this long

# Do the Optimization
outdata=NULL
fitParams_bySubj = NULL
for (sb in unique(ModelData$subjName)){
  tempfit = subset(ModelData, subjName==sb)
  soln = nloptr(x0 = guess,
                eval_f = Bm2,
                lb = LB,
                ub = UB,
                opts=options,
                ua = as.numeric(tempfit$Ua2),
                r = as.numeric(tempfit$VolumeRatio),
                vb = as.numeric(tempfit$Vb),
                vbp = tempfit$vbp,
                data = tempfit$probeSpeed)
  
  tempfit$kp=soln$solution
  
  fitParams_bySubj=c(fitParams_bySubj,soln$solution)
  
  outdata = rbind(outdata, tempfit) 
}
normalizeFactor <- 0.03889429
outdata$model_fit = sqrt((outdata$kp^2* outdata$vbp^2) +  .75*outdata$Vb^2)+ normalizeFactor
outdata$probeSpeed <- outdata$probeSpeed + normalizeFactor
#outdata$model_fit = outdata$kp*outdata$vbp +  (1-outdata$kp)*outdata$Vb


ggplot(data=outdata, aes(x=as.numeric(Vb), y=probeSpeed, color = as.factor(Ua), group = as.factor(Ua))) +   
  
  # Facet by VolumeRatio and apply custom labels
  facet_grid(.~VolumeRatio, 
             labeller = as_labeller(c(
               `0.3333333` = "1/3",  # Adjust based on actual values
               `1` = "1/1",
               `3` = "3/1"))) + 
  
  
  # Model fit points and dashed line with transparency (removed 'fill' aesthetic)
  stat_summary(fun.data = 'mean_se', aes(y = model_fit), geom = 'point', alpha = .3) +
  stat_summary(fun.data = 'mean_se', aes(y = model_fit), geom = 'line', linetype = "dashed") +
  
  # Empirical points, error bars, and solid line with higher opacity (removed 'fill' aesthetic)
  stat_summary(fun.data = 'mean_se', geom = 'point', alpha = .9) +   
  stat_summary(fun.data = 'mean_se', geom = 'errorbar', width = .05) +   
  stat_summary(fun.data = 'mean_se', geom = 'line', alpha = .9) +
  
  # Customize x and y axis labels with subscripts
  xlab(expression(Displayed~V[b])) +   
  ylab(expression(Probe~V[b])) +   
  
  # Apply theme_classic and customize font sizes
  theme_classic() +
  theme(
    aspect.ratio = 1,  # Adjust the aspect ratio (1 for square, adjust as needed)
    axis.title = element_text(size = 14),       # Font size for axis labels
    axis.text = element_text(size = 12),        # Font size for tick marks
    strip.text = element_text(size = 14),       # Font size for facet labels
    legend.text = element_text(size = 12),      # Font size for legend text
    legend.title = element_text(size = 14)      # Font size for legend title
  ) +
  
  # Use grayscale color palette
  scale_color_grey(start = 0.7, end = 0.1) +  # Lighter grayscale range
  
  # Customize the legend title to show u with a subscript a and unit cm/s
  labs(color = expression(u[a]~"(cm/s)")) +   # Legend title with subscript and unit
  
  # Add a title to the plot
  ggtitle("IC Model on Vb")

#fitting the Mass perception 
kpValues <- unique(outdata$kp)

load('MassJudgmentData_full')#Experiment 1 data. 
load('MassJudgmentData')  #Also Experiment 1, but less columns. 
P1 <- ggplot(newdata, aes(x=as.factor(MassRatio), 
                          y=as.numeric(ProportionMa), 
                          color=as.factor(VolumeRatio), 
                          group=as.factor(VolumeRatio))) + 
  geom_hline(yintercept = .5, lty='longdash', col='grey') + 
  stat_summary(fun.data='mean_se', geom='point', size = 2) + 
  stat_summary(fun.data='mean_se', geom='line', size = 1) + 
  stat_summary(fun.data='mean_se', geom='errorbar', width = .3, size = 1) + 
  theme_classic() + 
  coord_cartesian() + 
  ggtitle("Data") + 
  ylim(0,1) +
  
  # Customize axis labels with expression() for italic text and subscripts
  labs(
    y = expression(Proportion~choosing~italic(a)),  # Italic a for y-axis
    x = expression(M[a] / M[b])                     # M_a / M_b for x-axis
  ) +
  
  # Adjust font size for axis labels and tick marks
  theme(
    axis.title = element_text(size = 14),            # Axis labels font size
    axis.text = element_text(size = 10)              # Tick marks font size
  )

# Display the plot
P1

#Model mass ratio, using experiment 2 data. 
ModelData2 <- allData
#adjust ua
ModelData2$Ua2 <- ModelData2$Ua
ModelData2$Ua2 <- ModelData2$Ua2*-1
#adjust vb and va
ModelData2$Va <- ModelData2$Va *-1
ModelData2$Vb <- ModelData2$Vb *-1
#get volume ratio as number 
ModelData2$VolumeRatio <- ModelData2$VolumeRatioIndex
ModelData2$VolumeRatio[ModelData2$VolumeRatioIndex == 0] = round(1/3, digits = 2)
ModelData2$VolumeRatio[ModelData2$VolumeRatioIndex == 1] = 1
ModelData2$VolumeRatio[ModelData2$VolumeRatioIndex == 2] = 3

#add resposnes 

ModelData2$ProportionMa <-  newdata$ProportionMa

#get mass ratio as number 
ModelData2$MassRatio <- ModelData2$MassRatioIndex
ModelData2$MassRatio[ModelData2$MassRatioIndex == 1] = round(1/3, digits = 2)
ModelData2$MassRatio[ModelData2$MassRatioIndex == 2] = round(2/3, digits = 2)
ModelData2$MassRatio[ModelData2$MassRatioIndex == 3] = 3/2
ModelData2$MassRatio[ModelData2$MassRatioIndex == 4] = 3/1



#predict vb 
ModelData2$vbp <- ModelData2$VolumeRatio * (1 + 1) / (1 + ModelData2$VolumeRatio ) * ModelData2$Ua2  # predicted vb
#predict va
ModelData2$vap <- (ModelData2$VolumeRatio -1)/(1+ModelData2$VolumeRatio)*ModelData2$Ua2;   

#assign kp values
ModelData2$kp <- mean(kpValues)
#combine vbp + vb 
ModelData2$vb_IC <- sqrt(((ModelData2$kp* ModelData2$vbp)^2) +  ModelData2$Vb^2)
#combine vap + va 
ModelData2$va_IC <- 0
ModelData2$va_IC[sign(ModelData2$Va) == sign(ModelData2$vap) ] = sign(ModelData2$Va)* 
  sqrt(((ModelData2$kp* ModelData2$vap)^2) +  ModelData2$Va^2)
ModelData2$va_IC[sign(ModelData2$Va) != sign(ModelData2$vap) ] <- ModelData2$Va

#get a mass ratio 
ModelData2$PredMassRatio <- ModelData2$vb_IC / (ModelData2$Ua2 - ModelData2$va_IC)
ModelData2$delta = (ModelData2$PredMassRatio-1)/(ModelData2$PredMassRatio+1)

#fitting the decision noise
Bm2 = function(guess, PredMassRatio, data){
  # Parameters
  dN = guess[1] # we are fitting the reliability of the prediction
  PredMassRatio = PredMassRatio
  
  # Model Equation
  pa = pnorm(PredMassRatio, 1, dN) 
  simulation_error = sqrt(mean((pa-data)^2,na.rm=T)) #RMSE
  return (simulation_error)
}

# Setup for the optimizer
vec = NULL
guess = c(runif(1, min = 0.0, max = 0.2)) # initial guess for the fit parameter its a weight between 0 and 1 now. 
LB =  c(0.0)
UB =  c(0.2)
opt_print=0
options = list('algorithm'='NLOPT_LN_COBYLA',#'NLOPT_LN_BOBYQA',#'NLOPT_LN_NELDERMEAD',#'NLOPT_LN_SBPLX','NLOPT_LN_COBYLA'
               'print_level'=opt_print,
               'maxtime'=240) # better not take this long

# Do the Optimization
outdata=NULL
fitParams_bySubj = NULL
for (sb in unique(ModelData2$subjName)){
  tempfit = subset(ModelData2, subjName==sb)
  soln = nloptr(x0 = guess,
                eval_f = Bm2,
                lb = LB,
                ub = UB,
                opts=options,
                PredMassRatio = tempfit$PredMassRatio,
                data = tempfit$ProportionMa)
  
  tempfit$dN =soln$solution
  
  fitParams_bySubj=c(fitParams_bySubj,soln$solution)
  
  outdata = rbind(outdata, tempfit) 
}

ModelData2$dN <- outdata$dN
ModelData2$Pa <- pnorm(ModelData2$PredMassRatio, 1, ModelData2$dN ) 
#make model plot
P2 <- ggplot(data=ModelData2,aes(x= as.factor(MassRatio),y=Pa, group = as.factor(VolumeRatio), color = as.factor(VolumeRatio))) +
  geom_hline(yintercept = .5, lty='longdash',col='grey') +    
  stat_summary(fun.data='mean_se',geom='line', size = 1) + ylab("Proportion of Choosing Motor Object" ) + xlab( "Ma/Mb" )+
  theme_minimal()  + ggtitle("SPoD Model") + ylim(0,1)

P2 <- ggplot(ModelData2, aes(x=as.factor(MassRatio), 
                             y=as.numeric(Pa), 
                             color=as.factor(VolumeRatio), 
                             group=as.factor(VolumeRatio))) + 
  geom_hline(yintercept = .5, lty='longdash', col='grey') + 
  stat_summary(fun.data='mean_se', geom='line', size = 1) + 
  theme_classic() + 
  coord_cartesian() + 
  ggtitle("Model") + 
  ylim(0,1) +
  
  # Customize axis labels with expression() for italic text and subscripts
  labs(
    y = expression(Proportion~choosing~italic(a)),  # Italic a for y-axis
    x = expression(M[a] / M[b])                     # M_a / M_b for x-axis
  ) +
  
  # Adjust font size for axis labels and tick marks
  theme(
    axis.title = element_text(size = 14),            # Axis labels font size
    axis.text = element_text(size = 10)              # Tick marks font size
  )

# Display the plot
P2

ggarrange(P1, P2, common.legend = TRUE)





