### retrieving the data ###

## files from Google Books Ngram
pos = read.table(gzfile('pos.gz'))
h = read.table(gzfile('h.gz'))
s = read.table(gzfile('s.gz'))

## data for all pronouns
pron = pos[pos[,1]=='_PRON_',2:3]
names(pron) = c('yr','occurence')
pron = pron[pron[,1] >= 1750,]

## data for 'he'
he = h[h[,1] == 'he', 2:3]
he = he[he[,1] >=1750,]
# generating 'he' to all pronouns ratio
for (i in 1:length(pron[,1])){
  he[i,3] = he[i,2]/pron[i,2]
}
names(he) = c('yr','occurrence','ratio')


## data for 'she'
she = s[s[,1] == 'she', 2:3]
she = she[she[,1] >=1750,]
# generating 'she' to all pronouns ratio
for (i in 1:length(pron[,1])){
  she[i,3] = she[i,2]/pron[i,2]
}
names(she) = c('yr','occurrence','ratio')


### 'he' and 'she' time series ###
heratio = ts(he[,3], start = 1750, end = 2008)
sheratio = ts(she[,3], start = 1750, end = 2008)

## plot of 'he' ratio and 'she' ratio series on same plot
plot(heratio,type='l',ylim = c(0,0.11),main="Figure 1. Usage of 'he', 'she' over time",ylab='Ratio')
lines(sheratio,col = 'red')
legend(1750,0.06,c("'he' to all pronouns","'she' to all pronouns"),col=c('black','red'),lty='solid',cex=0.8)
# both series are not stationary


### examining ccf of heratio and sheratio ###
crosscor = ccf(heratio,sheratio,ylim=c(-.7,.2),lag.max = 30,main="Figure 2. Correlation between 'he' & 'she'")
# negative correlations --> implies inverse relation


### linear model ###
# 'he'
hr_linmod = lm(heratio ~ time(heratio))
plot(heratio)
abline(hr_linmod,col='red')
# 'she'
sr_linmod = lm(sheratio ~ time(sheratio))
plot(sheratio)
abline(sr_linmod,col='red')


## periodograms ##
# we want to determine whether finding a peak frequency
# and whether the harmonic model generated from that frequency would be a good fit for the data

## periodogram for the 'he' series
hr_pdgm = periodogram(heratio,main="Figure 3a. Periodogram for 'he' series",xlim=c(0,0.1))
hr_freq = hr_pdgm$freq
hr_spec = hr_pdgm$spec
hr_peakfreq = hr_freq[which(hr_spec == max(hr_spec[1:20]))]
# peak frequency is 0.007407407

## harmonic model with peak frequency
hr_harmod = lm(heratio ~ cos(2*pi*hr_peakfreq*time(heratio)) + sin(2*pi*hr_peakfreq*time(heratio)))
plot(heratio)
points(x=time(heratio),y=hr_harmod$fitted.values,col='red',type='l')
BIC(hr_harmod)
# -1841.854 --> lower than ARIMA model
summary(hr_harmod)
#Coefficients:
#                                            Estimate Std. Error t value Pr(>|t|)    
#(Intercept)                                0.0816697  0.0004146 197.004   <2e-16 ***
#cos(2 * pi * hr_peakfreq * time(heratio)) -0.0009079  0.0005927  -1.532    0.127    
#sin(2 * pi * hr_peakfreq * time(heratio))  0.0054830  0.0005797   9.458   <2e-16 ***

#Residual standard error: 0.006659 on 256 degrees of freedom
#Multiple R-squared: 0.262,  Adjusted R-squared: 0.2562 
#F-statistic: 45.43 on 2 and 256 DF,  p-value: < 2.2e-16 


## periodogram for the 'she' series
sr_pdgm = periodogram(sheratio,main="Figure 3b. Periodogram for 'she' series")
sr_freq = sr_pdgm$freq
sr_spec = sr_pdgm$spec
sr_peakfreq = sr_freq[which(sr_spec == max(sr_spec[1:20]))]
# peak frequency is 0.003703704

## harmonic model with peak frequency
sr_harmod = lm(sheratio ~ cos(2*pi*sr_peakfreq*time(sheratio)) + sin(2*pi*sr_peakfreq*time(sheratio)))
plot(sheratio)
points(x=time(sheratio),y=sr_harmod$fitted.values,col='red',type='l')
BIC(sr_harmod)
# -2055.447 --> lower than ARIMA model
summary(sr_harmod)
#Coefficients:
#                                            Estimate Std. Error t value Pr(>|t|)    
#(Intercept)                                0.0164374  0.0002745   59.88  < 2e-16 ***
#cos(2 * pi * sr_peakfreq * time(sheratio)) 0.0017366  0.0003955    4.39 1.66e-05 ***
#sin(2 * pi * sr_peakfreq * time(sheratio)) 0.0054381  0.0003807   14.29  < 2e-16 ***

#Residual standard error: 0.004409 on 256 degrees of freedom
#Multiple R-squared: 0.469,  Adjusted R-squared: 0.4649 
#F-statistic: 113.1 on 2 and 256 DF,  p-value: < 2.2e-16 


### ACF and PACF ###

## ACF and PACF for 'he' series
acf(heratio,lag.max = 50,main="Figure 4a. ACF of 'he' series")
# many nonzero components
pacf(heratio,lag.max = 50,main="Figure 4b. Partial ACF of 'he' series")
# one obviously nonzero component, a few less obvious nonzero components

## ACF and PACF for 'she' series
acf(sr_linmod$residuals,lag.max=50,main="Figure 4c. ACF for 'he' series")
# many nonzero components, uneven cyclic pattern, slowly approaching zero
pacf(sr_linmod$residuals,lag.max=50,main="Figure 4d. Partial ACF for 'he' series")
# three obviously nonzero components, a few less obvious nonzero components
# -->ARMA or ARIMA best?


###ARIMA###

# first, for xreg portion of ARIMA, create a time series for the feminism movements
# indicator: 1 if feminist movement is happening, 0 if not
x = rep(0,length(heratio))
# first-wave feminism
for(i in 1848:1920){
  x[which(time(heratio) == i)] = 1
}
# second-wave feminism
for(i in 1963:1982){
  x[which(time(heratio) == i)] = 1
}
# third-wave feminism
for(i in 1990:2008){
  x[which(time(heratio) == i)] = 1
}

## finding an ARIMA model for the 'he' series
hr_bicARIMA0 = matrix(nrow=6,ncol=6)
for(p in 0:5){ 
  for(q in 0:5){
    hr_modARIMA0 = arima(heratio,order=c(p,0,q),xreg=data.frame(time=time(heratio),fem=x))
    hr_bicARIMA0[p+1,q+1] = -2*hr_modARIMA0$loglik + length(hr_modARIMA0$coef) * log(length(heratio)) 
  }
}
hr_bicARIMA1 = matrix(nrow=6,ncol=6)
for(p in 0:5){ 
  for(q in 0:5){
    hr_modARIMA1 = arima(heratio,order=c(p,1,q),xreg=data.frame(time=time(heratio),fem=x))
    hr_bicARIMA1[p+1,q+1] = -2*hr_modARIMA1$loglik + length(hr_modARIMA1$coef) * log(length(heratio)) 
  }
}
hr_bicARIMA2 = matrix(nrow=6,ncol=6)
for(p in 0:5){ 
  for(q in 0:5){
    hr_modARIMA2 = arima(heratio,order=c(p,2,q),xreg=data.frame(time=time(heratio),fem=x))
    hr_bicARIMA2[p+1,q+1] = -2*hr_modARIMA2$loglik + length(hr_modARIMA2$coef) * log(length(heratio)) 
  }
}
hr_bicARIMA3 = matrix(nrow=6,ncol=6)
for(p in 0:5){ 
  for(q in 0:5){
    hr_modARIMA3 = arima(heratio,order=c(p,3,q),xreg=data.frame(time=time(heratio),fem=x))
    hr_bicARIMA3[p+1,q+1] = -2*hr_modARIMA3$loglik + length(hr_modARIMA3$coef) * log(length(heratio)) 
  }
}
min(hr_bicARIMA0)
#-2090.091
min(hr_bicARIMA1)
#-2086.165
min(hr_bicARIMA2)
#-2065.049
min(hr_bicARIMA3)
#-2035.445
which(hr_bicARIMA0 == min(hr_bicARIMA0))
# ARIMA with i = 0 has most negative BIC
# determine which order model has this BIC
# -->ARIMA(1,0,1)

## 'he' ARIMA model
hr_modARIMA = arima(heratio,order=c(1,0,1),xreg=data.frame(time=time(heratio),fem=x))
hr_modARIMA
#Coefficients:
#         ar1      ma1  intercept    time      fem
#      0.9406  -0.5528     0.1901  -1e-04  -0.0012
#s.e.     NaN   0.0519     0.0416   0e+00   0.0015

#sigma^2 estimated as 1.639e-05:  log likelihood = 1058.94,  aic = -2105.88

## 'he' ARIMA residual plots
plot(hr_modARIMA$residuals,main="F",ylab='')
abline(0,0,col='red',lty=2)
hist(hr_modARIMA$residuals,breaks=20,ylim=c(0,80),
     xlab='Residuals',main='Figure 5b. Histogram of ARIMA(1,0,1) residuals')
qqPlot(hr_modARIMA$residuals,ylab='Residuals',main='Figure 5c. QQ plot for ARIMA(1,0,1) Residuals')
acf(hr_modARIMA$residuals,main='Figure 5d. ACF of ARIMA(1,0,1) residuals')
pacf(hr_modARIMA$residuals,main='Figure 5e. Partial ACF of ARIMA(1,0,1) residuals')

## finding an ARIMA model for the 'she' series
sr_bicARIMA0 = matrix(nrow=6,ncol=6)
for(p in 0:5){ 
  for(q in 0:5){
    sr_modARIMA0 = arima(sheratio,order=c(p,0,q),xreg=data.frame(time=time(sheratio),fem=x))
    sr_bicARIMA0[p+1,q+1] = -2*sr_modARIMA0$loglik + length(sr_modARIMA0$coef) * log(length(sheratio)) 
  }
}
sr_bicARIMA1 = matrix(nrow=6,ncol=6)
for(p in 0:5){ 
  for(q in 0:5){
    sr_modARIMA1 = arima(sheratio,order=c(p,1,q),xreg=data.frame(time=time(sheratio),fem=x))
    sr_bicARIMA1[p+1,q+1] = -2*sr_modARIMA1$loglik + length(sr_modARIMA1$coef) * log(length(sheratio)) 
  }
}
sr_bicARIMA2 = matrix(nrow=6,ncol=6)
for(p in 0:5){ 
  for(q in 0:5){
    sr_modARIMA2 = arima(sheratio,order=c(p,2,q),xreg=data.frame(time=time(sheratio),fem=x))
    sr_bicARIMA2[p+1,q+1] = -2*sr_modARIMA2$loglik + length(sr_modARIMA2$coef) * log(length(sheratio)) 
  }
}
sr_bicARIMA3 = matrix(nrow=6,ncol=6)
for(p in 0:5){ 
  for(q in 0:5){
    sr_modARIMA3 = arima(sheratio,order=c(p,3,q),xreg=data.frame(time=time(sheratio),fem=x))
    sr_bicARIMA3[p+1,q+1] = -2*sr_modARIMA3$loglik + length(sr_modARIMA3$coef) * log(length(sheratio)) 
  }
}
min(sr_bicARIMA0)
#-2685.735
min(sr_bicARIMA1)
#-2678.591
min(sr_bicARIMA2)
#-2657.857
min(sr_bicARIMA3)
#-2627.043
which(sr_bicARIMA0 == min(sr_bicARIMA0))
# ARIMA model with i = 0 has most negative BIC
# -->ARIMA(4,0,1)

## 'she' ARIMA model
sr_modARIMA = arima(sheratio,order=c(4,0,1),xreg=data.frame(time=time(sheratio),fem=x))
sr_modARIMA
#Coefficients:
#         ar1      ar2     ar3      ar4      ma1  intercept   time     fem
#      1.3822  -0.1236  0.0220  -0.2883  -1.0000    -0.1132  1e-04  -1e-04
#s.e.  0.0596   0.1071  0.1118   0.0629   0.0105     0.0059  0e+00   3e-04

#sigma^2 estimated as 1.519e-06:  log likelihood = 1365.09,  aic = -2714.19

##'she' ARIMA residual plots
plot(sr_modARIMA$residuals,ylab='',main='Figure 6a. Residuals for ARIMA(4,0,1)')
abline(0,0,col='red',lty=2)
hist(sr_modARIMA$residuals,breaks = 20,xlab='Residuals',ylim=c(0,50),
     main='Figure 6b. Histogram of ARIMA(4,0,1) residuals')
qqPlot(sr_modARIMA$residuals,ylab='Residuals',main='Figure 6c. QQ plot for ARIMA(4,0,1) residuals')
acf(sr_modARIMA$residuals,lag.max=100,main='Figure 6d. ACF of ARIMA(4,0,1) residuals')
pacf(sr_modARIMA$residuals,lag.max=100,main='Figure 6e. Partial ACF of ARIMA(4,0,1) residuals')
