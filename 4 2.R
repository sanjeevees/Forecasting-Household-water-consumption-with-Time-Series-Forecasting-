library(car)
setwd("C:\\users\\Sanjeevees\\Desktop")

##
## Standard Regression Diagnositcs
##
setwd("C:\\users\\Sanjeevees\\Desktop")
concord<- foreign::read.spss("ConcordH2OTimeSeries.sav", to.data.frame = TRUE)
model<- lm(H2Ouse~temp+rain+educ+time,concord )
summary(model)

model.res<- ts(residuals(model))
acf(model.res)


## Generate harmonic variables
concord$sin6 <- sin(concord$time/6*2*pi)
concord$cos6 <- cos(concord$time/6*2*pi)
concord$cos12 <- cos(concord$time/12*2*pi)
concord$sin12 <- sin(concord$time/12*2*pi)

## Augmented Model
fourier.lm <- lm(H2Ouse~temp+rain+educ+time+sin6+cos6+cos12+sin12,concord)
summary(fourier.lm)  

vif(fourier.lm)

ts.res <- ts(residuals(fourier.lm))
acf(ts.res)
forecast::auto.arima(ts.res)

model2<-lm(H2Ouse~educ+time+sin6+cos6+cos12+sin12,concord)
model3<-lm(H2Ouse~educ+time,concord)
anova( model3,model)


res2<-residuals(model2)
acf(res2)
set.seed(NULL)
durbinWatsonTest(model2, max.lag=12, simulate=TRUE, alternative="positive")

ts.res2<-ts(residuals(model2),start=c(1970, 1))
forecast::auto.arima(ts.res2,)

library(nlme)
fourier.gls2 <- gls(H2Ouse~educ+time+sin6+cos6+cos12+sin12, data=concord, 
                    correlation=corARMA(p=1,q=1), method="ML")
summary(fourier.gls2)

ts.gls.res <- ts(residuals(fourier.gls2))

ts.gls.res <- ts(residuals(fourier.gls2, type="normalized"),      # type="normalized": AC adjusted residuals 
                 frequency = 12, start = c(1970, 1)) 

acf(ts.gls.res)
pacf(ts.gls.res)

sqrt((fourier.lm$coef["sin6"]^2+fourier.lm$coef["cos6"]^2))

sqrt((fourier.lm$coef["sin12"]^2+fourier.lm$coef["cos12"]^2))

deltaMethod()


ampYse <- deltaMethod(fourier.lm,"sqrt(cos6^2+sin6^2)")
ampHse <- deltaMethod(fourier.lm,"sqrt(sin12^2+cos12^2)")

fourierResult <- rbind(ampYse,ampHse)
row.names(fourierResult) <- c("Cycle 12 Months:","Cycle 6 Months:")
round(fourierResult,3)





## Durbin-Watson Test
durbinWatsonTest(fourier.lm, max.lag=6, simulate=TRUE, alternative="positive")

##
## Time Series Analysis. See Kleiber & Zeileis pp151-164
##
ts.avg <- ts(sandusky$avg7447, frequency = 12, start = c(1990, 1)) # convert to a ts-object
ts.avg <- ts(sandusky$avg7447, deltat = 1/12, start = c(1990, 1))  # equivalent alternative specification
plot(ts.avg)
spectrum(ts.avg)                     # spectral periodogram (cor^2 between series and sine/cosine-terms at given frequency)
abline(v=1, h=1, lty=3)
plot(stl(ts.avg,"periodic"))         # decomposition into periodic, trend, remainder with "stl"
                                     # see pp 155-156 KZ. Bar are of equal hight in user coordinates
## Moving average filter smoothing
ts.ma <- filter(ts.avg, filter=rep(1/3,3), sides=2)             # moving average smoothing with weights 1/3
(ts.both <- ts.union(ts.avg,ts.ma))                             # bind both series

## Time Series plot of original and smoothed series
plot(ts.both, plot.type = "single", ylab="Temperatur Residuals & Moving Average", 
     xlab="Month & Year", main="Smoothed (blue) and Original Series (red)",
     pch=20, type="b", lwd=c(1,1), col=c("red","blue"))
abline(h=mean(sandusky$avg7447),lty=2,lwd=2)                                      

## Dickey-Fuller stationary  test at lag k=1
ts.res <- ts(residuals(fourier.lm), start=c(1990, 1))
tseries::adf.test(ts.res, k=1)

## Automated ARMA order identification                                  
forecast::auto.arima(ts.res) 

pt()
##
## Full ARMA model estimation using GLS 
##
library(nlme)
## Plain model
fourier.gls0 <- gls(avg7447~time.idx+I(time.idx^2)+r.cos+r.sin, data=sandusky,
                    method="ML")
summary(fourier.gls0)

## AR(1) model
fourier.gls1 <- gls(avg7447~time.idx+I(time.idx^2)+r.cos+r.sin, data=sandusky, 
                    correlation=corARMA(p=1), method="ML")
summary(fourier.gls1)
anova(fourier.gls1, fourier.gls0)

## AR(1) & MA(1) model
fourier.gls2 <- gls(avg7447~time.idx+I(time.idx^2)+r.cos+r.sin, data=sandusky, 
                    correlation=corARMA(p=1,q=1), method="ML")
summary(fourier.gls2)
anova(fourier.gls2, fourier.gls1)

## Check residuals of corrected model
ts.gls.res <- ts(residuals(fourier.gls2, type="normalized"),      # type="normalized": AC adjusted residuals 
                 frequency = 12, start = c(1970, 1)) 
acf(ts.gls.res, xlab="lag in years", main="Adjusted Residuals AC-Function")
pacf(ts.gls.res, xlab="lag in years", main="Adjusted Residuals PAC-Function")

## Alternative estimation approach using function arima( )
mod.ARMA <- arima(ts.avg, order=c(1,0,1), xreg=model.matrix(fourier.lm), include.mean=F)
mod.ARMA$coef                                 # coefficients
mod.ARMA$coef/sqrt(diag(mod.ARMA$var.coef))   # t-values
