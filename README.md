# Cross-Validation Forecast (cvforecast)
Forecast time-series by cross-validation method and choose the best models automatically.
It contains functions to perform cross-validation and to help decision on the best models according to some
goodness of fit statistics, say: MAPE, MAE, RMSE, etc. Ideas like linearity and trend are used to help decision process.
Some R code are functions adapted from the code in https://github.com/zachmayer/cv.ts (thanks)

To install code run
```{R}
devtools::install_github('evandeilton/cvforecast')
```
Note: There are lots of undocumented code. This is a very unsafe experiment I'm not a R programming guy!

# Examples
##  cvforecast

This is the core function of the package. It computes multiple forecasts by the technique of Cross-Validation. The decision about the best models is based on linearity, trend, fit accuracy as for as residual analysis.

#### Usage
```{R}
cvforecast(tsdata, tsControl = cvForecastControl(), fcMethod = NULL, ...)
```

#### Run cvforecast example
Define cross validation parameters

```{R}
require("cvforecast")
myControl <- cvForecastControl(
 minObs = 14,                     # minimum of observations
 stepSize = 10,                   # step size for resampling
 maxHorizon = 30,                 # forecast horizon
 summaryFunc=tsSummary,           # function to sumarize cross-validation accuracy
 cvMethod="MAPE",                 # accuracy statistic for decicion
 tsfrequency='day',               # data points frequencies
 OutlierClean=FALSE,              # clean outlier in the data preparation
 dateformat='%d/%m/%Y %H:%M:%S')  # date format as it is char
```
Paralell execution improves the processing time

```{R}
require("doParallel")             # extra package for paralelization
cl <- makeCluster(4, type='SOCK') # 4 is the number os logical cores. Edit as your own!
registerDoParallel(cl)            # Register cluster
```
Load data and convert to 'ts'
```{R}
data(datasample, package="cvforecast")
dadosd <- ConvertData(datasample[,1:6], dateformat='%d/%m/%Y %H:%M:%S', tsfrequency = "day", OutType="ts")
table(sapply(dadosd, class))      # check class to confirm conversions
dim(dadosd)
```
Looping for several forecasts
```{R}
require("plyr")
FF <- llply(dadosd[,1:5], function(X) {
 fit <- try(cvforecast(X, myControl))
 if(class(fit) != "try-error") {
   #plot(fit)
   return(fit)
 } else NA
}, .progress = "time")

table(sapply(FF, class))
plot(FF[[1]])
sapply(FF, names)
stopCluster(cl)
```
