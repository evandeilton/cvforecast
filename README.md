# Crosss-Validation Forecast (cvforecast)
Forecast time-series by cross-validation method and choose the best models automatically.
It contains functions to perform cross-validation and to help decision on the best models according to some
goodness of fit statistics, say: MAPE, MAE, RMSE, etc. Ideas like linearity and trend are used to help decision process.
Some R code are functions adapted from the code in https://github.com/zachmayer/cv.ts (thanks)

To install code run

```{R}
devtools::install_github('evandeilton/cvforecast')
```

Note: There are lots of undocumented code. This is an very unsafe experiment I'm not a R programming guy!
