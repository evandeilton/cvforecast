#' Default Cross-validation control
#' @param stepSize size step for the cross-validation samples
#' @param maxHorizon forecasting horizon
#' @param minObs minumum number of observation. Default is two times cycle of data
#' @param fixedWindow keep fixed the sampling window, default is TRUE
#' @param summaryFunc extra function to compute statistics of the model
#' @param preProcess if TRUE does Box-Cox data transformation in to the data
#' @param ppMethod if 'preProcess' is TRUE make 'guerrero' or 'loglik' tranformation. See \code{\link{BoxCox.lambda}}
#' @param cvMethod accuracy method for best model choice. See \code{\link{accuracy}}
#' @param tsfrequency time series data frequency
#' @param OutlierClean if TRUE, remove outliers from the data. See \code{\link{tsclean}}
#' @param residlevel confidence level for residual tests
#' @param dateformat date format for charater dates
#' @return list of parameters
#' @author LOPES, J. E.
#' @examples
#' # Control
#' myControl <- cvForecastControl(
#' minObs = 14,
#' stepSize = 5,
#' maxHorizon = 30,
#' summaryFunc=tsSummary,
#' cvMethod="MAPE",
#' tsfrequency='day',
#' OutlierClean=FALSE)
#' myControl
#' @export
cvForecastControl <- function (stepSize = 1, maxHorizon = 1, minObs = 7, fixedWindow = TRUE, summaryFunc = tsSummary, preProcess = FALSE, ppMethod = "guerrero", cvMethod="MAPE", tsfrequency="month", OutlierClean = TRUE, residlevel = 0.10, dateformat='%d/%m/%Y %H:%M:%S') {
    list(stepSize = stepSize, maxHorizon = maxHorizon, minObs = minObs,
        fixedWindow = fixedWindow, summaryFunc = summaryFunc,
        preProcess = preProcess, ppMethod = ppMethod, cvMethod=cvMethod,
		tsfrequency=tsfrequency, OutlierClean=OutlierClean, residlevel=residlevel, dateformat=dateformat)
}

#' Function to cross-validate a time series
#'
#' Main function to perform cross-validation. This function was firstly created by Zach Mayer (https://github.com/zachmayer/cv.ts) thanks, and adapted by LOPES. J. E. It
#' @param x time series object of class 'ts'
#' @param FUN forecast wrapper function. These are some ones fc_sts, fc_hw, fc_tbats, fc_auto.arima, fc_ses, fc_mean, fc_holt, fc_bats, fc_ets, fc_lm, fc_theta, fc_rw, sfc_naive, fc_naive, fc_nnetar, fc_hws, fc_hwns and fc_hwes. This function works also in parallel very fast in multiple core computers.
#' @param tsControl Generic control for cross-validation process. See \code{\link{cvForecastControl}}.
#' @param progress if TRUE, show the progress bar.
#' @param packages extra R packages required by R. Default is NULL
#' @param ... extra args, if needed
#' @return list with information about then cross-validation like forecasts and accuracy
#' @examples
#' # Control
#' tsControl <- cvForecastControl(
#' minObs = 14,
#' stepSize = 10,
#' maxHorizon = 30,
#' summaryFunc=tsSummary,
#' cvMethod="MAPE",
#' OutlierClean=FALSE)
#' #cl <- makeCluster(4, type='SOCK')
#' #registerDoParallel(cl)
#' x <- AirPassengers
#' fit <- cvts2(x, fc_auto.arima)
#' #stopCluster(cl)
#' @export
cvts2 <- function(x, FUN, tsControl=cvForecastControl(), progress=TRUE, packages=NULL, ...) {

	stopifnot(is.ts(x))

	#Load parameters from the tsControl list
	stepSize <- tsControl$stepSize
	maxHorizon <- tsControl$maxHorizon
	minObs <- tsControl$minObs
	fixedWindow <- tsControl$fixedWindow
	summaryFunc <- tsControl$summaryFunc
	preProcess <- tsControl$preProcess
	ppMethod <- tsControl$ppMethod

	#Define additional parameters
	freq <- frequency(x)
	n <- length(x)
	st <- tsp(x)[1]+(minObs-2)/freq

	#Create a matrix of actual values.
	#X is the point in time, Y is the forecast horizon
	#http://stackoverflow.com/questions/8140577/creating-a-matrix-of-future-values-for-a-time-series
	formatActuals <- function(x,maxHorizon) {
		actuals <- outer(seq_along(x), seq_len(maxHorizon), FUN="+")
		actuals <- apply(actuals,2,function(a) x[a])
		actuals
	}

	actuals <- formatActuals(x,maxHorizon)
	actuals <- actuals[minObs:(length(x)-1),,drop=FALSE]

	#Create a list of training windows
	#Each entry of this list will be the same length, if fixed=TRUE
	steps <- seq(1,(n-minObs),by=stepSize)

	#Set progressbar
	combine <- rbind
	if (progress) {
	  f <- function(){
	    pb <- txtProgressBar(1,length(steps)-1,style=3)
	    count <- 0
	    function(...) {
	      count <<- count + length(list(...)) - 1
	      setTxtProgressBar(pb,count)
	      Sys.sleep(0.01)
	      flush.console()
	      rbind(...)
	    }
	  }
	  combine <- f()
	}

	#At each point in time, calculate 'maxHorizon' forecasts ahead
	forecasts <- foreach(i=steps, .combine=combine,
	                     .multicombine=FALSE,
	                     .packages=c('forecast', packages),
	                     .export=c('testObject', 'tsSummary', 'cvForecastControl')
	                     ) %dopar% {

		if (fixedWindow) {
			xshort <- window(x, start = st+(i-minObs+1)/freq, end = st+i/freq)
		} else {
			xshort <- window(x, end = st + i/freq)
		}
		if (preProcess) {
			if (testObject(lambda)) {
				stop("Don't specify a lambda parameter when preProcess==TRUE")
			}
			stepLambda <- BoxCox.lambda(xshort, method=ppMethod)
			xshort <- BoxCox(xshort, stepLambda)
		}
		out <- FUN(xshort, h=maxHorizon, ...)

		if (preProcess) {
			out <- InvBoxCox(out, stepLambda)
		}
		return(out)
	}

	#Extract the actuals we actually want to use
	actuals <- actuals[steps,,drop=FALSE]

	#Accuracy at each horizon
	out <- data.frame(
		plyr::ldply(1:maxHorizon,
			function(horizon) {
				P <- forecasts[,horizon,drop=FALSE]
				A <- na.omit(actuals[,horizon,drop=FALSE])
				P <- P[1:length(A)]
				P <- na.omit(P)
				A <- A[1:length(P)]
				summaryFunc(P,A)
			}
		)
	)

	#Add average accuracy, across all horizons

	overall <- na.omit(colMeans(out, na.rm = TRUE))
	Sd  <- sapply(out, sd, na.rm=TRUE)
	Md  <- sapply(out, median, na.rm=TRUE)
	out <- rbind(out, overall, Md, Sd)
	results <- data.frame(horizon=c(1:maxHorizon,'Mean','Median','Sd'), out)

	#Add a column for which horizon and output
	return(list(forecasts=forecasts, results=results))
}

#' Box and Cox tests and Ljung
#'
#' Performs Ljung and Box test and also Durbin-Watson both for autocorrelation on the residuals from the model. It returns also accuracy statistics. See \code{\link{accuracy}}, \code{\link{Box.test}}, \code{\link{dwtest}}.
#' @param model A forecast model. See \code{\link{cvts2}}.
#' @param ... extra args, if needed.
#' @author LOPES, J. E.
#' @examples
#' x <- AirPassengers
#' fit <- fc_auto.arima(x, h = 10, onlyfc=FALSE)
#' LjungBtest_Acuracia(fit)
#' plot(fit)
#' @export
LjungBtest_Acuracia <- function(model,...) {
	bt <- Box.test(residuals(model), lag=10, type="Ljung", fitdf=length(coef(model)))

	# Auto-correlation
	LJungBox <- data.frame(LJB_X_quad = round(bt[[1]], 5), LJB_gl = bt[[2]], LJB_p.valor = round(bt[[3]], 5))
	rownames(LJungBox) <- c()
	acc <- data.frame(round(accuracy(model), 5))

	## Serial correlation by Durbin-Watson
	D.Watson <- round(dwtest(model$x~residuals(model))$statistic, 5)

	res <- cbind(acc, LJungBox, D.Watson)
	res
}

#' Correlations by pairs between two class of variables
#'
#' Performs correlations by pairs between two class of variables.
#' @param tabela_y data.frame, vector or matrix for y variables
#' @param tabela_x data.frame, vector or matrix for x variables
#' @param method correlation method required. See \code{\link{cor}} for more.
#' @param digits output length of numeric data
#' @param verbose is TRUE shows information about the process
#' @param corsimples if TRUE, returns simple correlation only
#' @param ... extra args, if needed.
#' @author LOPES, J. E.
#' @return Results as a data.frame with correlation, R^2, adjusted R^2 for a lm model and PRESS statistics.
#' @examples
#' correlation(mtcars[,1:3])
#' correlation(mtcars[,1:3], mtcars[,4:6])
#' @export
correlation <- function(tabela_y, tabela_x=NULL, method = "pearson", digits = 4, verbose = FALSE, corsimples = FALSE, ...) {
	if(is.null(tabela_x)) {
		tab_cor <- round(cor(tabela_y, use = "complete.obs", method=method, ...), digits)
		dados <- tabela_y
	}
	else if(is.null(tabela_y)){
		tab_cor <- round(cor(tabela_x, use = "complete.obs", method=method, ...), digits)
		dados <- tabela_x
	} else {
		tab_cor <- round(cor(tabela_y, tabela_x, use = "complete.obs", method=method, ...), digits)
		dados <- cbind(tabela_y, tabela_x)
	}
	if (corsimples) {
	tab_cor <- round(cor(dados, use = "complete.obs", method=method, ...), digits)
	return(tab_cor)
	} else {
	tab_cor[lower.tri(tab_cor, diag=TRUE)] <- NA 
	tab_cor <- as.data.frame(as.table(tab_cor))
	tab_cor <- na.omit(tab_cor)
	names(tab_cor) <- c("Y","X", "Cor")
	tab_cor <- tab_cor[order(tab_cor$Y, -abs(tab_cor$Cor)),] 
	rownames(tab_cor) <- 1:nrow(tab_cor)
	r2 <- nm <- c()
	for (i in 1:nrow(tab_cor)) {
		formu <- as.formula(paste(as.character(unlist(tab_cor[i,1:2])), collapse="~"))
		mod <- lm(formu, data=dados)
		r2[[i]] <- round(lm_press_stat(mod), digits)
		nm[i]   <- deparse(formu)
		if(verbose) print(summary(mod))
	}
	dcal <- do.call("rbind", r2)
	out  <- data.frame(tab_cor, R2=dcal$r.squared, R2Ajustado=dcal$adj.r.squared, R2PRESS = dcal$press.r.squared)
	rownames(out) <- as.character(nm)
	return(out)
	}
}

#' Find best forecast methods based on trend and linearity
#'
#' Automatically estimates the most adequate forecasts method for a 'ts' object based on linearity and trend tests and return a list containing all functions names. See \code{\link{nparTrend}} and \code{\link{linearityTest}} for more.
#' @param x ts object, univariate time series
#' @author LOPES, J. E.
#' @examples
#' forecastMethod(lynx)
#' forecastMethod(AirPassengers)
#' @export
forecastMethod <- function(x) {
  if(is.null(x) | length(as.numeric(x)) < 8) {
    stop(cat("Need more data!", length(as.numeric(x)), "\n"))
  }
  if (class(x)[1]!='ts' & class(x)[1] == "numeric"){
    x <- ts(x, frequency=1)
  }
  # Short series
  if(length(as.numeric(x)) < 2*max(cycle(x)) | max(cycle(x))==1) {
    short <- TRUE
  } else {
    short <- FALSE
  }
  if (!short) {
    trend <- check_trend <- nlTests <- c()

    # Non-parametric trend tests
    trend <- nparTrend(x)
    check_trend <- unname(trend["trend_sign"] != 0)

    # Linearity tests
    nlTests <- list("terasvirta","white", "keenan", "tsay","tarTest")
    linear  <- Try_error(na.omit(sapply(nlTests,  function(n) linearityTest(x, n)$p.value)))

    if (class(linear)[1] == "numeric") {
      linear <- ifelse(sum(linear > 0.05) < 5, TRUE, FALSE)
    } else {
      linear <- TRUE
    }
  }

  ## Decision
  if (short) {
    ## Short time series or undefined frequencies 
    metodo <- list("fc_auto.arima","fc_ets","fc_lm")
  } else if(linear & !check_trend) {
    ## Linearity without trend
    metodo <- list("fc_naive","fc_rw","fc_sts","fc_theta","fc_hwns","fc_hwes")
  } else if(linear & check_trend) {
    ## Linearity with trend
    metodo <- list("fc_auto.arima","fc_ets","fc_hws","sfc_naive")
  } else {
    ## Non-linearity with or without trend
    metodo <- list("sfc_naive","fc_lm","fc_hws")
  }
  return(metodo)
}

#' Time series forecast residual analysis
#' 
#' Compute six residual test on a forecast object an return its p.values. The null hypothesis that the residual don't have the tested characteristic. If p.value is > 5\% we reject null. 
#' @param forecast object of class forecast
#' @seealso See \code{\link{Box.test}}, \code{\link{t.test}}, \code{\link{LB.test}}, \code{\link{jarque.bera.test}}, \code{\link{bptest}}   and \code{\link{dwtest}} for more information about the tests.
#' @examples
#' fcnm <- forecastMethod(lynx)
#' fun <- get(fcnm[[1]])
#' fit <- fun(AirPassengers, h=10, onlyfc = FALSE)
#' Mresid(fit)
#' @export
Mresid <- function(fit) {
	out <- c()
	x <- na.omit(as.numeric(fit$x))
	r <- na.omit(as.numeric(fit$residuals))

	l <- abs(length(x)-length(r))
	if (l > 0) x <- x[-seq(l)]

	# NULL: Residuals are independent
	independencia <-  Try_error(TSA::LB.test(fit, lag=10, type = "Box-Pierce", no.error=TRUE)$p.value)
	
	# NULL: Residuals are have zero mean
	media_zero <-  Try_error(t.test(r, alternative='two.sided', mu=0.0, conf.level=.95)$p.value)
	
	# NULL: Residuals are white noise
	ruido_branco <-  Try_error(TSA::LB.test(fit, lag=10, type = "Ljung-Box", no.error=TRUE)$p.value)
	
	# NULL: Residuals have Gaussian structure
	normalidade <-  Try_error(tseries::jarque.bera.test(r)$p.value)
	
	# NULL: Residuals are heterokedastic
	homocedasticidade <-  Try_error(lmtest::bptest(r ~ x)$p.value)
	
	# NULL: Residuals have auto-correlation structure
	autocorrelacao <-  Try_error(lmtest::dwtest(r ~ x)$p.value)

	p0 <- if (class(independencia) == "numeric") as.numeric(independencia) else NA
	p1 <- if (class(media_zero) == "numeric") 	 as.numeric(media_zero) else NA
	p2 <- if (class(ruido_branco) == "numeric")  as.numeric(ruido_branco) else NA
	p3 <- if (class(normalidade) == "numeric")   as.numeric(normalidade) else NA
	p4 <- if (class(homocedasticidade) == "numeric") as.numeric(homocedasticidade) else NA
	p5 <- if (class(autocorrelacao) == "numeric") as.numeric(autocorrelacao) else NA

	df.pvalor <- round(data.frame(p0, p1, p2, p3, p4, p5),5)
	colnames(df.pvalor) <- c("Box.test","t.test","LB.test","JB.test","BP.test","DW.test")
	rownames(df.pvalor) <- c("p-value")
	return(df.pvalor)
}

#' Test if an object exists
#' @param object any object present into the environment
#' @return TRE or FALSE
#' @export
testObject <- function(object){
  exists(as.character(substitute(object)))
}

#' Default summary function
#'
#' Return accuracy statistics for the forecast model. See \code{\link{accuracy}}
#' @param P forecast object
#' @param A data vector of sample for testing. Must have same length as P
#' @return accuracy data.frame
#' @examples
#' fcnm <- forecastMethod(lynx)
#' fun <- get(fcnm[[1]])
#' fit <- fun(AirPassengers, h=10, onlyfc = FALSE)
#' class(out <- tsSummary(fit))
#' @export
tsSummary <- function(P, A) {
	data.frame((as.data.frame(accuracy(P,A))))
}

#' Mean forecast wrapper
#' @param x 'ts' data
#' @param h forecast horizon
#' @param level confidence level. Default is 0.95
#' @param onlyfc if TRUE return only forecasts, otherwise returns a full forecast classed object
#' @param ... extra args, if needed.
#'
#' @return forecasts from ts data or an object of class 'forecast'
#' @examples
#' fit <- fc_mean(AirPassengers, h=10, level = 95, onlyfc = FALSE)
#' plot(fit)
#' Mresid(fit)
#' tsSummary(fit)
#' @export
fc_mean <- function(x, h, level = 95, onlyfc=TRUE, ...) {
  fc <- forecast::meanf(y = x, h = h, level = level, ...)
  if (onlyfc) fc$mean
  else fc
}

#' Naive forecast wrapper
#' @param x 'ts' data
#' @param h forecast horizon
#' @param level confidence level. Default is 0.95
#' @param onlyfc if TRUE return only forecasts, otherwise returns a full forecast classed object
#' @param ... extra args, if needed.
#' @return forecasts from ts data or an object of class 'forecast'
#' @examples
#' fit <- fc_naive(AirPassengers, h=10, onlyfc = FALSE)
#' plot(fit)
#' Mresid(fit)
#' tsSummary(fit)
#' @export
fc_naive <- function(x, h, level = 95, onlyfc=TRUE, ...) {
  fc <- forecast::naive(y = x, h = h, level=level, ...)
  if (onlyfc) fc$mean
  else fc
}

#' Seasonal naive forecast wrapper
#' @param x 'ts' data
#' @param h forecast horizon
#' @param level confidence level. Default is 0.95
#' @param onlyfc if TRUE return only forecasts, otherwise returns a full forecast classed object
#' @param ... extra args, if needed.
#' @return forecasts from ts data or an object of class 'forecast'
#' @examples
#' fit <- sfc_naive(AirPassengers, h=10, onlyfc = FALSE)
#' plot(fit)
#' Mresid(fit)
#' tsSummary(fit)
#' @export
sfc_naive <- function(x, h, level = 95, onlyfc=TRUE, ...) {
  fc <- forecast::snaive(y = x, h = h, level=level, ...)
  if (onlyfc) fc$mean
  else fc
}

#' Random walk forecast wrapper
#' @param x 'ts' data
#' @param h forecast horizon
#' @param level confidence level. Default is 0.95
#' @param onlyfc if TRUE return only forecasts, otherwise returns a full forecast classed object
#' @param ... extra args, if needed.
#' @return forecasts from ts data or an object of class 'forecast'
#' @examples
#' fit <- fc_rw(AirPassengers, h=10, onlyfc = FALSE)
#' plot(fit)
#' Mresid(fit)
#' tsSummary(fit)
#' @export
fc_rw <- function(x, h, level=95, onlyfc=TRUE, ...) {
  fc <- forecast::rwf(y = x, h = h, level=level, ...)
  if (onlyfc) fc$mean
  else fc
}

#' Theta forecast wrapper
#' @param x 'ts' data
#' @param h forecast horizon
#' @param level confidence level. Default is 0.95
#' @param onlyfc if TRUE return only forecasts, otherwise returns a full forecast classed object
#' @param ... extra args, if needed.
#' @return forecasts from ts data or an object of class 'forecast'
#' @examples
#' fit <- fc_theta(AirPassengers, h=10, onlyfc = FALSE)
#' plot(fit)
#' Mresid(fit)
#' tsSummary(fit)
#' @export
fc_theta <- function(x,h,level=95, onlyfc=TRUE, ...) {
  fc <- forecast::thetaf(y = x, h = h, level=level, ...)
  if (onlyfc) fc$mean
  else fc
}

#' Linear model forecast wrapper
#' @param x 'ts' data
#' @param h forecast horizon
#' @param level confidence level. Default is 0.95
#' @param onlyfc if TRUE return only forecasts, otherwise returns a full forecast classed object
#' @param xreg in liner model with covariates it's must be provide
#' @param newxreg in liner model with covariates it's must be provide for forecasting the covariates
#' @param ... extra args, if needed.
#' @return forecasts from ts data or an object of class 'forecast'
#' @examples
#' fit <- fc_lm(AirPassengers, h=10, onlyfc = FALSE)
#' plot(fit)
#' Mresid(fit)
#' tsSummary(fit)
#' @export
fc_lm <- function(x, h, level=95, onlyfc=TRUE, xreg=NULL,newxreg=NULL,...) {
  x <- data.frame(x)
  colnames(x) <- 'x'
  if (is.null(xreg) & is.null(newxreg)) {
    fit <- tslm(x ~ trend + season, data=x, ...)
	fc <- forecast(fit, h=h, level=level)
	if (onlyfc) {
		fc$mean
	}  else {
		fc
	}

  } else if ((!is.null(xreg)) & !(is.null(newxreg))) {
    newnames <- c('x',colnames(xreg))
    x <- cbind(x,xreg)
    colnames(x) <- newnames
    fmla <- as.formula(paste("x ~ trend + season +", paste(colnames(xreg), collapse= "+")))
    fit <- tslm(fmla, data=x, ...)

	fc <- forecast(fit, h=h, level=level, newdata=newxreg)
	fc_mean <- fc$mean
	if (onlyfc) {
		fc$mean
	} else {
		fc
	}
  } else {
    stop('xreg and newxreg must both be NULL or both be provided')
  }
}

#' Structural time series forecast wrapper
#' @param x 'ts' data
#' @param h forecast horizon
#' @param level confidence level. Default is 0.95
#' @param onlyfc if TRUE return only forecasts, otherwise returns a full forecast classed object
#' @param ... extra args, if needed.
#' @return forecasts from ts data or an object of class 'forecast'
#' @examples
#' fit <- fc_sts(AirPassengers, h=10, onlyfc = FALSE)
#' plot(fit)
#' Mresid(fit)
#' tsSummary(fit)
#' @export
fc_sts <- function(x, h, level=95,  onlyfc=TRUE, ...) {
  fit <- StructTS(x, ...)
  fc <- forecast::forecast(fit, h=h, level=level)
  if (onlyfc) fc$mean
  else fc
}

#' Stl forecast wrapper
#' @param x 'ts' data
#' @param h forecast horizon
#' @param level confidence level. Default is 0.95
#' @param method Method to use for forecasting the seasonally adjusted series.
#' @param onlyfc if TRUE return only forecasts, otherwise returns a full forecast classed object
#' @param ... extra args, if needed.
#' @return forecasts from ts data or an object of class 'forecast'
#' @examples
#' fit <- fc_stl(AirPassengers, h=10, onlyfc = FALSE)
#' plot(fit)
#' Mresid(fit)
#' tsSummary(fit)
#' @export
fc_stl <- function(x, h, level=95, method='arima', onlyfc=TRUE, ...) {
  fc <- forecast::stlf(y = x, h = h, method = method, level=level, ...)
  if (onlyfc) fc$mean
  else fc
}

#' auto.arima forecast wrapper
#' @param x 'ts' data
#' @param h forecast horizon
#' @param level confidence level. Default is 0.95
#' @param onlyfc if TRUE return only forecasts, otherwise returns a full forecast classed object
#' @param xreg in ARIMA model with covariates it's must be provide
#' @param newxreg in ARIMA model with covariates it's must be provide for forecasting the covariates
#' @param ... extra args, if needed.
#' @return forecasts from ts data or an object of class 'forecast'
#' @examples
#' fit <- fc_auto.arima(AirPassengers, h=10, onlyfc = FALSE)
#' plot(fit)
#' Mresid(fit)
#' tsSummary(fit)
#' @export
fc_auto.arima <- function(x, h, level=95, onlyfc=TRUE, xreg=NULL,newxreg=NULL,...) {
  fit <- forecast::auto.arima(y = x, xreg = xreg, ...)
  fc <- forecast::forecast(fit, h=h, level=level, xreg=newxreg)
  if (onlyfc) fc$mean
  else fc
}

#' Ets forecast wrapper
#' @param x 'ts' data
#' @param h forecast horizon
#' @param level confidence level. Default is 0.95
#' @param onlyfc if TRUE return only forecasts, otherwise returns a full forecast classed object
#' @param ... extra args, if needed.
#' @return forecasts from ts data or an object of class 'forecast'
#' @examples
#' fit <- fc_ets(AirPassengers, h=10, onlyfc = FALSE)
#' plot(fit)
#' Mresid(fit)
#' tsSummary(fit)
#' @export
fc_ets <- function(x, h, level=95, onlyfc=TRUE, ...) {
  fit <- forecast::ets(y = x, ...)
  fc <- forecast::forecast(fit, h=h, level=level)
  if (onlyfc) fc$mean
  else fc
}

#' BATS forecast wrapper
#' @param x 'ts' data
#' @param h forecast horizon
#' @param level confidence level. Default is 0.95
#' @param onlyfc if TRUE return only forecasts, otherwise returns a full forecast classed object
#' @param ... extra args, if needed.
#' @return forecasts from ts data or an object of class 'forecast'
#' @examples
#' fit <- fc_bats(AirPassengers, h=10, onlyfc = FALSE)
#' plot(fit)
#' Mresid(fit)
#' tsSummary(fit)
#' @export
fc_bats <- function(x, h, level=95, onlyfc=TRUE, ...) {
  fit <- forecast::bats(y =x, ...)
  fc <- forecast::forecast(fit, h=h, level=level)
  if (onlyfc) fc$mean
  else fc
}

#' TBATS forecast wrapper
#' @param x 'ts' data
#' @param h forecast horizon
#' @param level confidence level. Default is 0.95
#' @param onlyfc if TRUE return only forecasts, otherwise returns a full forecast classed object
#' @param ... extra args, if needed.
#' @return forecasts from ts data or an object of class 'forecast'
#' @examples
#' fit <- fc_tbats(AirPassengers, h=10, onlyfc = FALSE)
#' plot(fit)
#' Mresid(fit)
#' tsSummary(fit)
#' @export
fc_tbats <- function(x, h, level=95, onlyfc=TRUE, ...) {
  fit <- forecast::tbats(y = x, ...)
  fc <- forecast::forecast(fit, h=h, level=level)
  if (onlyfc) fc$mean
  else fc
}

#' NNetar forecast wrapper
#' @param x 'ts' data
#' @param h forecast horizon
#' @param level confidence level. Default is 0.95
#' @param onlyfc if TRUE return only forecasts, otherwise returns a full forecast classed object
#' @param nn_p Embedding dimension for non-seasonal time series. Number of non-seasonal lags used as inputs. For non-seasonal time series, the default is the optimal number of lags (according to the AIC) for a linear AR(p) model. For seasonal time series, the same method is used but applied to seasonally adjusted data (from an stl decomposition) - From Rob J Hyndman.
#' @param ... extra args, if needed.
#' @return forecasts from ts data or an object of class 'forecast'
#' @examples
#' fit <- fc_nnetar(AirPassengers, h=10, onlyfc = FALSE)
#' plot(fit)
#' Mresid(fit)
#' tsSummary(fit)
#' @export
fc_nnetar <- function(x, h, level=95,  onlyfc=TRUE, nn_p=1, ...) {
  fit <- forecast::nnetar(y = x, p = nn_p, ...)
  fc <- forecast::forecast(fit, h=h, level=level)
  if (onlyfc) fc$mean
  else fc
}

#' ses forecast wrapper
#' @param x 'ts' data
#' @param h forecast horizon
#' @param level confidence level. Default is 0.95
#' @param onlyfc if TRUE return only forecasts, otherwise returns a full forecast classed object
#' @param ... extra args, if needed.
#' @return forecasts from ts data or an object of class 'forecast'
#' @examples
#' fit <- fc_ses(AirPassengers, h=10, onlyfc = FALSE)
#' plot(fit)
#' Mresid(fit)
#' tsSummary(fit)
#' @export
fc_ses <- function(x, h, level=95, onlyfc=TRUE, ...) {
  fit <- forecast::ses(y = x, h = h, level = level, ...)
  fc  <- forecast::forecast(fit, h = h, level = level)
  if (onlyfc) fc$mean
  else fc
}

#' Holt forecast wrapper
#' @param x 'ts' data
#' @param h forecast horizon
#' @param level confidence level. Default is 0.95
#' @param onlyfc if TRUE return only forecasts, otherwise returns a full forecast classed object
#' @param ... extra args, if needed.
#' @return forecasts from ts data or an object of class 'forecast'
#' @examples
#' fit <- fc_holt(AirPassengers, h=10, onlyfc = FALSE)
#' plot(fit)
#' Mresid(fit)
#' tsSummary(fit)
#' @export
fc_holt <- function(x, h, level=95,  onlyfc=TRUE, ...) {
  fit <- forecast::holt(y = x, h = h, level=level, ...)
  fc  <- forecast::forecast(fit, h=h, level=level)
  if (onlyfc) fc$mean
  else fc
}

#' HoltWinters forecast wrapper
#' @param x 'ts' data
#' @param h forecast horizon
#' @param level confidence level. Default is 0.95
#' @param onlyfc if TRUE return only forecasts, otherwise returns a full forecast classed object
#' @param ... extra args, if needed.
#' @return forecasts from ts data or an object of class 'forecast'
#' #@examples
#' #fit <- fc_hw(lynx, h=10, onlyfc = FALSE)
#' #plot(fit)
#' #Mresid(fit)
#' #tsSummary(fit)
#' #@export
fc_hw <- function(x, h, level=95,  onlyfc=TRUE, ...) {
	## Check Seasonality, if YES, use Holt, else use HW
	if (sum(cycle(x)) == length(x)) {
		fit <- forecast::holt(y = x, h = h, level = level, ...)
	} else {
		fit <- forecast::hw(y = x, h = h, level = level, ...)
	}

  fc  <- forecast::forecast(fit, h=h, level=level)
  if (onlyfc) fc$mean
  else fc
}



#' Seasolaity detection test
#' 
#' This function perform a simple test for seasonality test based on http://robjhyndman.com/hyndsight/detecting-seasonality/.
#' The idea is to perform a chi-squared test for the deviance on two forecst models of type ETS. The first model is an automated ETS model and the second is a ETS model with fixed seasonal component. The NULL is that the data has a seasonal component. If chisq p.value is greater than (1-level) then reject NULL.
#' 
#' @note Due this function performs two ETS forecast inside it. It can be slow.
#' 
#' @param x 'ts' data to be tested
#' @param level significance value to test the null hypothesis that the data has a seasonal pattern.
#' @return bolean TRUE or FALSE
#' @examples 
#' par(mfrow = c(2,1))
#' plot(AirPassengers)
#' plot(lynx)
#' # don't run
#' #DetectSeasonality(AirPassengers)
#' #DetectSeasonality(lynx)
#' par(mfrow = c(1,1))
#' @author LOPES, J. E.
#' @export
#' 
DetectSeasonality <- function(x, level = 0.95){
  f0 <- ets(x)
  f1 <- ets(x, model = "ANN")
  de <- 2*c(logLik(f0) - logLik(f1))
  df <- attributes(logLik(f0))$df - attributes(logLik(f1))$df 
  oo <- 1-pchisq(deviance,df)
  return(oo < 1-level)
}

#' Meanf forecast wrapper
#' @param x 'ts' data
#' @param h forecast horizon
#' @param level confidence level. Default is 0.95
#' @param onlyfc if TRUE return only forecasts, otherwise returns a full forecast classed object
#' @param ... extra args, if needed.
#' @return forecasts from ts data or an object of class 'forecast'
#' @examples
#' fit <- fc_mean(AirPassengers, h=10, onlyfc = FALSE)
#' plot(fit)
#' Mresid(fit)
#' tsSummary(fit)
#' @export
fc_mean <- function(x, h, level=95,  onlyfc=TRUE, ...) {
  fit <- forecast::meanf(y = x, h = h, level=level,  ...)
  fc  <- forecast::forecast(fit, h=h, level=level)
  if (onlyfc) fc$mean
  else fc
}


#' HoltWinters Sazonal forecast wrapper by stats
#' @param x 'ts' data
#' @param h forecast horizon
#' @param level confidence level. Default is 0.95
#' @param onlyfc if TRUE return only forecasts, otherwise returns a full forecast classed object
#' @param ... extra args, if needed.
#' @return forecasts from ts data or an object of class 'forecast'
#' @examples
#' fit <- fc_hws(AirPassengers, h=10, onlyfc = FALSE)
#' plot(fit)
#' Mresid(fit)
#' tsSummary(fit)
#' @export
fc_hws <- function(x, h, level=95, onlyfc=TRUE, ...) {
  fit <- HoltWinters(x, ...)
  fc <- forecast::forecast(fit, h=h, level=level)
  if (onlyfc) fc$mean
  else fc
}

#' HoltWinters Non Seazonal forecast wrapper
#' @param x 'ts' data
#' @param h forecast horizon
#' @param level confidence level. Default is 0.95
#' @param onlyfc if TRUE return only forecasts, otherwise returns a full forecast classed object
#' @param ... extra args, if needed.
#' @return forecasts from ts data or an object of class 'forecast'
#' @examples
#' fit <- fc_hwns(AirPassengers, h=10, onlyfc = FALSE)
#' plot(fit)
#' Mresid(fit)
#' tsSummary(fit)
#' @export
fc_hwns <- function(x, h, level=95, onlyfc=TRUE, ...) {
  fit <- HoltWinters(x, gamma = FALSE, ...)
  fc <- forecast::forecast(fit, h=h, level=level)
  if (onlyfc) fc$mean
  else fc
}

#' HoltWinters Exponential Smoothing forecast wrapper
#' @param x 'ts' data
#' @param h forecast horizon
#' @param level confidence level. Default is 0.95
#' @param onlyfc if TRUE return only forecasts, otherwise returns a full forecast classed object
#' @param ... extra args, if needed.
#' @return forecasts from ts data or an object of class 'forecast'
#' @examples
#' fit <- fc_hwes(AirPassengers, h=10, onlyfc = FALSE)
#' plot(fit)
#' Mresid(fit)
#' tsSummary(fit)
#' @export
fc_hwes <- function(x, h, level=95, onlyfc=TRUE, ...) {
  fit <- HoltWinters(x, gamma = FALSE, beta = FALSE, ...)
  fc <- forecast::forecast(fit, h=h, level=level)
  if (onlyfc) fc$mean
  else (fc)
}

#' Switcher of time series forecast methods
#'
#' Auxiliary switcher that helps the forecasting process. If you pass it a forecast wrapper it returns the forecast model.
#' @param x 'ts' data
#' @param nmodelo string with forecast wraper, eg. fc_ets
#' @param h forecast horizon
#' @param level confidence level. Default is 0.95
#' @param onlyfc if TRUE return only forecasts, otherwise returns a full forecast classed object
#' @param ... extra args, if needed.
#' @author LOPES, J. E.
#' @return forecasts from ts data or an object of class 'forecast'
#' @examples
#'fit <- switch.cvforecast(AirPassengers, "fc_ets", h=20, level=95)
#'class(fit)
#'plot(fit)
#'@export
switch.cvforecast <- function(x, nmodelo, h, level=95, onlyfc=FALSE) {
  switch(nmodelo,
		fc_sts    = fc_sts(x, h=h, level=level, onlyfc=onlyfc),
		fc_hw     = fc_hw(x, h=h, level=level, onlyfc=onlyfc),
		fc_tbats  = fc_tbats(x, h=h, level=level),
		fc_auto.arima = fc_auto.arima(x, h=h, level=level, onlyfc=onlyfc),
		fc_stl    = fc_stl(x, h=h, level=level),
		fc_ses    = fc_ses(x, h=h, level=level, onlyfc=onlyfc),
		fc_holt   = fc_holt(x, h=h, level=level, onlyfc=onlyfc),
		fc_bats   = fc_bats(x, h=h, level=level, onlyfc=onlyfc),
		fc_ets    = fc_ets(x, h=h, level=level, onlyfc=onlyfc),
		fc_lm     = fc_lm(x, h=h, level=level, onlyfc=onlyfc),
		fc_theta  = fc_theta(x, h=h, level=level, onlyfc=onlyfc),
		fc_rw     = fc_rw(x, h=h, level=level, onlyfc=onlyfc),
		sfc_naive = sfc_naive(x, h=h, level=level, onlyfc=onlyfc),
		fc_naive  = fc_naive(x, h=h, level=level, onlyfc=onlyfc),
		fc_mean   = fc_mean(x, h=h, level=level, onlyfc=onlyfc),
		fc_nnetar = fc_nnetar(x, h=h, level=level, onlyfc=onlyfc),
		fc_hws    = fc_hws(x, h=h, level=level, onlyfc=onlyfc),
		fc_hwns   = fc_hwns(x, h=h, level=level, onlyfc=onlyfc),
		fc_hwes   = fc_hwes(x, h=h, level=level, onlyfc=onlyfc)
	)
}

#' Box and Cox tests and Ljung
#'
#' Use Mann-Kendall test (MK) and the Seasonal and the Regional Kendall Tests for trend (SKT and RKT) and Theil-Sen's slope estimator for checking trend
#' @param x 'ts' data
#' @param npoints data points to test. Default twelve points.
#' @param ... extra args, if needed.
#' @author LOPES, J. E.
#' @return trend analysis statistics. See \code{\link{rkt}} for more information about tests
#' @examples
#' nparTrend(AirPassengers)
#' @export
nparTrend <- function (x, npoints = 12, ...) {
  #cat("Check trend in the last", npoints, "data points!\n")
  if(length(x) > 12 & npoints < length(x)) x <- x[(length(x)-(npoints-1)):length(x)]

  rk <- rkt::rkt(time(x), x)
  if (round(rk$B, 6) != 0) {
    if (rk$tau > 0.1) {
      if (rk$tau <= 0.95)
        sinal <- 1L
      else sinal <- 2L
    }
    else if (rk$tau < -0.1) {
      if (rk$tau > -0.95)
        sinal <- -1L
      else sinal <- -2L
    }
    else {
      sinal <- 0L
    }
  }
  else {
    sinal <- 0L
  }
  out <- c(slope = rk$B, tau = rk$tau, score = rk$S, p.value = rk$sl,
           trend_sign = sinal)
  round(out, 4)
}


#' Error handler improved
#'
#' Improved verions of R error handler
#' @param code any stuff you want
#' @param silent if TRUE doesn't print results
#' @author LOPES, J. E.
#' @export
Try_error <- function(code, silent = TRUE) {
	W <- NULL
	w.handler <- function(w){
		W <<- w
		invokeRestart("muffleWarning")
	}
	withCallingHandlers(tryCatch(code, error = function(c) {
    msg <- conditionMessage(c)
    if (!silent) message(c)
    invisible(structure(msg, class = "try-error"))
  }), warning = w.handler)
}

#' Linearity tests
#'
#' Compute six linearity tests among Terasvirta, White, Keenan, McleodLi, Tsay and Tar.
#' @seealso See \code{\link{terasvirta.test}}, \code{\link{white.test}},\code{\link{Keenan.test}},\code{\link{Tsay.test}} and \code{\link{tlrt}}
#' @param x 'ts' data
#' @param Test string containing test's names. They are "terasvirta", "white", "keenan", "tsay", "tarTest".
#' @author LOPES, J. E.
#' @return data.frame with tests statistics and p.values
#' @examples
#' linearityTest(AirPassengers)
#' linearityTest(AirPassengers, Test="tsay")
#' @export
linearityTest <- function(x, Test, all = FALSE) {
	if (class(x)[1] != "ts") x <- ts(x)
	if (missing(Test)) Test <- "keenan"
	else {
		Test <- match.arg(Test, c("terasvirta","white", "keenan", "tsay","tarTest"))
	}

	test <- switch(Test,
			terasvirta = tseries::terasvirta.test(x = x, type = "Chisq"),
			white = tseries::white.test(x),
			keenan = TSA::Keenan.test(x),
			#mcleodLi = TSA::McLeod.Li.test(y = x, plot = FALSE),
			tsay = TSA::Tsay.test(x),
			tarTest = TSA::tlrt(x)
	)

	if (Test == "terasvirta") {
		out <- data.frame(statistic = test$statistic, p.value = test$p.value)
	} else if (Test == "white") {
		out <- data.frame(statistic = test$statistic, p.value = test$p.value)
	} else if (Test == "keenan") {
		out <- data.frame(statistic = test$test.stat, p.value = test$p.value)
	} else if (Test == "tsay") {
		out <- data.frame(statistic = test$test.stat, p.value = test$p.value)
	} else if (Test == "tarTest") {
		out <- data.frame(statistic = test$test.stat, p.value = test$p.value)
	} else {
		cat("Nenhum dos testes se aplica!\n\n")
	}

	if(all) {## recursive sapply
	  dd <- do.call("rbind",
	                as.data.frame(sapply(c("terasvirta","white", "keenan", "tsay","tarTest"), 
	                                       function(i, ...) {
	                                         o <- try(linearityTest(x = x, Test = i))
	                                         if(class(o)[1] != "try-error") o else NULL
	                                       }))
	  )
	  return(suppressMessages(suppressWarnings(dd)))
	 }
	
	out <- round(out, 4)
	rownames(out) <- Test
	return(out)
}

#' Generic function to convert time series data
#'
#' Suport several data types included, data.frame(date, value), xts, zoo, ts or numeric. If numeric data are input, R try to build a date pattern for output in case of  \code{outType} is defined as "xts" or "df", otherwise it returns a ts object.
#' @param Data input data set
#' @param tsfrequency data frequency. It can be "year", "month", "day", "hour", "min" or "sec"
#' @param dateformat format for dates in case of data.frame data. See \code{\link{strptime}} for more information and examples
#' @param outType output format os converted data. It can be "ts" "xts" or "df" data.frame
#' @param OutlierClean if TRUE compute outliers clean. See \code{\link{tsclean}}
#' @param tz timezone. Default is the system tz
#' @param ... extra args
#' @author LOPES, J. E.
#' @examples
#' data(datasample)
#' ConvertData(datasample[,1:3], tsfrequency="day", dateformat='%d/%m/%Y %H:%M:%S',outType = "ts")
#' ConvertData(datasample[,1:3], tsfrequency="day", dateformat='%d/%m/%Y %H:%M:%S',outType = "xts")
#' #ConvertData(AirPassengers, tsfrequency = "month")
#' #ConvertData(rnorm(30), tsfrequency = "month")
#' #ConvertData(rnorm(30), tsfrequency = "month", outType = "xts")
#' @export
ConvertData <- function(Data, tsfrequency = "day", dateformat='%d/%m/%Y %H:%M:%S', OutType = "ts", OutlierClean = TRUE, tz = Sys.getenv("TZ"), ...) {

	## Check data, dates and frequencies
	if(is.null(Data)) stop("Empty dataset!\n")
	if(is.null(dateformat) && any(class(Data)=="data.frame")) stop("data.frame needs first column as date format. ex: '%Y/%m/%d', '%Y-%m-%d', '%Y/%m/%d %H:%M:%S', '%Y-%m-%d %H:%M:%S', etc.!\n")
	#check data frequancies
	tsfrequency <- match.arg(tsfrequency, c("year","month","day","hour","min","sec"))

	OutType <- match.arg(OutType, c("ts","xts","df"))

	# check ts frequency
	if (tsfrequency %in% c("year","month")){
	freq <- 12
	} else if (tsfrequency=="day") {
	freq <- 7
	} else if (tsfrequency=="hour") {
	freq <- 24
	} else if (tsfrequency=="min"){
	freq <- 60
	} else {
	freq <- 1
	}
	## Rules for converting process
	  if (any(class(Data)=="data.frame")){
		date  <- as.character(Data[,1])
		value <- data.frame(Data[,-1, drop=FALSE])
		fmt <- dateformat
		#out <- TimeSeries(date, fmt, value)
	} else if(any(class(Data) %in% c("zoo","xts"))) {
		stopifnot(inherits(time(Data), "Date") || inherits(time(Data), "POSIXt"))
		fmt <- if (inherits(time(Data), "Date")) "%Y-%m-%d" else "%Y-%m-%d %H:%M:%S"
		if (length(dim(Data)) < 2) {
			date <- time(Data)
			value <- data.frame(zoo::coredata(Data))
			names(value) <- deparse(substitute(Data))
			#out <- TimeSeries(date, fmt, value)
		} else {
			date <- time(Data)+0.01
			value <- as.data.frame(zoo::coredata(Data))
			if (length(dim(Data)) < 2) names(value) <- deparse(substitute(Data))
			#out <- TimeSeries(date, fmt, value)
		}
	} else if(any(class(Data) %in% c("ts","mts"))) {
		fmt <- if (inherits(time(Data), "Date")) "%Y-%m-%d" else "%Y-%m-%d %H:%M:%S"
		date  <- as.POSIXlt(as.Date(time(Data)))+0.01
		value <- data.frame(zoo::coredata(Data))
		if (length(dim(Data)) < 2) names(value) <- deparse(substitute(Data))
		#out <- TimeSeries(date, fmt, value)
	} else if(any(class(Data) %in% c("numeric","integer"))){
		Data <- ts(Data)
		fmt <- if (inherits(time(Data), "Date")) "%Y-%m-%d" else "%Y-%m-%d %H:%M:%S"
		date  <- as.POSIXlt(as.Date(time(Data)))+0.01
		value <- data.frame(value = zoo::coredata(Data))
		#out <- TimeSeries(date, fmt, value)
	} else {
		stop("Check your dataset, dates and/or class!\n")
	}

	if (OutlierClean) {
		value <- sapply(value, function(X) {
			tmp <- Try_error(tsclean(X))
			if(class(tmp) !="try-error") tmp else X
		})

		value <- as.data.frame(value)
	}

	out <- TimeSeries(date, fmt, value)

	if(OutType=='ts') {
		O <- ts(value, start = Start(tsfrequency, out), frequency=freq)
	} else if (OutType == "xts"){
		O <- xts::as.xts(zoo::zoo(x=out[-seq(7)], out$dates), tzone = tz)
	} else {
		O <- out
	}
	return(invisible(O))
}

#' Generate forecast horizon time sequence
#'
#' Suport several data types included xts, zoo, ts or numeric. If numeric data are input, R try to build a date pattern based on sysdate().
#' @param XtsData input data set
#' @param tsfrequency data frequency. It can be "year", "month", "day", "hour", "min" or "sec"
#' @param horizon number representing total amount of points to generate the sequency of dates for the forecasts
#' @author LOPES, J. E.
#' @return A lista conatining data transformed to ts and a date sequency for the forecasts
#' @examples
#' x <- ForecastHorizon(rnorm(100), 'day',20)
#' str(x)
#' data(datasample, package = "cvforecast")
#' y <- ConvertData(datasample[,1:2], tsfrequency = "day", OutType = "xts")
#' y <- ForecastHorizon(y, 'day', 20)
#' @export
ForecastHorizon <- function(XtsData, tsfrequency, horizon) {

	if(!class(XtsData)[1] %in% c("zoo","xts","ts","numeric")) stop("Data must be of 'zoo', 'xts', 'ts' or 'numeric' classes!\n")

	#check date arg
	tsfrequency <- match.arg(tsfrequency, c("year","month","day","hour","min","sec"))

	# check ts frequency
	if (tsfrequency %in% c("year","month")){freq <- 12
	} else if (tsfrequency=="day") {freq <- 7
	} else if (tsfrequency=="hour") {freq <- 24
	} else {freq <- 1}

	if (class(XtsData)[1] %in% c("ts","numeric")) {

		if (class(XtsData)[1] == "numeric") {
			x <- ts(XtsData, frequency=freq)
			SE <- seq(as.POSIXct(Sys.Date(), tz = "GMT"), by = tsfrequency, l=horizon)
		} else {
		  x <- XtsData
		  SE <- seq(max(as.POSIXct(as.Date(time(x))+1,tz = "GMT")), by = tsfrequency, l=horizon)
		}
		return(list(x, FCHorizon = SE))
	} else {
		IndexXTS <- index(XtsData)
		IndexXTS <- IndexXTS[!is.na(IndexXTS)] # Remove missing cado haja
		SE <- seq(from=max(IndexXTS), by=tsfrequency, length.out = horizon+1)[-1]
		value  <- as.numeric(XtsData)

		tdata <- TimeSeries(as.character(IndexXTS), "%Y-%m-%d", value)

		if(class(XtsData)[1] == "xts") {
			x <- ts(value, start = Start(tsfrequency, tdata), frequency = freq)
		} else {
			x <- XtsData
		}
	}
	return(list(x, FCHorizon = SE))
}

#' Compute starting values for ts transformation
#'
#' This switcher is a internal function that chooses frequencies based on char dates. It's used to help transforming data into ts objects in R
#' @param freq data frequency for the dates. It works with  "year", "month", "day", "hour" and "min".
#' @param tdata data came from \code{TimeSeries} function
#' @author LOPES, J. E.
#' @export
Start <- function(freq, tdata) {
  if (!freq %in% c("year","month","day","hour","min")) {
	return(c(1))
  } else {
	  switch(freq,
			year  = c(tdata$year[1], tdata$month[1]),
			month = c(tdata$year[1], tdata$month[1]),
			day   = c(tdata$week[1], tdata$day[1]),
			hour  = c(tdata$day[1],  tdata$hour[1]),
			min   = c(tdata$hour[1], tdata$minute[1])
			)
	   }
}

#' @title PRESS and other statistics from lm model
#'
#' @description Returns the PRESS statistic (predictive residual sum of squares).
#' Useful for evaluating predictive power of regression models.
#' @param obj A linear regression model (class 'lm'). Required.
#' @author LOPES, J. E.
#' @examples
#' fit <- step(lm(mpg~., data = mtcars), trace = 0)
#' lm_press_stat(fit)
#' @export
lm_press_stat <- function(obj) {
	if(class(obj)[1]!="lm") stop(paste("Only lm models are allowed!"))
	PRESS <- function(obj) {
	  #' calculate the predictive residuals
	  pr <- residuals(obj)/(1-lm.influence(obj)$hat)
	  #' calculate the PRESS
	  PRESS <- sum(pr^2)
	  return(PRESS)
	}
	pred_r_squared <- function(obj) {
	  #' Use anova() to get the sum of squares for the linear model
	  lm.anova <- anova(obj)
	  #' Calculate the total sum of squares
	  tss <- sum(lm.anova$'Sum Sq')
	  # Calculate the predictive R^2
	  pred.r.squared <- 1-PRESS(obj)/(tss)

	  return(pred.r.squared)
	}

	r.sqr <- summary(obj)$r.squared
	adj.r.sqr <- summary(obj)$adj.r.squared
	pre.r.sqr <- pred_r_squared(obj)
	PRESS <- PRESS(obj)

	return.df <- data.frame(press.r.squared = pre.r.sqr, press = PRESS, r.squared = r.sqr, adj.r.squared = adj.r.sqr)
	return(return.df)
}


#' Convert char dates in Dates keeping data
#'
#' Convert a data.frame with dates and numeric values into a new one including year, month, day, minute and second generated trough the \code{dateformat} applied to the original dates.
#' @param dates char vector of dates same length of data
#' @param dateformat date format. See \code{\link{strptime}} for more information and examples
#' @param data data numeric vector, matrix or data.frame containing time series data
#' @param tz time zone. It gets the system date
#' @return data.frame with split dates and values
#' @examples
#' dates <- as.character(as.POSIXct(as.Date(time(AirPassengers))+1))
#' value <- as.numeric(AirPassengers)
#' dateformat <- "%Y-%m-%d"
#' y <- TimeSeries(dates, dateformat, value)
#' @export
TimeSeries<-function(dates, dateformat, data = NULL, tz = "GMT")
{
  dimensions<-dim(data)
  clss<-lapply(data,class)
  if (!is.null(data)){
    if (!is.null(dimensions[1])){
      if (length(dates)!=length(data[,1]))
        stop("Data and Date lengths differ")
    }else{
      if (length(dates)!=length(data))
        stop("Data and Date lengths differ")
    }
  }
  dates <- (strptime(paste(dates), dateformat, tz = tz))
  minute <- minute(dates)
  hour <- hour(dates)
  day <- day(dates)
  week <- week(dates)
  month <- month(dates)
  year <- year(dates)
  if (is.null(data)) {
    results <- data.frame(dates, minute, hour, day, week,
                          month, year)
  }
  else {
    results <- data.frame(dates, minute, hour, day, week,
                          month, year, data)
  }
  return(results)
}

#' Days Extra function
#' @export
daysAgg<-function(data,process,multiple=NULL,na.rm=FALSE){
		if(is.null(multiple)){
		multiple=1
	}
		if(multiple==1){
		day<-aggregate(data[,8:length(data)],list(day=data$day,month=data$month,year=data$year),process,na.rm=na.rm)
	days<- ymd(paste(day$year, day$month, day$day))
	data2<-data.frame(date=days,data=day[,4:length(day)])
	names(data2)<-c("Date",names(data[8:length(data)]))
	return(data2)
		}
	temp<-data
	day<-aggregate(list(data[,8:length(data)],count=1),list(day=data$day,month=data$month,year=data$year),process,na.rm=na.rm)
	days<- ymd(paste(day$year, day$month, day$day))
	data<-data.frame(date=days,day[,5:length(day)-1],count=day[length(day)])
	days=paste(multiple,"days")
	all.dates<-seq.Date(as.Date(data$date[1]),as.Date(data$date[length(data[,1])]),by="day")
	dates<-data.frame(date=all.dates)
	aggreGated<-merge(dates,data,by="date",all.x=TRUE)
	aggreGated$date<-rep(seq.Date(as.Date(data$date[1]),as.Date(data$date[length(data[,1])]),by=days),each=multiple,length=length(all.dates))
#	data<-subset(aggreGated,!is.na(count))
	results<-aggregate(list(aggreGated[2:length(aggreGated)]),list(date=aggreGated$date),process,na.rm=TRUE)
	results<-subset(results,results$count!=0)
	results<-results[,-length(results)]
	names(results)<-c("Date",names(temp[8:length(temp)]))
	return(results)
}

#' Hours Extra function
#' @export
hoursAgg<-function
(data, process, multiple = 1, na.rm = FALSE, tz = "GMT")
{
	gap <- as.numeric(difftime(data$dates,data$dates[1],tz="GMT", units = "hours"))
	agg.gap<-gap-(gap%%multiple)
	data$dates<-data$dates[1] + agg.gap * 60 * 60
	data <- TimeSeries(data$dates, "%Y-%m-%d %H:%M:%S",
                     data[, 8:length(data)], tz = tz)
	result <- aggregate(data[, 8:length(data)],
                      list(day = data$day,month = data$month,
                           year = data$year, hour = data$hour,
                           minute=data$minute),
                      process, na.rm = na.rm)
  data2 <- data.frame(Date = strptime(
    paste(result$minute, result$hour,result$day, result$month, result$year),
    "%M %H %d %m %Y",tz = tz), result[, 6:length(result)])
  sorted <- data2[order(data2$Date), ]
  final <- data.frame(Date = as.POSIXlt(sorted$Date),
                      sorted[, 2:length(sorted)])
  names(final) <- c("Date", names(data[8:length(data)]))
  return(final)
}

#' Months Extra function
#' @export
monthsAgg<-function(data,
		    process,
	    	    multiple=NULL,
		    na.rm=FALSE){

	if(is.null(multiple)){
		multiple=1
	}

	d.cols<-length(data)
	month<-aggregate(data[,8:d.cols],list(month=data$month,year=data$year),process,na.rm=na.rm)
	data.cols<-length(month)


	if(multiple>1){

month.gap<-month[,1]
for(i in 1:length(month[,1])){
	month.gap[i]=(month[i,2]%%month[1,2])*12+month[i,1]}
month.gap<-month.gap-month.gap%%multiple
month.gap<-month.gap%%12+1

year<-month[,2]
if(sum(month.gap)==length(month.gap)){
year<-year[1]+(year-year[1])-(year-year[1])%%(multiple/12)
}else{

for(i in 2:length(month.gap)){
	if(month.gap[i]==month.gap[i-1])
		year[i]=year[i-1]
	else
		year[i]=month[i,2]
		}
	}
		date=strptime(paste(01,month.gap,year),"%d %m %Y")
		results<-data.frame(date,data=month[,3:data.cols])
		final<-aggregate(results[,2:length(results)],list(date=results$date),process,na.rm=na.rm)
		names(final)<-c("Date",names(data[8:length(data)]))
		return(final)
}
	else{
		date=strptime(paste(01,month$month,month$year),"%d %m %Y")
		results<-data.frame(date,data=month[,3:data.cols])
		names(results)<-c("Date",names(data[8:length(data)]))
		return(results)
	}
}

#' Years Extra function
#' @export
yearsAgg<-function(data,
		    process,
	    	    multiple=NULL,
		    na.rm=FALSE,
		    from.first.obs=TRUE){

	if(is.null(multiple)){
		multiple=1
	}
	#Test the last gap to make sure it's complete.
	start<-min(data$year)
	end<-max(data$year)
	if((end-start)%%multiple!=0)warning("Last Gap does not contain ",multiple," years. There is only ",((end-start)%%multiple)," year(s) in this gap.",call.=FALSE)

	if(from.first.obs==TRUE){
	years<-(as.numeric(difftime(data$date,data$date[1],units=c("hours")))-as.numeric(difftime(data$date,data$date[1],units=c("hours")))%%8765.81277)/8765.81277
	#8765.81277 hours in each year.
	years.again<-years-years%%multiple
	data$year<-data$year[1]+years.again
	date.1<-as.Date(strptime(paste(data$day[1],data$month[1],data$year),"%d %m %Y"))
	new.1<-data.frame(date=date.1,data[,8:length(data)])
	result<-aggregate(new.1[,2:length(new.1)],list(date=new.1$date),process,na.rm=na.rm)
	names(result)<-c("Date",names(data[8:length(data)]))
	}else{
		years<-data$year-data$year[1]
		years.again<-years-years%%multiple
	data$year<-data$year[1]+years.again
	date.1<-as.Date(strptime(paste(1,1,data$year),"%d %m %Y"))
	new.1<-data.frame(date=date.1,data[,8:length(data)])
	result<-aggregate(new.1[,2:length(new.1)],list(date=new.1$date),process,na.rm=na.rm)
	names(result)<-c("Date",names(data[8:length(data)]))
		}



		#Everything that need to be done to calculate the annual agg.
	return(result)
		}

#' Extra function
#'
#' TimeSeries to zoo
#' @export
TimeSeries2zoo <- function(x, ...) {
	zoo(x[-seq(7)], x$dates)
}

#' Extra function
#'
#' Zoo to TimeSeries
#' @export
zoo2TimeSeries <- function(x, ...) {
	stopifnot(inherits(time(x), "Date") || inherits(time(x), "POSIXt"))
	fmt <- if (inherits(time(x), "Date")) "%Y-%m-%d" else "%Y-%m-%d %H:%M:%S"
	if (length(dim(x)) < 2) {
		DF <- data.frame(coredata(x))
		names(DF) <- deparse(substitute(x))
		TimeSeries(time(x), fmt, DF)
	} else TimeSeries(time(x), fmt, as.data.frame(coredata(x)))
}

#' Best forecast among several ones
#'
#' This function receives a list of forecast models from the class 'cvforecst' and returns a list of their characteristics.
#' @param objcvfore 'cvforecst' object
#' @param cvMethod gof statistic, MAPE, MASE, MAE, etc. \code{\link{accuracy}}
#' @param residlevel confidence level for residuals tests
#' @author LOPES, J. E.
#' @examples
#' #Define cross validation parameters
#' myControl <- cvForecastControl(
#' minObs = 14,
#' stepSize = 1,
#' maxHorizon = 30,
#' residlevel = 0.05)
#'
#' #Paralell execution improves the processing time
#' #cl <- makeCluster(4, type='SOCK')
#' #registerDoParallel(cl)
#'
#' #Load data
#' x <- AirPassengers
#' fit <- cvforecast(x, myControl)
#' cvBestModel(fit)
#' @export

cvBestModel <- function(objcvfore, cvMethod = "MAPE", residlevel = 0.10, ...) {

	#if(!class(objcvfore) %in% c("forecast","cvforecast")) stop("Objeto deve ser de classe 'cvforecast'")

	## Análise de residuos dos modelos
    Resid <- try(plyr::ldply(objcvfore, function(X) {
    if (class(X)[1] != "try-error") {
      sm <- sum(Mresid(X) > residlevel)
      sm[sm %in% c(Inf, -Inf, NA, NULL)] <- 99
      sm
    } else NULL
	}))

	if (class(Resid) != "try-error") {
		# Remove entradas missing (se houver)
		Resid <- Filter(Negate(function(X) is.null(unlist(X))), Resid)
		names(Resid) <- c(".id","RESID.PVALUE")

		# Cria rank de resíduo
		Resid$RES.RANK <- rank(Resid$RESID.PVALUE, na.last = TRUE, ties.method = "first")
		Resid <- Resid[,-2]
	} else {
		Resid <- data.frame(.id = NA, RES.RANK = NA)
	}
	## Estatisticas de bondade
	STATS <- Try_error(plyr::ldply(objcvfore, function(X) {
		if (class(X)[1] != "try-error") {
			st <- round(tsSummary(X), 4)
			st[st %in% c(Inf, -Inf, NA, NULL)] <- 999999999
			st
		} else 999999999 # Trata casos de erros com 999.999.999
	}))
	if(class(STATS) != "try-error") {
		STATS <- Filter(Negate(function(X) is.null(unlist(X))), STATS)
	    ## Cria rank de estatística de bondade
		STATS$RESID.GOF <- rank(STATS[,cvMethod], na.last = TRUE, ties.method = "first")
	} else {
		STATS <- data.frame(.id=NA, ME=NA, RMSE=NA, MAE=NA, MPE=NA, MAPE=NA, MASE=NA, ACF1=NA, RESID.GOF=NA)
	}

	if (!all(is.na(STATS)) & !all(is.na(Resid))) {
		ESTFINAL <- merge(STATS, Resid, by.x = ".id")
		ESTFINAL$SOMA.RK <- rowSums(ESTFINAL[,c("RESID.GOF", "RES.RANK")], na.rm = FALSE)
		ESTFINAL <- ESTFINAL[order(ESTFINAL$RESID.GOF),]
		CV_names <- ESTFINAL[,".id"]
	} else {
		# Without residual tests
		ESTFINAL <- STATS[order(STATS$RESID.GOF),]
	}
	return(ESTFINAL)
}

#' Descriptive statistics
#'
#' Receives a data.frame, vector or matrix and computes descriptive statistics.
#' @param dados data.frame, matrix or numeric vector
#' @param nivel significance level for bootstrap non-paramentric ci for the mean
#' @param tipoci type of confidence interval test. See \code{\link{boot.ci}}
#' @param nsimu bootstrap simulations number
#' @param dig number of digits for output
#' @author LOPES, J. E.
#' @importFrom boot boot boot.ci
#' @return data.frame with min, mean, median, max, sd, sd + mean, coefficeint of variation and bootstrap non-parametric confidence interval for the mean
#' @examples
#' Desc(cars)
#' Desc(rnorm(100, 2, 3))
#' @export
Desc <- function(dados, nivel=0.95, tipoci = "basic", nsimu = 500, dig=2) {

 tipoci <- match.arg(tipoci, c("norm","basic", "stud", "perc", "bca"))
 nivel <- nivel

  aux <- function(x,...){
    x <- na.omit(x)
    minimo  <- media <- mediana <- maximo <- desviop <- meddesv <- na <- null <- B <- cv <- as.numeric(NA)
    minimo  <- min(x, na.rm=TRUE)
    media   <- mean(x, na.rm=TRUE)
    mediana <- median(x, na.rm=TRUE)
    maximo  <- max(x, na.rm=TRUE)
    desviop <- sd(x, na.rm=TRUE)
    meddesv <- media+desviop
    na      <- x[is.na(x)]
    null    <- x[is.null(x)]
    cv <- desviop/media

	out <- as.data.frame(t(round(c(Mediana = mediana, DesvP = desviop, "Media+DesvP" = meddesv, Nulos = na + null, Min = minimo, Max = maximo, CV = cv), dig)))

	# Funcao para intervalode confiança da média
    cimean<- function(x, i, ...) {
		m <- mean(x[i])
		n <- length(i)
		v <- (n-1)*var(x[i])/n^2
		c(m, v)
    }
	da.boot <- boot(x, cimean, R = nsimu)

	# A função sink permite omitir textos de print() e também de cat()
	sink(tempfile())
	    B <- boot.ci(da.boot, conf = nivel, type = tipoci)
	sink()
	if (!is.null(B)){
		if (tipoci == "norm") {Li <- B$normal[2];Ls <- B$normal[3]
		} else if(tipoci == "basic") {Li <- B$basic[4]; Ls <- B$basic[5]
		} else if(tipoci == "stud") {Li <- B$student[4]; Ls <- B$student[5]
		} else if(tipoci == "perc") {Li <- B$percent[4]; Ls <- B$percent[5]
		} else {Li <- B$bca[4]; Ls <- B$bca[5]
		}
	} else {Li <- Ls <- media}

    media <- paste(round(media,dig), "(", round(Li,dig), ";", round(Ls,dig),")", sep="")

	nm <- c(paste("Media(IC", 100*nivel, "%)", sep=""), names(out))

	out <- as.data.frame(cbind(media, out))
	colnames(out) <- nm
    return(out)
  }
  dados <- as.data.frame(dados)
  out <- do.call("rbind",
	lapply(dados, function(X) {if(!is.numeric(X)) NULL else aux(X)}))
  return(out)
}


#' Aggregate xts object
#'
#' Aggregate xts object by a given function. It can be any function defined by the user. Actually it's possible to aggregate data by day, month, week, quarter and year.
#' @param daxts object of calss \code{xts} with simple or multiple columns
#' @param FUN a function defined by the user.
#' @param freq a char string to the aggregation. It can be daily, weekly, monthly, quarterly or yearly
#' @param dig number of output digits
#' @author LOPES, J. E.
#' @return a object from class \code{xts} with the data aggregated
#' @seealso \code{\link{apply.daily}}
#' @examples
#' data(datasample)
#' d1 <- cvforecast::ConvertData(datasample,  tsfrequency = "day", OutType = "xts", OutlierClean = FALSE)
#' fun1 <- function(x) {c(mean(x, na.rm=TRUE) + sd(x,na.rm=TRUE))}
#' aggreg(d1[,1:5], fun1, "monthly")
#' aggreg(d1[,1:5], mean, "monthly")
#'
#' @export
aggreg <- function(daxts, FUN, freq = "daily", dig = 4, ...) {
	stopifnot(require(xts))
	freq <- match.arg(freq, c("hourly","daily","weekly","monthly","quarterly","yearly"))

	apply.hourly <- function(x, FUN, ...) {
		ends <- endpoints(x, 'hours', 1)
		period.apply(x, ends, FUN=FUN)
	}

	fun <- function(freq, ...) {
	  switch(freq,
		daily     = apply.daily,
		weekly    = apply.weekly,
		monthly   = apply.monthly,
		quarterly = apply.quarterly,
		yearly    = apply.yearly,
		hourly    = apply.hourly)
	}

	Fun <- fun(freq)
	La <- lapply(daxts, function(X) Fun(X, FUN=FUN))
	Dc <- do.call("cbind.xts", La)
	return(round(Dc, dig))
}
