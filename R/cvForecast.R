#' Core function to compute multiple forecasts by Cross-Validation
#'
#' This the core function of the package. It computes multiple forecasts by the technique of Cross-Validation. The decision about the best models are based on tests as linearity, trend and fit accuracy.
#'
#' @param tsdata data.frame type date-value, ts, mts or xts time series objects
#' @param tsControl generic contol with several args for the modelling process. See
#'    \code{\link{cvForecastControl}}.
#' @param fcMethod accept the forecast method fefined by the user. This argument can be a string or a list, eg. fcMethod = "fc_ets" or a list as fcMethod = list("fc_ets", "fc_hws"). If NULL, decision is made automatically.
#' @param ... other arguments
#' @author LOPES, J. E.
#' @return A list of class 'cvforecast' containing several objetcts from the forecasting process. It includes: betso models (less tahn 6), crossValidation statistics for all models, accuracy of all models, the control, etc. As below.
#' @keywords crossValidation time series
#' @examples
#'#Define cross validation parameters
#'myControl <- cvForecastControl(
#'  minObs = 14,
#'  stepSize = 10,
#'  maxHorizon = 30,
#'  summaryFunc=tsSummary,
#'  cvMethod="MAPE",
#'  tsfrequency='day',
#'  OutlierClean=FALSE,
#'  dateformat='%d/%m/%Y %H:%M:%S')
#'
#'#Paralell execution improves the processing time
#'#require(doParallel)
#' #cl <- makeCluster(4, type='SOCK')
#' #registerDoParallel(cl)
#'
#'#Load data
#'require(plyr)
#'data(datasample, package="cvforecast")
#'dadosd <- ConvertData(datasample[,1:6], dateformat='%d/%m/%Y %H:%M:%S', tsfrequency = "day", OutType="ts")
#'table(sapply(dadosd, class))
#'dim(dadosd)
#'
#'#Looping example
#'
#'FF <- llply(dadosd[,1:2], function(X) {
#'  fit <- try(cvforecast(X, myControl))
#'  if(class(fit) != "try-error") {
#'    #plot(fit)
#'    return(fit)
#'  } else NA
#'}, .progress = "time")
#'
#'table(sapply(FF, class))
#'plot(FF[[1]])
#'sapply(FF, names)
#'#stopCluster(cl)
#'@export
cvforecast <- function(tsdata, tsControl=cvForecastControl(), fcMethod = NULL, ...) {
	h <- tsControl$maxHorizon
	cvMethod <- toupper(tsControl$cvMethod)
	tsfrequency  <- tsControl$tsfrequency
	OutlierClean <- tsControl$OutlierClean
	dateformat <- tsControl$dateformat

	if(class(tsdata)[1] %in% c("data.frame")) {
		# To format the forecast output
		flag_forecast_format <- TRUE
		# Work with data and date formats

		x <- ConvertData(tsdata, tsfrequency=tsfrequency, dateformat=dateformat, OutType = "ts", OutlierClean=OutlierClean)

		XTS <- ConvertData(tsdata, tsfrequency=tsfrequency, dateformat=dateformat, OutType = "xts", OutlierClean=OutlierClean)
		x <- x[,1]
		# Dates in nice format for forecast output
		FCH <- ForecastHorizon(XTS, tsfrequency=tsfrequency, horizon=h)
		ForecastDates <- FCH$FCHorizon
		ForecastDates <- ForecastDates[!is.na(ForecastDates)]
	} else if (class(tsdata)[1] %in% c("xts","ts")){
		x <- tsdata
		flag_forecast_format <- FALSE
	} else if (class(tsdata)[1] == "numeric") {
		x <- ts(tsdata, frequency = 1)
		flag_forecast_format <- FALSE
	} else {
		stop("Check data!\n")
	}
	if(is.null(x) | length(as.numeric(x)) < 8) {
		stop(cat("Cross-Validation with", length(as.numeric(x)), "data points? Are you kidding! I need at least 8 data points!!!\n"))
	}
	if (class(x)[1]!='ts'){
		x <- ts(x, frequency=1)
	}
	if (is.null(fcMethod)) { #Automatically
		nmforecast <- Try_error(forecastMethod(x))
		if(length(as.numeric(x)) < 2*max(cycle(x)) | max(cycle(x))==1 | class(nmforecast) == "try-error") {
			cat("Series with less than 2 cycles or non-periodic. Try ARIMA and ETS!\n")
			tsControl$minObs <- 8
			tsControl$stepSize <- 1
			nmforecast <- list("fc_auto.arima","fc_ets")
		}
	} else if (fcMethod=="all"){
		nmforecast <- list("fc_auto.arima","fc_ets", "fc_lm","cf_naive","fc_rw","fc_sts","fc_theta","fc_hwns","fc_hwes","fc_hws",  "scf_naive")
	} else {
		nmforecast <- fcMethod
	}
	if(length(as.numeric(x)) < 2*max(cycle(x)) | max(cycle(x))==1 | class(nmforecast) == "try-error") {
		cat("Series with less than 2 cycles or non-periodic. Try ARIMA and ETS!\n")
		tsControl$minObs <- 8
		tsControl$stepSize <- 1
		nmforecast <- list("fc_auto.arima","fc_ets")
	}
	out <- plyr::llply(nmforecast, function(X, ...) {
		temp <- Try_error(cvts2(x=x, FUN=get(X), tsControl))
			if (class(temp)[1]!="try-error") return(temp)
			else NA
			}
		)
	names(out) <- nmforecast

	# Remove missing values
	out <- out[!is.na(out)]
	out <- Filter(Negate(function(X) is.null(unlist(X))), out)
	## Statistical tables
	STATS <- lapply(out,
		function(X) {
			if (class(X)[1]!="try-error") X$results
			else NA
		}
	)
	STATS <- Filter(Negate(function(X) is.null(unlist(X))), STATS)
	## Cross-Validation statistics extraction
	STATS1 <- lapply(STATS,
		function(X) {
			if (is.null(dim(X)) | mean(as.numeric(X[, cvMethod]), na.rm=TRUE)%in% c(NA, 0, -Inf, Inf)) {
				tsControl$cvMethod <- cvMethod <- "MAE"
				return(as.numeric(X[, cvMethod]))
			} else { # Is MAPE is wrong, try MAE
				return(as.numeric(X[, cvMethod]))
			}
		}
	)
	## Make data.frame and remove missing elements
	estatisticas <- as.data.frame(do.call("cbind", STATS1))
	estatisticas <- estatisticas[sapply(estatisticas, function(X) length(unique(na.omit(X)))) > 1];
	estatisticas <- na.omit(estatisticas)

	## Coefficient of variation plus median for the accuracy statistics. It works better for decision due data variation in the Cross-Validation process.
	CV <- sapply(estatisticas, function(X) {
			cv = sd(X, na.rm=TRUE)/mean(X, na.rm=TRUE)
			#DAMM = Median Absolute Deviation of the median
			#damm = median(abs(X-median(X, na.rm=TRUE)), na.rm=TRUE)
			#cvmed = damm/median(X, na.rm=TRUE)
			##mpcv <- (cv+cvmed-cv*cvmed) ## Union AUB = A+B-AC
			return(abs(median(X, na.rm=TRUE)) + cv)
		}
	)

	## Remove strange values if exists and sort form the best to worst model
	CV <- sort(CV[!CV %in% c(Inf, -Inf, 0, NA)])
	CV_names <- as.list(names(CV))

	best_models <- plyr::llply(CV_names, function(J) {
			temp <- try(switch.cvforecast(x=x, nmodelo=J, h=h, onlyfc=FALSE))
			if(class(temp)!="try-error") temp else NA
		}
	)
	names(best_models) <- CV_names

	best_models <- best_models[!is.na(best_models)]

	Bm <- cvBestModel(best_models, cvMethod=cvMethod, residlevel = 0.10)
	CV_names <- as.list(Bm[,1])

	## If more than 5 models, keep the 6 firsts
	if (length(CV_names) > 5) CV_names <- CV_names[1:6]

	## Reorder model list according with best models
	temp <- c()
	for(i in CV_names) {
		temp[[i]] <- best_models[[i]]
	}
	best_models <- temp
	
	## Work with best format forecast output. This works only if input dates are in char format.
	if(flag_forecast_format) {
		if(!is.na(best_models[[1]][1])) {
			values <- as.data.frame(best_models[[1]])
			names(values) <- tolower(names(values))
			rownames(values) <- c()
		} else if(!is.na(best_models[[2]][1])){
			values <- as.data.frame(best_models[[2]])
			names(values) <- tolower(names(values))
			rownames(values) <- c()
		} else {
			stop("No best models found!\n")
		}
		bestForecast <- data.frame(date = ForecastDates,values)

		dfx <- tsdata
		#dfx$dates <- strptime(dfx$dates, '%d/%m/%Y %H:%M:%S', tz = Sys.getenv("TZ"))
		dfx$li <- dfx$ls <- dfx$forecast <- as.numeric(NA)
		names(dfx) <- c("data", "realizado","li","forecast","ls")

		dfy <- bestForecast
		dfy$date <- strftime(dfy$date, '%d/%m/%Y %H:%M:%S', tz = Sys.getenv("TZ"))
		dfy$realizado <- as.numeric(NA)

		dfy <- data.frame(data=dfy$date, realizado = dfy$realizado, li=dfy$lo.95, forecast=dfy$point.forecast, ls=dfy$hi.95)

		Tt <- rbind(dfx, dfy)
		Tt <- ConvertData(Tt, dateformat='%d/%m/%Y %H:%M:%S', tsfrequency = "day", OutType="xts", OutlierClean = FALSE)

		D <- data.frame(data = index(Tt), Tt)
		row.names(D) <- NULL
		bestForecast <- list(bestForecast, D)
	} else {
		if(!is.na(best_models[[1]][1])) {
			bestForecast <- best_models[[1]]
		} else if(!is.na(best_models[[2]][1])) {
			bestForecast <- best_models[[2]]
		}  else {
			stop("No best models found!\n")
		}
	}
	# Results
	out1 <- c()
	out1$bestForecast <- bestForecast
	out1$melhores <- best_models
	out1$cv_stat  <- estatisticas[, names(best_models), drop=FALSE]
	out1$acuracia <- STATS
	out1$myControl<- tsControl
	class(out1) <- "cvforecast"
	return(invisible(out1))
}


#' Default summary for Cross-validation forecasts
#'
#' Summarizes \code{cvforecast} classed objects. The output is statistics about the best models fited by \code{cvforecast} function. It may include residual tests results and simple accuracy, cross-validation accuracy and a plot for the models.
#' @param obj \code{cvforecast} object
#' @param digits number of digits for the default print output
#' @param plot if TRUE, plot forecast for best models
#' @param accuracy.best if TRUE, compute simple accuracy statistics for the best final models
#' @param cv.accuracy if TRUE, compute accuracy of the cross-validation points from 1 to horizon
#' @author LOPES, J. L.
#' @seealso \code{\link{cvforecast}}, \code{\link{accuracy}}, \code{\link{Mresid}}
#' @export
summary.cvforecast <- function(obj, digits = 4, plot=TRUE, accuracy.best = TRUE, cv.accuracy = FALSE) {
  if(class(obj) != "cvforecast") stop("'cvforecast' required!\n")
  
  res <- ldply(obj$melhores, function(X) {
    tmp <- Try_error(Mresid(X))
    if(class(tmp)!="try-error") tmp else NULL
  })
  cat("-----------------------------------------------------------------------\n")
  cat("Residual Analysis for the best models!\n")
  print(res, digits = digits)
  cat("-----------------------------------------------------------------------\n")
  
  cv.best <- if (accuracy.best) {
    tmp <- ldply(obj$melhores, function(X) {
      tmp <- Try_error(accuracy(X))
      if(class(tmp) !="try-error") tmp else NULL
    })
    cat("-----------------------------------------------------------------------\n")
    cat("Accuracy for the best models!\n")
    print(tmp, digits = digits)
    cat("-----------------------------------------------------------------------\n")
  }

  cv.ac <- if (cv.accuracy) {
    cat("-----------------------------------------------------------------------\n")
    cat("Cross-Validation Accuracy for the best models!\n")
    print(obj$acuracia, digits = digits)
    cat("-----------------------------------------------------------------------\n")
    } else NULL
  
  if(plot) plot(obj)
  return(invisible(list(cv.best=cv.best, cv.ac=cv.ac, residuals = res)))
}

#setMethod("summary", "cvforecast", summary.cvforecast)

#' Default plot for Cross-validation forecast
#'
#' Plot \code{cvforecast} classed objects
#' @param obj \code{cvforecast} object. See \code{\link{cvforecast}}
#' @param ... extra args, if needed
#' @author LOPES, J. E.
#' @export
plot.cvforecast <- function(obj, ...) {
	# Color palette, blues
	mypalette <- c("#F7FBFF","#DEEBF7","#C6DBEF","#9ECAE1","#6BAED6","#4292C6","#2171B5","#08519C","#08306B")

	estatisticas <- obj$cv_stat
	cvMethod <- obj$myControl$cvMethod
	melhores <- obj$melhores
	#check_trend <- obj$tendencia
	linear <- obj$linear

	## Plot CrossValidation
	estatisticas <- estatisticas[, names(melhores), drop=FALSE]
	Plot <- data.frame(id = 1:nrow(estatisticas), estatisticas)
	P <- ts(Plot[,-1, drop=FALSE])

	## Sequencia com impar para gera??o dos gr?ficos
	## Plot forecasts
	if (length(melhores) > 3) {
		lfor <- c(1, 1, 2, 2, 3:(length(melhores)))
		lrep <- rep(length(melhores) + 1, 2)
		if(length(lfor)%%2 != 0) {
			lfor[length(lfor)+1] <- 0
		}
		nf <- layout(matrix(c(lfor, lrep), byrow=TRUE, ncol=2))
	} else if (length(melhores) == 3){
		lfor <- c(1, 1, 2, 2, 3, 3, 4,4)
		nf <- layout(matrix(c(lfor), byrow=TRUE, ncol=2))
	} else if (length(melhores) == 2){
	  lfor <- c(1, 1, 2, 2, 3, 3)
	  nf <- layout(matrix(c(lfor), byrow=TRUE, ncol=2))
	} else if (length(melhores) == 1){
	  lfor <- c(1, 1, 2, 2)
	  nf <- layout(matrix(c(lfor), byrow=TRUE, ncol=2))
	}

	par(nf, mar = c(2.5, 4, 1.5, 1.5))
	lapply(melhores, function(X) {
		if (class(X)[1] != "try-error" | !is.null(X))
			try(plot(X, las=1))
		}
	)
	plot(P, col=mypalette[9:5], lwd=2, plot.type = "single", xlab="", ylab="", main=paste(toupper(cvMethod), "Statistic vs Cross-Validation Horizon", sep=""), las=1)
	legend("topleft", colnames(P), lwd=2, col=mypalette[9:5], cex = 1, bty="n")
}
#setMethod("plot", "cvforecast", plot.cvforecast)
