#' Generic function to perform multiple forecasts by CrossValidation
#'
#' This function receives data and performes multiple forecasts by the technique of CrossValidation. The decision about the best model is based on tests as linearity, trend and fit accuracy.
#'
#'
#' @param tsdata data.frame type date-value, ts, mts or xts time series objects
#' @param tsControl generic contol with several args for the modelling process. See
#'    \code{\link{cvForecastControl}}.
#' @param ... other arguments
#' @return A list of class 'cvforecast' containing several objetcts from the forecasting process. It includes: betso models (less tahn 6), crossValidation statistics for all models, accuracy of all models, the control, etc. As below.
#' @keywords crossValidation time series
#' @examples
#'#Define cross validation parameters
#'myControl <- cvForecastControl(
#'  minObs = 14,
#'  stepSize = 5,
#'  maxHorizon = 30,
#'  fixedWindow=TRUE,
#'  preProcess=FALSE,
#'  ppMethod='guerrero',
#'  summaryFunc=tsSummary,
#'  cvMethod="MAPE",
#'  tsfrequency='day',
#'  OutlierClean=FALSE,
#'  dateformat='%d/%m/%Y %H:%M:%S')
#'
#'#Paralell execution improves the processing time
#'
#' #cl <- makeCluster(4, type='SOCK')
#' #registerDoParallel(cl)
#'
#'#Load data
#'
#'data(diario, package="cvforecast")
#'dadosd <- ConvertData(diario, dateformat='%d/%m/%Y %H:%M:%S', tsfrequency = "day", OutType="ts")
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
cvforecast <- function(tsdata, tsControl=cvForecastControl(), ...) {

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

	nmforecast <- Try_error(forecastMethod(x))

	if(length(as.numeric(x)) < 2*max(cycle(x)) | max(cycle(x))==1 | class(nmforecast) == "try-error") {
		cat("Series with less than 2 cycles or non-periodic. Try ARIMA and ETS!\n")
		tsControl$minObs <- 8
		tsControl$stepSize <- 1
		nmforecast <- list("auto.arimaForecast","etsForecast")
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

	## monta tabela com as estat
	STATS <- lapply(out,
		function(X) {
			if (class(X)[1]!="try-error") X$results
			else NA
		}
	)

	STATS <- Filter(Negate(function(X) is.null(unlist(X))), STATS)

	## Extra??o da estta?stica de CrossValida??o (MAPE)
	STATS1 <- lapply(STATS,
		function(X) {
			if (is.null(dim(X)) | mean(as.numeric(X[, cvMethod]), na.rm=TRUE)%in% c(NA, 0, -Inf, Inf)) {
				tsControl$cvMethod <- cvMethod <- "MAE"
				return(as.numeric(X[, cvMethod]))
			} else { # Utiliza a estatistica MAE para decisao
				return(as.numeric(X[, cvMethod]))
			}
		}
	)

	## Montar data.frame e remover elementos missing
	estatisticas <- as.data.frame(do.call("cbind", STATS1))
	estatisticas <- estatisticas[sapply(estatisticas, function(X) length(unique(na.omit(X)))) > 1];
	estatisticas <- na.omit(estatisticas)

	## Calcular Coeficiente de Varia??o da estat?stica STATS1 no processo Cross
	## A decis?o n?o ? baseada apenas no menor valor, mas sim naqueque que aprenseta melhor comportamento ao longo do CrossVall. Ex: menor coeficiente de varia??o m?dio e coeficiante de varia??o mediano.
	CV <- sapply(estatisticas, function(X) {
			cv = sd(X, na.rm=TRUE)/mean(X, na.rm=TRUE)
			#DAMM = Desvio Absoluto Mediano da Mediana
			#damm = median(abs(X-median(X, na.rm=TRUE)), na.rm=TRUE)
			#cvmed = damm/median(X, na.rm=TRUE)
			##mpcv <- (cv+cvmed-cv*cvmed) ## Uniao dos valores AUB = A+B-AC
			return(abs(median(X, na.rm=TRUE)) + cv)
		}
	)

	## Remover estat?sticas com valores discrepantes e reordenar os modelos pelo menor
	CV <- sort(CV[!CV %in% c(Inf, -Inf, 0, NA)])
	CV_names <- as.list(names(CV))

	best_models <- plyr::llply(CV_names, function(J) {
			temp <- try(switch.cvforecast(x=x, nmodelo=J, h=h, onlyfc=FALSE))
			if(class(temp)!="try-error") temp else NA
		}
	)

	names(best_models) <- CV_names

	best_models <- best_models[!is.na(best_models)]

	Bm <- cvBestModel(best_models, cvMethod=cvMethod,residlevel = 0.10)
	CV_names <- as.list(Bm[,1])

	## Seleciona apenas os cinco melhores
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
	out1$myControl<- myControl

	class(out1) <- "cvforecast"
	return(invisible(out1))
}

