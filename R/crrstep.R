crrstep <-
function(formula, scope.min= ~1, etype, ... , subset, data, direction = c("backward", "forward"), criterion = c("AIC", "BICcr", "BIC"), crr.object=FALSE, trace = TRUE, steps = 100) {
#
#################################################
# Input variables:
#
# formula = formula object where LHS is failure time and RHS is linear predictors; intercept `1' should always be included
# For example:  formula <- y ~ 1 + x1 + as.factor(x2)
# scope.min = formula object denoting final model for backward selection and starting model for forward selection
# etype = integer variable that denotes type of failure for each person 
# ... = variables passed to `crr' function; two key variables are failcode and cencode (see following)
# failcode = integer value that denotes primary failure
# cencode = integer denoting censoring
# other `etype' values denote competing events
# subset = subset of data to be used for model selection
# data = data-frame containing all the variables
# criteria = selection criteria; default is AIC; BICcr uses log(n*) as the penalty where n* is the total number of Type I events
# crr.object = logical variable indicating whether to return final `crr' object
# trace = logical indicating whether to display stepwise model selection process
# steps = maximum number of steps in stepwise selection
#################################################
# All the required functions are defined first
#
    crr.extractAIC <- function (object, cov1, k = 2){
	round(-2*object$loglik + k*ncol(cov1), 2) 
    }

    crr.makecov1 <- function(formula, data){
      cov1 <- model.matrix(formula, data = as.data.frame(data))
      return(cov1)
    }

    crr.dropterm <- function (scope.max, scope.min, ftime, fstatus, data) { 
      Terms <- terms(scope.max)
      scope <- attr(Terms, "term.labels")
      if (missing(scope.min)) scope.min <- as.formula(~1)
      else if (is.null(scope.min)) scope.min <- as.formula(~1)
      Terms1 <- terms(scope.min)
      scope1 <- attr(Terms1, "term.labels")
      cov1 <- crr.makecov1(formula=scope.max, data=data)
      cov11 <- crr.makecov1(formula=scope.min, data=data)
      covtmp <- cov1
      cov1 <- cov1[,-1, drop = FALSE]
      attr(cov1, "assign") <- attr(covtmp, "assign")[-1]
      rm(covtmp)
      covtmp1 <- cov11
      cov11 <- cov11[,-1, drop = FALSE]
      attr(cov11, "assign") <- attr(covtmp1, "assign")[-1]
      rm(covtmp1)
      ns <- length(scope)
      ns1 <- length(scope1)
      if(ns1 == 0) {
      	ans <- matrix(nrow = ns+1, ncol = 1, dimnames = list(c("<none>", scope), "AIC"))	
      	cnames2 <- attributes(cov1)$assign
        cov1.2 <- cov1
        init <- crr(ftime, fstatus, cov1.2, ...)
        ans[1,] <- crr.extractAIC(init, cov1.2,k)
        
        if (ncol(cov1.2) == 0){
       		ans[1,] <- -2*init$loglik.null
       	} 
        else if (ncol(cov1.2) == 1){
          for (i in seq(ns)){
            remove.x <- (1:length(cnames2))[cnames2 == i]
            cov1.3 <- cov1.2[, -remove.x, drop=FALSE]
            ans[i+1,] <- round(-2*init$loglik.null, 2)
          }
        } else { 
          for (i in seq(ns)){
            remove.x <- (1:length(cnames2))[cnames2 == i]
            cov1.3 <- cov1.2[, -remove.x, drop=FALSE]

            if (ncol(cov1.3) == 0) {
              ans[i+1,] <- round(-2*init$loglik.null, 2)
              object2 <- crr(ftime, fstatus, cov1.2, ...)
            } else {
              object2 <- crr(ftime, fstatus, cov1.3, ...)
              ans[i+1,] <- crr.extractAIC(object2, cov1.3,k)
            }
          }
        }
     } else { 
       scope.new2 <- scope[!(scope %in% scope1)]
       ns.new <- length(scope.new2)
       ans <- matrix(nrow = ns.new+1, ncol = 1, dimnames = list(c("<none>", scope.new2), "AIC"))
       init <- crr(ftime, fstatus, cov1, ...)

       ans[1,] <- crr.extractAIC(init, cov1,k)
       cov1.2 <- cov1
       for(i in seq(length(scope.new2))){
         remove.x <- scope.new2[i]
         remove.x0 <- grep("\\(*\\)", remove.x)
         if (length(remove.x0) == 0) remove.x <- remove.x
         else if(length(remove.x0) > 0) {
           remove.x1 <- gsub("\\(", "\\\\\\(", remove.x)
           remove.x1 <- gsub("\\)", "\\\\\\)", remove.x1)
           remove.x <- remove.x1
         }
         cov1.3 <- cov1.2[, -c(grep(remove.x, colnames(cov1.2))), drop=FALSE]
         object2 <- crr(ftime, fstatus, cov1.3, ...)
         ans[i+1,] <- crr.extractAIC(object2, cov1.3,k)
       }
     }
     aod <- data.frame(AIC = ans[,1])
     head <- c("Single term deletions", "\nModel:", deparse(as.vector(init$call)))
     class(aod) <- "data.frame"
     attr(aod, "heading") <- head
     aod
 }

    crr.addterm <- function(scope.max, scope.min, ftime, fstatus, data){
      Terms <- terms(scope.max)
      scope <- attr(Terms, "term.labels")
      if (missing(scope.min)) scope.min <- as.formula(~1)
      else if (is.null(scope.min)) scope.min <- as.formula(~1)
      Terms1 <- terms(scope.min)
      scope1 <- attr(Terms1, "term.labels")
      cov1 <- crr.makecov1(scope.max, data)
      cov11 <- crr.makecov1(scope.min, data)
      covtmp <- cov1
      cov1 <- cov1[,-1, drop = FALSE]
      attr(cov1, "assign") <- attr(covtmp, "assign")[-1]
      rm(covtmp)
      covtmp1 <- cov11
      cov11 <- cov11[,-1, drop = FALSE]
      attr(cov11, "assign") <- attr(covtmp1, "assign")[-1]
      rm(covtmp1)
      ns <- length(scope)
      ns1 <- length(scope1)
 
      if(ns1 == 0) { # if scope.min == ~1
      	ans <- matrix(nrow = ns+1, ncol = 1, dimnames = list(c("<none>", scope), "AIC"))
      	cnames2 <- attributes(cov1)$assign
        cov1.2 <- cov1 
        init <- crr(ftime, fstatus, cov1, ...)
        ans[1,] <- round(-2*init$loglik.null, 2)
        for (i in seq(ns)){
          add.x <- (1:length(cnames2))[cnames2 == i]
          cov1.3 <- cov1.2[, c(add.x), drop=FALSE]
          object2 <- crr(ftime, fstatus, cov1.3, ...)
          ans[i+1,] <- crr.extractAIC(object2, cov1.3,k)
        }
      } else { 
        scope.new2 <- scope[!(scope %in% scope1)]
        ns.new <- length(scope.new2)
        ans <- matrix(nrow = ns.new+1, ncol = 1, dimnames = list(c("<none>", scope.new2), "AIC"))	
      	init <- crr(ftime, fstatus, cov11, ...)
        ans[1,] <- crr.extractAIC(init, cov11,k)
        cov1.2 <- cov1
        for(i in seq(length(scope.new2))){
          add.x <- scope.new2[i]
          add.x0 <- grep("\\(*\\)", add.x)
          if (length(add.x0) == 0) add.x <- add.x
          else if(length(add.x0) > 0) {
             add.x1 <- gsub("\\(", "\\\\\\(", add.x)
             add.x1 <- gsub("\\)", "\\\\\\)", add.x1)
             add.x <- add.x1
          }
          cov1.3 <- cbind(cov11, cov1.2[,c(grep(add.x, colnames(cov1.2)))])
          colnames(cov1.3) <- c(colnames(cov11), colnames(cov1.2)[c(grep(add.x, colnames(cov1.2)))])
          object2 <- crr(ftime, fstatus, cov1.3, ...)
          ans[i+1,] <- crr.extractAIC(object2, cov1.3,k)
        }
      }
     aod <- data.frame(AIC = ans[,1])
     head <- c("Single term additions", "\nModel:", deparse(as.vector(init$call)))
     class(aod) <- "data.frame"
     attr(aod, "heading") <- head
     aod
   }
 
  cut.string <- function(string) {
	if(length(string) > 1)
		string[-1] <- paste("\n", string[-1], sep= "")
		string	
	}
###########################################################
### Main body starts here ###

stopifnot(identical(formula[[1]], as.name("~")))
lhs <- if(length(formula)==3) formula[[2]] else NULL
scope.max <- as.formula(formula[-2])

stopifnot(identical(scope.min[[1]], as.name("~")))
scope.min <- if(length(scope.min)==3) as.formula(scope.min[-2]) else as.formula(scope.min)

data <- data[subset, ]

Xmat <- model.frame(formula, data, na.action=na.pass)
nomiss <- apply(Xmat, 1, function(x) !any(is.na(x)))
data <- data[nomiss, ]
ftime <- Xmat[nomiss ,1]
fstatus <- eval(substitute(etype), data)
fstatus <- fstatus[nomiss]
criterion <- match.arg(criterion, choices=c("AIC", "BICcr", "BIC"))
if (criterion == "AIC") k <- 2
if (criterion == "BICcr") k <- log(sum(fstatus == 1))
if (criterion == "BIC") k <- log(nrow(data))
rm(Xmat)
invisible(gc())

	call <- match.call()
	if (trace) print(call)
    
   direction <- match.arg(direction)
   backward <- direction == "backward"
   forward <- direction == "forward"
   formula <- scope.max
   cov1 <- crr.makecov1(formula=formula, data=data)[, -1, drop=FALSE]
   object <- crr(ftime, fstatus, cov1, ...)
   fit <- object
   if (backward) cov3.2 <- cov1
   if (forward) {
     cov3.1 <- vars.in <- NULL
     formula1 <- scope.min
     cov2 <- crr.makecov1(scope.min, data = data)[, -1, drop= FALSE]
   }
   while (steps > 0) {
     steps <- steps - 1	
     aod <- NULL
     change <- NULL
     if (backward) {
       if(formula==scope.min) break
       if(ncol(cov3.2) == 0) break

       aod <- crr.dropterm(scope.max=formula, scope.min=scope.min, ftime=ftime, fstatus=fstatus, data=data) 
       rn <- rownames(aod)
       row.names(aod) <- c(rn[1], paste("-", rn[-1], sep = ""))
       if(is.null(aod) || ncol(aod) == 0) break
       o <- order(aod[,1])
       aod.o <- aod[o, ,drop=FALSE]
       if (trace) {
         print(aod.o)
         utils::flush.console()
       }
       if (o[1] ==1) {
         std.error <- rep(0, ncol(cov3.2))
           for (i in 1:ncol(cov3.2)){
             std.error[i] <- sqrt(fit$var[i,i])
           }
         output <- list("variables" = colnames(cov3.2), "coefficients" = fit$coef, "std. errors" = std.error, "log.lik" = round(fit$loglik,2))
	       outmat <- cbind(signif(fit$coef,3), signif(std.error,3), round(abs(fit$coef)/std.error, 3))
        colnames(outmat) <- c("estimates", "std.error", "t-stat")
	       rownames(outmat) <- colnames(cov3.2)
         break         
       }
       change <- rownames(aod)[o[1]]
       s.split <- strsplit(change, "-")
       var <- s.split[[1]][2]
       if (trace) print(var) 
       var0 <- grep("\\(*\\)", var)
       if (length(var0) == 0) var <- var
       else if(length(var0) > 0) {
         var1 <- gsub("\\(", "\\\\\\(", var)
         var1 <- gsub("\\)", "\\\\\\)", var1)
         var <- var1
       }
       ans <- sub(paste("\\+ ", var, sep=""), "", formula)
       formula <- as.formula(paste(ans[1], ans[2], sep=""))
       if (!is.na(var)){
         cov3.2 <- cov3.2[,-c(grep(var, colnames(cov3.2))), drop=FALSE]
         if (ncol(cov3.2) != 0) {
           fit <- update(object, cov1 = cov3.2)
           std.error <- rep(0, ncol(cov3.2))
           for (i in 1:ncol(cov3.2)){
             std.error[i] <- sqrt(fit$var[i,i])
           }
           output <- list("variables" = colnames(cov3.2), "coefficients" = fit$coef, "std. errors" = std.error, "log.lik" = round(fit$loglik,2))
         } else { #if ncol(cov3.2) == 0
           fit2 <- fit
           fit <- NULL
           null.output <- paste("The best model is the NULL model and the likelihood is:", round(fit2$loglik.null,2), sep="")
           print(null.output)
           output <- list("variables" = NULL, "coefficients" = NA, "std. errors" = NA, "log.lik" = round(fit2$loglik.null,2))
	outmat <- NULL
	} ########################## Here is where I think the problem is #########
       } else if (is.na(var)) { # if is.na(var); deleting no variables gives lowest AIC
       	 scope.max <- as.formula(paste(ans[2], ans[1], ans[3], sep="")) #revised this from line below
         #scope.max <- as.formula(paste(ans[1], ans[2], sep="")) #this would only give ~Y, which is not the correct formula
         cov3.2 <- crr.makecov1(scope.max, data)[,-1, drop = FALSE]
         if (ncol(cov3.2) != 0){ #cov1 contains >= 1 variable 
           fit <- update(object, cov1 = cov3.2)	
           std.error <- rep(0, ncol(cov3.2))
           for (i in 1:ncol(cov3.2)){
             std.error[i] <- sqrt(fit$var[i,i])
           }
           output <- list("variables" = colnames(cov3.2), "coefficients" = fit$coef, "std. errors" = std.error, "log.lik" = round(fit$loglik,2))
         } else{ # if no variables in model
           fit2 <- fit
           fit <- NULL
           null.output <- paste("The best model is the NULL model and the likelihood is:", round(fit2$loglik.null,2), sep="")
           print(null.output)
           output <- list("variables" = NULL, "coefficients" = NA, "std. errors" = NA, "log.lik" = round(fit2$loglik.null,2))
         }
       }
 if (!is.null(fit)) {
	outmat <- cbind(signif(fit$coef,2), signif(std.error,2), round(abs(fit$coef)/std.error, 2))
	  colnames(outmat) <- c("estimates", "std.error", "t-stat")
	  rownames(outmat) <- colnames(cov3.2)
	} else fit <- fit2

       if (trace){
         cat("\nStep: AIC = ", format(round(crr.extractAIC(fit, cov3.2, k),2)), "\n", cut.string(deparse(as.vector(fit$call))), "\n\n", sep="")
         utils::flush.console()
       }
### End of Backward selection ###
     }
##### Forward selection ######
	if (forward) {
		aod <- crr.addterm(formula, formula1, ftime, fstatus, data)
      	var <- NULL
      	rn <- rownames(aod)
      	if (trace) print(vars.in)
      	invar <- if (is.null(vars.in)) rep(FALSE, length(rn)) else (rn %in% vars.in)
      	aod <- aod[!invar, , drop=FALSE]
      	rn <- rownames(aod)
      	row.names(aod) <- c(rn[1], paste("+", rn[-1], sep = ""))
      	if(is.null(aod) || ncol(aod) == 0) break
      	o <- order(aod[,1])
      	aod.o <- aod[o, ,drop=FALSE]
      	if (trace) {
       	 print(aod.o)
       	 utils::flush.console()
      	}
      
      	if (o[1] ==1) {                                         
      		formula1 <- as.formula(paste(scope.min[1], scope.min[2], sep=""))
        	if (formula1 != ~1 ) {
          	cov3.1 <- crr.makecov1(formula1, data)[,-1, drop = FALSE]
          	fit <- update(object, cov1 = cov3.1)
          	std.error <- rep(0, ncol(cov3.1))
          	for (i in 1:ncol(cov3.1)){
            	std.error[i] <- sqrt(fit$var[i,i])
          	}
          	output <- list("variables" = colnames(cov3.1), "coefficients" = fit$coef, "std.errors" = std.error, "log.lik" = round(fit$loglik, 2))
        	} else { # formula == ~1
          	null.output <- paste("The best model is the NULL model and the likelihood is:", round(fit$loglik.null,2), sep="")
          	if (trace) print(null.output)
          	fit2 <- fit
          	fit <- NULL
          	output <- list("variables" = NULL, "coefficients" = NA, "std. errors" = NA, "log.lik" = round(fit2$loglik.null,2))
		outmat <- NULL
	  	fit$loglik <- fit2$loglik.null

        	}
      		break
      	}
      	change <- rownames(aod)[o[1]]
      	s.split <- strsplit(change, "\\+")
      	var <- s.split[[1]][2]
      	print(var)
      	var0 <- grep("\\(*\\)", var)
      	if (length(var0) == 0) var1 <- var 
      	else if(length(var0) > 0) {
       	 	var1 <- gsub("\\(", "\\\\\\(", var)
        	var1 <- gsub("\\)", "\\\\\\)", var1)
      	}
      	ans <- paste("+", var, sep = "")
        formula1 <- as.formula(paste(scope.min[1], scope.min[2], ans, sep=""))
        scope.min <- formula1
        cov2 <- cbind(cov2, cov1[, c(grep(var1, colnames(cov1))), drop=FALSE])
        fit <- update(object, cov1 = cov2)
        vars.in <- c(vars.in, var)
        std.error <- rep(0, ncol(cov2))
        for (i in 1:ncol(cov2)){
        	std.error[i] <- sqrt(fit$var[i,i])
        }
           output <- list("variables" = NULL, "coefficients" = NA, "std. errors" = NA, "log.lik" = round(fit$loglik.null,2))

	if (!is.null(fit)) {
	outmat <- cbind(signif(fit$coef,2), signif(std.error,2), round(abs(fit$coef)/std.error, 3))
	  colnames(outmat) <- c("estimates", "std.error", "t-stat")
	  rownames(outmat) <- colnames(cov2)
	} else fit <- fit2

        Terms <- terms(formula)
        scope <- attr(Terms, "term.labels")
        Terms1 <- terms(formula1)
        scope1 <- attr(Terms1, "term.labels")
        scope <- sort(scope)
        scope1 <- sort(scope1)
        scope.new2 <- scope[!(scope %in% scope1)]
        if(length(scope.new2) == 0) break
              
       	if (trace){
         	cat("\nStep: AIC = ", format(round(crr.extractAIC(fit, cov2, k),2)), "\n", cut.string(deparse(as.vector(fit$call))), "\n\n", sep="")
         	utils::flush.console()
       	} 
##### End of Forward selection #####
    	} 
### End of while loop ###
   	}
	if (trace) print(output)
   	if (crr.object == TRUE) return(fit)
    return(list(coefficients=outmat, log.likelihood=round(fit$loglik, 2)))
}

