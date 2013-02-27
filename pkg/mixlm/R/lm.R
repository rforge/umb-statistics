# This function replaces stats::lm, organizes fixed and random effects, removes r() from formula and parses to lm or lmer.

# ANOVA (sequential SS)
anova.lm <- function(object,...){
	if(is.null(object$random)){
		return(stats::anova.lm(object,...))
	} else {
		stop("Only type III error Anova() supported in mixed/random models")
	}
}


Anova.lm <- function(mod,type="III", ...){
	if(!is.null(mod$random)){
		if(type=="II" || type=="2"){
			warning("Using (default) type III error supported for mixed/random models")
		}
		return(AnovaMix(mod,...))
	} else {
		a <- car:::Anova.lm(mod, type=type, ...)
		d <- cbind(a[2],a[1],a[1]/a[2],a[3],a[4]) # Add mean squares in Anova for linear model
		attr(d,"names")[3] <- "Mean Sq"
		class(d) <- class(a)
		return(d)
	}
}

## FIXME: Her antas modellen å kun inneholde faktorer, ingen andre effekter. Kan gi rare resultater!
# Mixed model ANOVA
AnovaMix <- function(object){
	formula         <- formula(object)
	formula.text    <- as.character(formula)
	all.effects     <- object$random$all							  # All model effects and interactions
	fixed.effects   <- object$random$fixed							  # All fixed effects
	random.effects  <- object$random$random						  # All random effects
	main.rands.only.inter <- object$random$main.rands.only.inter     # Random effects only present in interactions
	restrictedModel <- !object$random$unrestricted
	data    <- object$model
	n.effects    <- length(all.effects)
	main.effects <- fparse(formula)							  # All main effects (even though only included in interactions)
	n.levels     <- numeric(length(main.effects))
	for(i in 1:length(main.effects)){
		n.levels[i] <- length(levels(data[,main.effects[i]])) # Number of levels per main effect
	}
	names(n.levels) <- main.effects
	N <- dim(data)[1]
	
	ind.randoms <- numeric()
	ind.randoms <- match(random.effects,all.effects) # Placement of random effects in "all.effects"
	ind.fixed   <- match(fixed.effects,all.effects)  # Placement of fixed effects in "all.effects"
	ind.fixed <- setdiff(1:n.effects,ind.randoms)										
	n.randoms    <- length(ind.randoms)
					
	# Estimate fixed effect Anova
	fixed.model <- as.data.frame(car:::Anova.lm(object, type='III'))
	fixed.model <- fixed.model[-1,] # Remove intercept
	fixed.model <- cbind(fixed.model[,c(1,2)], "Mean Sq"=fixed.model[,1]/fixed.model[,2], fixed.model[,c(3,4)])
	fixed.model <- fixed.model[c(all.effects,"Residuals"),] # Sort according to all.effects

	# Check which effects should use interactions as denominators instead of error
	approved.interactions <- list()
	approved.interactions.fixed <- list()
	for(i in 1:n.effects){
		this.effect <- strsplit(all.effects[i],":")[[1]]
		which.contains <- numeric()
		for(j in 1:n.effects){ # Find all other effects containing this.effect
			effect.names <- is.element(strsplit(all.effects[j],":")[[1]],this.effect)
			if(i!=j && sum(effect.names)==length(this.effect) && length(effect.names)>length(this.effect)){
				which.contains <- union(which.contains,j)}
		}
		which.contains <- sort(which.contains)
		if(length(which.contains)>0){
			approved.interaction <- numeric(length(which.contains))
			approved.interaction.fixed <- numeric(length(which.contains))
			for(j in 1:length(which.contains)){
				if(restrictedModel){
					approved.interaction[j] <- prod(is.element(setdiff(strsplit(all.effects[which.contains],":")[[j]],strsplit(all.effects[i],":")[[1]]),c(random.effects,main.rands.only.inter)))}
				else{
					if(any(is.element(ind.fixed,i))){
						approved.interaction.fixed[j] <- prod(is.element(setdiff(strsplit(all.effects[which.contains],":")[[j]],strsplit(all.effects[i],":")[[1]]),fixed.effects))}
					approved.interaction[j] <- 1-prod(!is.element(strsplit(all.effects[which.contains],":")[[j]],c(random.effects,main.rands.only.inter)))}
			}
			if(length(which(approved.interaction==1))>0){
				approved.interactions[[i]] <- which.contains[which(approved.interaction==1)]}
			else{
				approved.interactions[[i]] <- FALSE}
			if(length(which(approved.interaction.fixed==1))>0){
				approved.interactions.fixed[[i]] <- which.contains[which(approved.interaction.fixed==1)]}
			else{
				approved.interactions.fixed[[i]] <- FALSE}
		}
		else{
			approved.interactions[[i]] <- FALSE
			approved.interactions.fixed[[i]] <- FALSE}
	}

	# Find variance components (except MSerror), 
	# and find linear combinations needed to produce denominators of F-statistics
	mix.model.attr <- list()
	denom.df <- numeric(n.effects+1)
	exp.mean.sq <- rep(paste("(",n.effects+1,")", sep=""), n.effects+1)
	var.comps <- numeric(n.effects+1)*NA
	var.comps[n.effects+1] <- fixed.model[n.effects+1,3]
	for(i in 1:n.effects) {
		if(!is.logical(approved.interactions[[i]])){
			# Set up matrix A and vector b to find linear combinations of effects to use as denominators in F statistics
			## This is probably where unbalancedness should be included !!!!!!!!
			lap <- length(approved.interactions[[i]])
			A <- matrix(0,lap+1,n.effects+1)
			b <- rep(1,lap+1)
			for(j in 1:lap){
				A[j,approved.interactions[[approved.interactions[[i]][j]]]] <- 1
				A[j,approved.interactions[[i]][j]] <- 1
				k <- length(approved.interactions[[i]])+1-j
				exp.mean.sq[i] <- paste(exp.mean.sq[i], " + ", N/prod(n.levels[strsplit(all.effects[approved.interactions[[i]][k]],":")[[1]]]), " (",which(all.effects==all.effects[approved.interactions[[i]][k]]),")", sep="")
			}
			A[, n.effects+1] <- 1
			A <- A[,apply(A,2,sum)>0]
			denominator <- solve(t(A),b)
			denominator.id <- c(approved.interactions[[i]],n.effects+1)
			denominator.id <- denominator.id[denominator!=0]
			mix.model.attr[[i]] <- denominator <- denominator[denominator!=0]
			names(mix.model.attr[[i]]) <- denominator.id
			if(length(denominator)==1){ # Original df
				denom.df[i] <- fixed.model[denominator.id,2]}
			else{ # Satterthwaite's df correction
				denom.df[i] <- sum(fixed.model[denominator.id,3]*denominator)^2/sum((fixed.model[denominator.id,3]*denominator)^2/fixed.model[denominator.id,2])} 
		} else{
			denominator.id <- n.effects+1
			mix.model.attr[[i]] <- 1
			names(mix.model.attr[[i]]) <- denominator.id
			denom.df[i] <- fixed.model[denominator.id,2]
			denominator <- 1
		}
		if(sum(ind.randoms==i)>0){
			exp.mean.sq[i] <- paste(exp.mean.sq[i], " + ", N/prod(n.levels[strsplit(all.effects[i],":")[[1]]]), " (",i,")", sep="")
			var.comps[i] <- (fixed.model[i,3]-fixed.model[denominator.id,3]%*%denominator)/(N/prod(n.levels[strsplit(all.effects[i],":")[[1]]]))}
		else{
			if(!is.logical(approved.interactions.fixed[[i]])){
				ex.ind <- paste(",", paste(approved.interactions.fixed[[i]], sep="", collapse=","),sep="")}
			else{
				ex.ind <- ""}
			exp.mean.sq[i] <- paste(exp.mean.sq[i], " + ", N/prod(n.levels[strsplit(all.effects[i],":")[[1]]]), " Q[",i,ex.ind,"]", sep="")
		}
		fixed.model[i,4] <- fixed.model[i,3]/(fixed.model[denominator.id,3]%*%denominator)
		if(fixed.model[i,4]<0){
			fixed.model[i,4] <- NA}
		fixed.model[i,5] <- 1-pf(fixed.model[i,4],fixed.model[i,2],denom.df[i])
	}
	names(denom.df) <- rownames(fixed.model)
	object <- list(lm=object, anova=fixed.model, err.terms=c(mix.model.attr,NA), denom.df=denom.df, restricted=restrictedModel,
		exp.mean.sq=exp.mean.sq, var.comps=var.comps, random.effects=random.effects, ind.randoms=ind.randoms, formula.text=formula.text)
	class(object) <- "AnovaMix"
	object
}


## Print method for object from AnovaMix
print.AnovaMix <- function(x,...){
	object <- x
	N <- length(object$err.terms)
	output1 <- object$anova
	Fs <- PrF <- character(N)
	PrF[!is.na(output1$"Pr(>F)")] <- format(round(output1$"Pr(>F)"[!is.na(output1$"Pr(>F)")],4), digits=1, scientific=FALSE, nsmall=4)
	PrF[is.na(output1$"Pr(>F)")] <- "-"
	output1$"Pr(>F)" <- PrF
	Fs[!is.na(output1$"F value")] <- format(output1$"F value"[!is.na(output1$"F value")], digits=1, scientific=FALSE, nsmall=2)
	Fs[is.na(output1$"F value")] <- "-"
	output1$"F value" <- Fs
	output1$"Sum Sq" <- format(output1$"Sum Sq", digits=1, scientific=FALSE, nsmall=2)
	output1$"Mean Sq" <- format(output1$"Mean Sq", digits=1, scientific=FALSE, nsmall=2)

	err.terms <- character(length(object$err.terms))
	for(i in 1:N){
		if(length(object$err.terms[[i]])==1 && is.na(object$err.terms[[i]])){
			err.terms[i] <- "-"
		}
		else{
			err.terms[i] <- paste(ifelse(object$err.terms[[i]][1]>1,paste(object$err.terms[[i]][1],"*",sep=""),""),"(",names(object$err.terms[[i]][1]),")",sep="")
			if(length(object$err.terms[[i]])>1){
				for(j in 2:length(object$err.terms[[i]])){
					if(object$err.terms[[i]][j]<0){
						err.terms[i] <- paste(err.terms[i], paste(ifelse(object$err.terms[[i]][j]<(-1),paste(abs(object$err.terms[[i]][j]),"*",sep=""),""), "(", names(object$err.terms[[i]][j]), ")", sep=""), sep=" - ")
					} else {
						err.terms[i] <- paste(err.terms[i], paste(ifelse(object$err.terms[[i]][j]>1,paste(object$err.terms[[i]][j],"*",sep=""),""), "(", names(object$err.terms[[i]][j]), ")", sep=""), sep=" + ")}
				}
			}
		}
	}
	var.comps <- format(object$var.comps, digits=3)
	var.comps[setdiff(1:(N-1), object$ind.randoms)] <- "fixed"

	denom.df <- character(N)
	denom.df[!is.na(object$denom.df)&round(object$denom.df)==object$denom.df] <- format(object$denom.df[!is.na(object$denom.df)&round(object$denom.df)==object$denom.df], digits=3)
	denom.df[!is.na(object$denom.df)&round(object$denom.df)!=object$denom.df] <- format(object$denom.df[!is.na(object$denom.df)&round(object$denom.df)!=object$denom.df], digits=3)
	denom.df[object$denom.df==0] <- "-"
	output2 <- data.frame("Err.terms"=err.terms, "Denom.df"=denom.df, "VC(SS)"=var.comps)
	colnames(output2) <- c("Err.term(s)", "Err.df", "VC(SS)")
	output3 <- data.frame("E(MS)"=format(object$exp.mean.sq))
	colnames(output3) <- "Expected mean squares"
	rownames(output2) <- paste(1:N," ",rownames(object$anova), sep="")
	rownames(output3) <- rownames(object$anova)
	if(!object$restricted){
		un <- "un"}
	else{
		un <- ""}
	cat("Analysis of variance (", un, "restricted model)\n", sep="")
	cat("Response: ", object$formula.text[2], "\n", sep="")
	print(format(output1, digits=3))
	cat("\n")
	print(output2)
	cat("(VC = variance component)\n\n")
	print(output3)
	if(!is.balanced(object$lm)){
		cat("\nWARNING: Unbalanced data may lead to poor estimates\n")
	}
}

# lm from stats, edited to use treatment names in sum contrasts
# and enable limited classical least squares mixed models
lm <- function (formula, data, subset, weights, na.action,
		method = "qr", model = TRUE, x = FALSE, y = FALSE,
		qr = TRUE, singular.ok = TRUE, contrasts = NULL,
		offset, unrestricted = TRUE, REML = NULL, ...)
{
    ret.x <- x
    ret.y <- y
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
	## Edited by KHL
	mfd <- match(c("formula","data"), names(mf), 0L)
	if(length(mfd)==2){ # Has formula and data
		is.random <- TRUE
		if( any(grepl("r(",formula,fixed=TRUE)) ){
			rw <- random.worker(formula, data, REML)
		} else {
			rw <- list(0)
		}
		if(length(rw) == 1){
			is.random <- FALSE
		} else { # Removed r() from formula
			formula <- rw$formula
			mf$formula <- rw$formula
			rw$unrestricted <- unrestricted
			if(is.logical(REML)){ # Perform 
				if(!require(lme4))
					stop("Package lme4 required")
				cl[[1]] <- as.name("lmer")
				cl[["formula"]] <- rw$reml.formula
				object <- eval(cl,parent.frame())
				# object <- lme4::lmer(rw$reml.formula, data, REML = REML, contrasts = contrasts, subset = subset, weights = weights, na.action = na.action, model = model, ...)
				object@call <- cl
				return(object)
			}
		}
	} else {
		is.random <- FALSE
	}
	## End of edit
    m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"),
	       names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    if (method == "model.frame")
	return(mf)
    else if (method != "qr")
	warning(gettextf("method = '%s' is not supported. Using 'qr'", method),
                domain = NA)
    mt <- attr(mf, "terms") # allow model.frame to update it
    y <- model.response(mf, "numeric")
    ## avoid any problems with 1D or nx1 arrays by as.vector.
    w <- as.vector(model.weights(mf))
    if(!is.null(w) && !is.numeric(w))
        stop("'weights' must be a numeric vector")
    offset <- as.vector(model.offset(mf))
    if(!is.null(offset)) {
        if(length(offset) != NROW(y))
            stop(gettextf("number of offsets is %d, should equal %d (number of observations)",
                          length(offset), NROW(y)), domain = NA)
    }

    if (is.empty.model(mt)) {
	x <- NULL
	z <- list(coefficients = if (is.matrix(y))
                  matrix(,0,3) else numeric(), residuals = y,
		  fitted.values = 0 * y, weights = w, rank = 0L,
		  df.residual = if(!is.null(w)) sum(w != 0) else
                  if (is.matrix(y)) nrow(y) else length(y))
        if(!is.null(offset)) {
            z$fitted.values <- offset
            z$residuals <- y - offset
        }
    }
    else {
	x <- model.matrix(mt, mf, contrasts)
	## Edited by KHL
	if(is.null(contrasts) && (options("contrasts")[[1]][1]!="contr.treatment" || options("contrasts")[[1]][1]!="contr.poly") && !missing(data)){
		col.names   <- effect.labels(mt,data) # mt er "terms" fra formula, x er model.matrix
		if(length(col.names)==length(colnames(x))){
			colnames(x) <- effect.labels(mt,data)
		}
	}
	## End edit
	z <- if(is.null(w)) lm.fit(x, y, offset = offset,
                                   singular.ok=singular.ok, ...)
	else lm.wfit(x, y, w, offset = offset, singular.ok=singular.ok, ...)
    }
    class(z) <- c(if(is.matrix(y)) "mlm", "lm")
    z$na.action <- attr(mf, "na.action")
    z$offset <- offset
    z$contrasts <- attr(x, "contrasts")
    z$xlevels <- .getXlevels(mt, mf)
    z$call <- cl
    z$terms <- mt
    if (model)
	z$model <- mf
    if (ret.x)
	z$x <- x
    if (ret.y)
	z$y <- y
    if (!qr) z$qr <- NULL
	## Edited by KHL
	if( is.random ){
		z$random <- rw
		if(!all(grepl("factor",attr(mt,"dataClasses")[-1]))){
			stop("Mixed models containing continuous effects not supported")
		}
	}
	## End edit
    z
}

## Collect and extract randomness
random.worker <- function(formula, data, REML = NULL){
	formula <- formula(formula)
	terms <- terms(formula)
	effsr <- attr(terms,"term.labels")
	effs  <- attr(terms(rparse(formula)),"term.labels")
	if(length(effs)==0){
		return( list(0) )
	}

	has.intercept <- attr(terms,"intercept")==1
	rands <- sort(unique(c(grep("[:]r[(]",effsr),   # Match random in interaction
						   grep("^r[(]",  effsr),   # Match random in the beginning
						   grep("[(]r[(]",effsr)))) # Match random inside function

	# which.rands <- match(rands,effsr)
    eff.splits <- list()
	for(i in 1:length(effs)){ # Split effect to look for hidden random interactions
		eff.splits[[i]] <- fparse(formula(paste("1~", effs[i],sep="")))
	}
	eff.lengths <- lapply(eff.splits,length)
	main.effs   <- effs[eff.lengths==1]
	main.rands  <- main.effs[main.effs%in%effs[rands]]
	main.rands.only.inter <- character(0)
	for(i in rands){
		main.rands.only.inter <- c(main.rands.only.inter, setdiff(eff.splits[[i]],main.effs)) # Random main effects only present in interactions
	}
	inter.rands <- which(unlist(lapply(eff.splits,function(i) any(main.rands%in%i))))
	# Check if any interactions containing random effects are not labeled as random
	if(any(is.na(match(inter.rands,rands)))){
		extra.randoms <- inter.rands[which(is.na(match(inter.rands,rands)))]
		warning(paste(paste(effs[extra.randoms],sep="",collapse=", "), " included as random interaction",ifelse(length(extra.randoms)==1,"","s"),sep=""))
		rands <- cbind(rands,extra.randoms)
		effs  <- effs[!(extra.randoms%in%effs)]
	}
	if(length(rands)==0){
		return( list(0) ) 
	} else {
		if(is.logical(REML)){
			remleffs     <- c(effs[setdiff(1:length(effs),rands)],paste("(1|",effs[rands],")",sep=""))
			reml.formula <- formula(paste(formula[[2]],"~",paste(remleffs,collapse="+"),ifelse(has.intercept,"","-1"),sep=""))

			return( list(formula = formula(paste(formula[[2]],"~",paste(effs,collapse="+"),ifelse(has.intercept,"","-1"),sep="")), random = effs[rands], main.rands.only.inter = main.rands.only.inter, fixed = effs[setdiff(1:length(effsr),rands)], all = effs, has.intercept = has.intercept, remleffs = remleffs, reml.formula = reml.formula))
		} else {
			return( list(formula = formula(paste(formula[[2]],"~",paste(effs,collapse="+"),ifelse(has.intercept,"","-1"),sep="")), random = effs[rands], main.rands.only.inter = main.rands.only.inter, fixed = effs[setdiff(1:length(effsr),rands)], all = effs, has.intercept = has.intercept))
		}
	}
}
