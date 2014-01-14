###########
#Libraries#
###########
require(ggplot2)
require(pwr)

#Creates array of power values given f^2-df pairs (power of .03 with 50; .03 with 51, etc)
genPwrArray.f2 <- function(explodedEff,explodedDfs,dfNumerator,sig.level=.05){
	if(is.null(dfNumerator)){
		stop("Must specify the numerator for the F-test (k-1 groups, or k predictors)")
	}
	pwr <- pwr.f2.test(u=dfNumerator,v=explodedDfs,f2=explodedEff,sig.level=sig.level,power=NULL)$power
	return(pwr)
}
#Creates array of power values given r-df values
genPwrArray.r <- function(explodedEff,explodedDfs,sig.level=.05,alternative="two.sided"){
	pwr <- pwr.r.test(n=explodedDfs,r=explodedEff,sig.level=sig.level,alternative=alternative)$power
	return(pwr)
}
genPwrArray.r2 <- function(explodedEff,explodedDfs,sig.level=.05,alternative="two.sided"){
	explodedEff <- sqrt(explodedEff)
	pwr <- pwr.r.test(n=explodedDfs,r=explodedEff,sig.level=sig.level,alternative=alternative)$power
	return(pwr)
}
#Creates array of power values given d-df values
genPwrArray.d <- function(explodedEff,explodedDfs,sig.level=.05,alternative="two.sided",type="two.sample"){
	pwr <- pwr.t.test(n=explodedDfs,d=explodedEff,sig.level=sig.level,alternative=alternative)$power
	return(pwr)
}


#Returns list of 'exploded' effectSize array and df array.
#Namely, given an array of possible ef's and df's, it will duplicate ef's for each possible df,
#and repeat the array of df options for each ef.
#Ex: c(.03,.05),c(10,100,1000) -> .03,.03,.03,.05,.05,.05 and 10,100,1000,10,100,1000
#Used to create columns in power data.frame below.
explodeEffDfs <- function(effectSizes,dfs){
	effectEx <- rep(effectSizes,each=length(dfs))
	dfsEx <- rep(dfs,length(effectSizes))
	exploded <- list(effect=effectEx,df=dfsEx)
	return(exploded)
}

#Takes unexploded effect sizes and df's, generates a data.frame
#data.frame has |effectSize|df|power
#Each effect size has all of the dfs. Think of effect sizes as factors.
#E.x., effectSize	df	power
#	.02		10	.25
#	.02		100	.80
#	.02		1000	.99
createPwrFrame <- function(effectSizes,dfs=1000,dfNumerator=NULL,sig.level=.05,test="f2",alternative="two.sided",type="two.sample",standard=FALSE){
	if(length(dfs) == 1){ ##Allow it to take one argument, and consider it max sample size.
		dfs <- seq(4,dfs,1)
	}
	if(standard == TRUE){
		es <- cohenES()
		effectSizes <- c(effectSizes,es$es[es$test == test])
	}
	exploded <- explodeEffDfs(effectSizes,dfs)
	effects <- exploded$effect #create an ef for each df
	dfs <- exploded$df #create a df for each ef
	rm(exploded) #kill the exploded list to free mem
	if(test=="f2"){
	pwr <- genPwrArray.f2(effects,dfs,dfNumerator,sig.level)
	}
	else if(test == "r2"){
	pwr <- genPwrArray.r2(effects,dfs,sig.level,alternative)
	}
	else if(test == "d"){
	pwr <- genPwrArray.d(effects,dfs,sig.level=sig.level,alternative=alternative,type=type)
	}
	else if(test == 'r'){
	pwr <- genPwrArray.r(effects,dfs,sig.level,alternative)
	}
	else{
		stop("Test must be specified. Possible values: r,f2, r2, d")
	}
	ds <- data.frame(effectSize=effects,df=dfs,power=pwr)
	attr(ds,which="test") <- test
	return(ds)	
}
#Takes a pwrFrame and cutoff value.
#Splits the frame by effect size,
#returns the first value of x that intersects with the given y-intercept
getXIntercepts <- function(pwrFrame,power){
	require(plyr)
	xs <- daply(.data=pwrFrame,.(effectSize),.fun=function(x){
		intercepts <- x[x$power >= power & x$power <= (power + .05),"df"][1]
		return(intercepts)
	})
	return(xs)
}

getYIntercepts <- function(pwrFrame,dfs){
	require(plyr)
	ys <- daply(.data=pwrFrame,.(effectSize),.fun=function(x){
		intercepts <- x[x$df %in% dfs, "power"]
		return(intercepts)
	})
	if(!is.null(dim(ys))){
		ys <- as.numeric(ys)
	}
	return(round(ys,digits=2))
}
#Functions to add power or df guides to the plot.
addDfGuide <- function(powerPlot,power,pwrFrame){
	p <- powerPlot
	p <- p + geom_hline(yintercept=power,linetype=2) #Add horizontal line
	p <- p + annotate(x=-20,y=power,geom="text",label=as.character(power),angle=45) #Add power label
	p <- p + geom_vline(xintercept = getXIntercepts(pwrFrame,power),linetype=2) #Add vertical lines
	p <- p + annotate(x=getXIntercepts(pwrFrame,power),geom="text",label=as.character(getXIntercepts(pwrFrame,power)),y=.1,angle=45)
	return(p)
}
addPowerGuide <- function(powerPlot,dfs,pwrFrame){
	p <- powerPlot
	p <- p + geom_vline(xintercept=dfs,linetype=2)
	p <- p + annotate(x=dfs,geom="text",label=as.character(dfs),y=1,angle=45)
	p <- p + geom_hline(yintercept=getYIntercepts(pwrFrame,dfs),linetype=2)
	p <- p + annotate(x=-20,y=getYIntercepts(pwrFrame,dfs),geom="text",label=as.character(getYIntercepts(pwrFrame,dfs)),angle=45)
}

#Plots the power data.frame generated from createPwrFrame()
#Will treat effectSize as factor if fewer than 8 effectsizes are given
#Otherwise, the whole thing is a scatterplot.
##v2: Removed scatterplot, because it's useless.
##: Will take either a power value and give dfs, or dfs and give power values
plotPower <- function(pwrFrame,guides=TRUE,power=NULL,df=NULL){
	pwrFrame$effectSize <- factor(pwrFrame$effectSize)
	p <- qplot(data=pwrFrame,x=df,y=power,color=effectSize,group=effectSize,geom="line")
	if(guides==TRUE){
		if(!is.null(power)){
		p <- addDfGuide(p,power,pwrFrame)
		}
		if(!is.null(df)){
		p <- addPowerGuide(p,df,pwrFrame)
		}
	}
	p <- p + theme_classic()
	p <- p + labs(x="Degrees of freedom",y='Power',title='Power Analysis',colour=paste0("Effect Size ",'(',attr(pwrFrame,which='test'),')'))
	return(p)
}

#Converts r2 to f2; f2 is the ratio of [relatively] explained variance to unexplained variance
r2tof2 <- function(r2,r2baseModel=0){
	f2 <- (r2 - r2baseModel)/(1-r2)
	return(f2)
}

cohenES <- function(){
	es <- data.frame(test=rep(c("f2","r2","r","d"),each=3),es=c(.02,.15,.35,.01,.09,.25,.1,.3,.5,.2,.5,.8))
	return(es)
}