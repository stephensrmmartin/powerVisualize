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
createPwrFrame <- function(effectSizes,dfs,dfNumerator=NULL,sig.level=.05,test="f2"){
	if(length(dfs) == 1){ ##Allow it to take one argument, and consider it max sample size.
		dfs <- seq(1,dfs,1)
	}
	exploded <- explodeEffDfs(effectSizes,dfs)
	effects <- exploded$effect #create an ef for each df
	dfs <- exploded$df #create a df for each ef
	rm(exploded) #kill the exploded list to free mem
	if(test=="f2"){
	pwr <- genPwrArray.f2(effects,dfs,dfNumerator,sig.level)
	}
	else if(test == "r2"){
	pwr <- genPwrArray.r2(effects,dfs,sig.level)
	}
	else if(test == "d"){
	pwr <- genPwrArray.d(effects,dfs,sig.level)
	}
	else{
		stop("Test must be specified. Possible values: f2, r2, d")
	}
	ds <- data.frame(effectSize=effects,df=dfs,power=pwr)
	return(ds)	
}
#Takes a pwrFrame and cutoff value.
#Splits the frame by effect size,
#returns the first value of x that intersects with the given y-intercept
getXIntercepts <- function(pwrFrame,cutoff){
	require(plyr)
	daply(.data=pwrFrame,.(effectSize),.fun=function(x){
		intercepts <- x[x$power >= cutoff & x$power <= (cutoff + .05),"df"][1]
		return(intercepts)
	})
}

getYIntercepts <- function(pwrFrame,dfs){
	require(plyr)
	daply(.data=pwrFrame,.(effectSize),.fun=function(x){
		intercepts <- x[x$df %in% dfs, "power"]
	})
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

#Plots the power data.frame generated from createPwrFrame()
#Will treat effectSize as factor if fewer than 8 effectsizes are given
#Otherwise, the whole thing is a scatterplot.
##v2: Removed scatterplot, because it's useless.
##: Will take either a power value and give dfs, or dfs and give power values
plotPower <- function(pwrFrame,guides=TRUE,power=NULL,df=NULL){
	if(length(unique(pwrFrame$effectSize)) < 8){
		pwrFrame$effectSize <- factor(pwrFrame$effectSize)
		p <- qplot(data=pwrFrame,x=df,y=power,color=effectSize,group=effectSize,geom="line")
		if(guides==TRUE & !is.null(power)){
			p <- addDfGuide(p,power,pwrFrame)
		}
	}
	else{
		p <- qplot(data=pwrFrame,x=df,y=power,color=effectSize,geom="point") + scale_color_gradient2(low='red',mid='blue',high='black')
	}
	p <- p + theme_classic()
	p <- p + labs(x="Degrees of freedom",y='Power',title='Power Analysis',colour='Effect Size')
	return(p)
}

#Converts r2 to f2; f2 is the ratio of [relatively] explained variance to unexplained variance
r2tof2 <- function(r2,r2baseModel=0){
	f2 <- (r2 - r2baseModel)/(1-r2)
}