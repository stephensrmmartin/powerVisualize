###########
#Libraries#
###########
require(ggplot2)
require(pwr)

#Creates array of power values given effectsize-df pairs (power of .03 with 50; .03 with 51, etc)
genPwrArray <- function(explodedEff,explodedDfs,dfNumerator,sig.level=.05){
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
createPwrFrame <- function(effectSizes,dfs,dfNumerator=4,sig.level=.05){
	exploded <- explodeEffDfs(effectSizes,dfs)
	effects <- exploded$effect
	dfs <- exploded$df
	rm(exploded)
	pwr <- genPwrArray(effects,dfs,dfNumerator,sig.level)
	ds <- data.frame(effectSize=effects,df=dfs,power=pwr)
	return(ds)	
}
#Takes a pwrFrame and cutoff value.
#Splits the frame by effect size,
#returns the first value of x that intersects with the given y-intercept
getXIntercepts <- function(pwrFrame,cutoff){
	require(plyr)
	daply(.data=pwrFrame,.(effectSize),.fun=function(x){
		return(x[x$power >= cutoff & x$power <= (cutoff + .05),"df"][1])
	})
}

#Plots the power data.frame generated from createPwrFrame()
#Will treat effectSize as factor if fewer than 8 effectsizes are given
#Otherwise, the whole thing is a scatterplot.
plotPower <- function(pwrFrame,guides=TRUE,cutoff=.8){
	if(length(unique(pwrFrame$effectSize)) < 8){
		pwrFrame$effectSize <- factor(pwrFrame$effectSize)
		p <- qplot(data=pwrFrame,x=df,y=power,color=effectSize,group=effectSize,geom="line")
		if(guides==TRUE){
			p <- p + geom_hline(yintercept=cutoff,linetype=2)
			p <- p + annotate(x=-20,y=cutoff,geom="text",label=as.character(cutoff),angle=45)
			p <- p + geom_vline(xintercept = getXIntercepts(pwrFrame,cutoff),linetype=2)
			p <- p + annotate(x=getXIntercepts(pwrFrame,cutoff),geom="text",label=as.character(getXIntercepts(pwrFrame,cutoff)),y=.1,angle=45)
		}
	}
	else{
		p <- qplot(data=pwrFrame,x=df,y=power,color=effectSize,geom="point") + scale_color_gradient2(low='red',mid='blue',high='black')
	}
	p <- p + theme_classic()
	p <- p + labs(x="Degrees of freedom",y='Power',title='Power Analysis',colour='Effect Size')
	return(p)

}