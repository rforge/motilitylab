library( MotilityLab )
library( beeswarm )

par( bty="L" )

par( mfrow=c(3,3) )

tracks <- c( BCells, TCells )
labels <- rep( c("BCells","TCells"), c(length(BCells),length(TCells) ) )

beebx <- function( x, g, ylab ){
	beeswarm( x ~ g, ylim=c(0,max(x)),
		xlab="", ylab=ylab )
	bxplot( x ~ g, add=TRUE, width=0.5 )
}

beebx( sapply( tracks, speed ), labels,
	expression("mean track speed ("*mu*"m/sec)")  )

beebx( sapply( tracks, straightness ), labels,
	expression("straightness index")  )

with( aggregate(TCells,squareDisplacement,FUN="mean.se",max.overlap=0), {
	plot( mean ~ i,
		xlab="time", ylab="mean square displacement", col=2, type='l' )
	segments( i, lower, y1=upper, col=2 )
} )
with( aggregate(BCells,squareDisplacement,FUN="mean.se",max.overlap=0), {
	lines( mean ~ i, col=3 )
	segments( i, lower, y1=upper, col=3 )
} )

with( aggregate(TCells,displacement,FUN="mean.se",max.overlap=0), {
	plot( mean ~ sqrt(i),
		xlab=expression("time"^"1/2"), ylab="mean square displacement", col=2, type='l' )
	segments( sqrt(i), lower, y1=upper, col=2 )
} )
with( aggregate(BCells,displacement,FUN="mean.se",max.overlap=0), {
	lines( mean ~ sqrt(i), col=3 )
	segments( sqrt(i), lower, y1=upper, col=3 )
} )

with( (aggregate(TCells,overallAngle,FUN="mean.se",na.rm=TRUE)),{
	plot( cos(mean) ~ i, ylim=c(min(0,min(cos(mean))),1), xlab="time", 
		ylab="correlation of orientation", type="l", col=2 )
	segments( i, cos(lower), y1=cos(upper), col=2 )
})

with( (aggregate(BCells,overallAngle,FUN="mean.se",na.rm=TRUE)),{
	lines( cos(mean) ~ i, xlab="time", 
		ylab="correlation of orientation", type="l", col=3 )
	segments( i, cos(lower), y1=cos(upper), col=3 )
})

with( (aggregate(Neutrophils,overallAngle,FUN="mean.se",na.rm=TRUE)),{
	lines( cos(mean) ~ i, xlab="time", 
		ylab="correlation of orientation", type="l", col=4 )
	segments( i, cos(lower), y1=cos(upper), col=4 )
})

dx <- function(d) diff(d[,"x"])
plot( t(sapply( subtracks( TCells, 2 ), dx ) ), col=2 )
plot( t(sapply( subtracks( BCells, 2 ), dx ) ), col=3 )
plot( t(sapply( subtracks( Neutrophils, 2 ), dx ) ), col=4 )
