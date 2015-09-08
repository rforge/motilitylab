library(MotilityLab)

cols <- c('red','orange','blue')

TCells <- projectDimensions( TCells, c("x","y") )
BCells <- projectDimensions( BCells, c("x","y") )
Neutrophils <- projectDimensions( Neutrophils, c("x","y") )

pl <- function( f, ylim=c(0,1), filter.subtracks=NULL ){
	ft <- aggregate( TCells, f, subtrack.length=1:20, na.rm=TRUE, 
		filter.subtracks=filter.subtracks )[,2]
	fb <- aggregate( BCells, f, subtrack.length=1:20, na.rm=TRUE,
		filter.subtracks=filter.subtracks )[,2]
	fn <- aggregate( Neutrophils, subtrack.length=1:20, f, na.rm=TRUE,
		filter.subtracks=filter.subtracks )[,2]
	
	plot( ft, ylim=c(0,max(c(ft,fb,fn))), type='l', 
		xlab="subtrack length (positions)", col=cols[1],
		ylab=deparse( substitute( f ) ) )
	lines( fb, ylim=c(0,1), col=cols[2] )
	lines( fn, ylim=c(0,1), col=cols[3] )
}

par( bty="L", mfrow=c(4,4) )

ch <- function(t) nrow(t)<3 || all( sapply( subtracks( t, 1, -nrow(t)+3 ), 
	trackLength ) > 2.0 )

pl( asphericity )
pl( displacement )
pl( displacementRatio )
pl( squareDisplacement )
pl( straightness )
pl( maxDisplacement )
pl( outreachRatio )
pl( overallAngle, filter.subtracks=ch )
pl( overallDot )
pl( speed )
pl( straightness )

## fractal dimension
## hurst exponent