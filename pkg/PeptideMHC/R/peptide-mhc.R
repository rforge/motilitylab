
#' Supported MHC Molecules and Peptide Lengths.
#'
#' @return A named list. The names are the MHC molecules and the list entries
#' are integer vectors containing the supported peptide lengths for each MHC.
#' @examples
#' ## Which MHC molecules are supported?
#' names( supportedMHCs() )
#' ## Which peptide lengths are supported for HLA-A02:01?
#' supportedMHCs()[["HLA-A02:01"]]
#'
#' @export
supportedMHCs <- function(){
	r <- read.table( system.file( "extdata", "model_list.txt", package="PeptideMHC" ),
		as.is=TRUE )$V1
	mhcname <- gsub( '-[0-9][0-9]?$', '', r )
	mhcname <- gsub( '-([0-9][0-9])([0-9][0-9])', '\\1:\\2', mhcname )
	peptidelength <- as.integer(gsub( ".*-([0-9][0-9]?)$", "\\1", r ))
	split( peptidelength, mhcname )
}

.smm.cache <- new.env()

.readSMMmatrix <- function( mhc="HLA-A02:01", l=9 ){
	matrixname <- paste0(gsub(":","",mhc),"-",l)
	if( exists( matrixname, .smm.cache ) ){
		return( get( matrixname, envir=.smm.cache ) )
	}
	matrixfile <- system.file( "extdata", paste0(matrixname,".txt"), package="NetMHCpan" )
	if( matrixfile == "" ){
		stop( paste( "SMM matrix for MHC",mhc,"and peptide length",l,"not found!") )
	} else {
		M <- as.matrix( read.table(matrixfile, skip=1, nrows=20, row.names=1) )
		colnames(M) <- NULL
		assign( paste0(matrixname,"-c"),  as.numeric( read.table(matrixfile, skip=21) ),
			.smm.cache )
		assign( matrixname, M, .smm.cache )
		return( M )
	}
}

.smmConstant <- function( mhc, l ){
	get( paste0(gsub(":","",mhc),"-",l,"-c"), .smm.cache )
}

#' Peptide-MHC Binding Prediction
#'
#' Predicts peptide-MHC binding using the stabilized matrix method (SMM) 
#' algorithm with a specifically constructed amino acid substitution matrix. 
#' See \code{citation(package='PeptideMHC')} for the reference.
#'
#' @param peptides vector of strings containing the peptides for which to predict
#' binding.
#' @param mhc string or vector of strings identifying the MHC molecules. See
#' \code{supportedMHCs} for allowed values.
#' @param output.IC50 whether to output the IC50 value itself (default) or its 
#' base10-logarithm.
#' @export
smmPMBEC <- function( peptides=c("SLYNTVATL","SYFPEITHI"), mhc="HLA-A02:01", 
	output.IC50=TRUE ){
	if( length(mhc)>1 && ( length(mhc) != length(peptides) ) ){
		stop( "If 'mhc' is a vector, it must have the same length as 'peptides'!" )
	}
	if( length(mhc) == 1 && length(peptides) > 1 ){
		mhc <- rep.int( mhc, length( peptides ) )
	}
	pred <- function( i ){
		l <-  nchar(peptides[i])
		M <- .readSMMmatrix( mhc[i] )
		.smmConstant( mhc[i], l ) + 
			sum( sapply( seq_len(l), function(j) M[substring(peptides[i],j,j),j] ) )
	}
	v <- sapply( seq_along( peptides ), pred )
	if( output.IC50 ){
		return( 10^v )
	} else {
		return( v )
	}
}

smmNonAnchorAAs <- function( peptides=c("SLYNTVATL","SYFPEITHI"), mhc="HLA-A02:01", n=3 ){
	if( length(mhc)>1 && ( length(mhc) != length(peptides) ) ){
		stop( "If 'mhc' is a vector, it must have the same length as 'peptides'!" )
	}
	if( length(mhc) == 1 && length(peptides) > 1 ){
		mhc <- rep.int( mhc, length( peptides ) )
	}
	pred <- function( i ){
		M <- .readSMMmatrix( mhc[i] )
		unc <- apply( M, 2, sd )
		anchors <- c(which( unc >= sort(unc,decreasing=TRUE)[3] ))
		paste( strsplit( peptides[i], "" )[[1]][-anchors], collapse='' )
	}
	sapply( seq_along( peptides ), pred )
}

