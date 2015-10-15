#' @import V8 jsonlite
NULL

#' Get Bundled Examples
#'
#' Provides access to the builtin examples of the dagitty website.
#'
#' @param x name of the example, or part thereof. Supported values are:
#' \itemize{
#'  \item{"M-bias"}{The M-bias graph}
#'  \item{"confounding"}{An extended confounding triangle}
#'  \item{"mediator"}{A small model with a mediator}
#'  \item{"paths"}{A graph with many variables but few paths}
#'  \item{"Sebastiani"}{A small part of a genetics study (Sebastiani et al., 2005)}
#'  \item{"Polzer"}{DAG from a dentistry study (Polzer et al., 2012)}
#'  \item{"Schipf"}{DAG from a study on diabetes (Schipf et al., 2010)}
#'  \item{"Shrier"}{DAG from a classic sports medicine example (Shrier & Platt, 2008)}
#'  \item{"Thoemmes"}{DAG with unobserved variables (communicated by F. Thoemmes, 2013)}
#'  \item{"Kampen"}{DAG from a psychiatry study (van Kampen, 2014)}
#' }
#' @references
#' Sabine Schipf, Robin Haring, Nele Friedrich, Matthias Nauck, Katharina Lau,
#' Dietrich Alte, Andreas Stang, Henry Voelzke, and Henri Wallaschofski (2011),
#' Low total testosterone is associated with increased risk of incident
#' type 2 diabetes mellitus in men: Results from the study of health in
#' pomerania (SHIP). \emph{The Aging Male} \bold{14}(3):168--75.
#'
#' Paola Sebastiani, Marco F. Ramoni, Vikki Nolan, Clinton T. Baldwin, and
#' Martin H. Steinberg (2005), Genetic dissection and prognostic modeling of overt 
#' stroke in sickle cell anemia. \emph{Nature Genetics}, \bold{37}:435--440.
#' 
#' Ian Shrier and Robert W. Platt (2008), 
#' Reducing bias through directed acyclic graphs.
#' \emph{BMC Medical Research Methodology}, \bold{8}(70).
#'
#' Ines Polzer, Christian Schwahn, Henry Voelzke, Torsten Mundt, and Reiner
#' Biffar (2012), The association of tooth loss with all-cause and circulatory
#'  mortality. Is there a benefit of replaced teeth? A systematic review and
#'  meta-analysis. \emph{Clinical Oral Investigations}, \bold{16}(2):333--351.
#'
#' Dirk van Kampen (2014),
#' The SSQ model of schizophrenic prodromal unfolding revised: An
#'  analysis of its causal chains based on the language of directed graphs.
#' \emph{European Psychiatry}, \bold{29}(7):437--48.
#'
#' @export
getExample <- function( x ){
	ct <- .getJSContext()
	xv <- .getJSVar()
	tryCatch({
		.jsassign( xv, x )
		r <- ct$eval( paste0("DagittyR.findExample(global.",xv,")") )
	}, finally={.deleteJSVar(xv)})
	if( r != "undefined" ){
		structure(r,class="dagitty")
	} else {
		stop("Example ",x," could not be found!")
	}
}

#' Change Status of Variables in Graph
#'
#' First removes the given status from all variables in the graph that had it, and then 
#' gives it to the given variables.
#' For instance, if  \code{status="exposure"}  and \code{var.names="X"} are given, then
#' \code{X} will be the only exposure in the resulting graph.
#'
#' @param x the input graph.
#' @param var.names names of variables to change.
#' @param status character, one of "exposure" or "outcome".
#' @export
setVariableStatus <- function( x, var.names, status ) {
	allowed.statuses <-  c("exposure","outcome")
	if( !(status %in% allowed.statuses) ){
		stop( "Status must be one of ", allowed.statuses )
	}
	
	xv <- .getJSVar()
	vv <- .getJSVar()
	tryCatch({
		.jsassign( xv, as.character(x) )
		.jsassign( xv, .jsp("GraphParser.parseGuess(global.",xv,")") )
		if( status == "exposure" ){
			.jseval( paste0( "global.",xv,".removeAllSources()" ) )
			for( n in var.names ){
				.jsassign( vv, n )
				.jseval( paste0( "global.",xv,".addSource(global.",vv,")" ) )
			}
		} else if (status == "outcome" ){
			.jseval( paste0( "global.",xv,".removeAllTargets()" ) )
			for( n in var.names ){
				.jsassign( vv, n )
				.jseval( paste0( "global.",xv,".addTarget(",vv,")" ) )
			}
		}
		r <- .jsget( paste0( xv,".toString()" ) )
	},finally={ 
		.deleteJSVar(vv)
		.deleteJSVar(xv)
	})
	structure(r,class="dagitty")
}

#' Replace Bidirected Edges by Latent Variables
#' @param x the input graph.
#' @export
canonicalGraph <- function( x ){
	xv <- .getJSVar()
	r <- list()
	tryCatch({
		.jsassign( xv, as.character(x) )
		.jsassign( xv, .jsp("GraphParser.parseGuess(global.",xv,")") )
		.jsassign( xv, .jsp("GraphTransformer.canonicalGraph(global.",xv,")") )
		r$g <- structure( .jsget( paste0(xv,".g.toString()") ), class="dagitty" )
		r$L <- .jsget( paste0(xv,".L") )
		r$S <- .jsget( paste0(xv,".S") )
	},finally={
		.deleteJSVar(xv)
	})
	r	
}

#' Get Graphical Descendants
#'
#' Finds all variables that are reachable from \eqn{v} via a directed path. By definition,
#' this includes \eqn{v} itself.
#'
#' @param x the input graph.
#' @param v character, name of a variable. 
#' @export
descendants <- function( x, v ){
	.kins( x, v, "descendants" )
}

#' Get Graphical Ancestors
#'
#' Finds all variables from which \eqn{v} is reachable via a directed path. 
#' By definition, this includes \eqn{v} itself.
#'
#' @param x the input graph.
#' @param v character, name of a variable. 
#' @export
ancestors <- function( x, v ){
	.kins( x, v, "ancestors" )
}

#' Get Graphical Children
#'
#' Finds all variables \eqn{w} that are connected to \eqn{v} 
#' by an edge \eqn{v} -> \eqn{w}. 
#'
#' @param x the input graph.
#' @param v character, name of a variable.
#' @export
children <- function( x, v ){
	.kins( x, v, "children" )
}

#' Get Graphical Parents
#'
#' Finds all variables \eqn{w} to which \eqn{v} is connected
#' by an edge \eqn{w} -> \eqn{v}. 
#'
#' @param x the input graph.
#' @param v character, name of a variable.
#' @export
parents <- function( x, v ){
	.kins( x, v, "parents" )
}

#' Get Exposure Variables
#' @param x the input graph.
#' @export
exposures <- function( x ){
	.nodesWithProperty( x, "source" )
}

#' Get Outcome Variables
#' @param x the input graph.
#' @export
outcomes <- function( x ){
	.nodesWithProperty( x, "target" )
}

#' Names of Variables in Graph
#' @param x the input graph.
#' @export
#' @examples
#' ## A "DAG" with Greek and Portuguese variable names. These are input using 
#' ## URL-encoding.
#' g <- dagitty( "graph {
#'   %CE%BA%CE%B1%CF%81%CE%B4%CE%AF%CE%B1 [pos=\"0.297,0.502\"]
#'   cora%C3%A7%C3%A3o [pos=\"0.482,0.387\"]
#'   %CE%BA%CE%B1%CF%81%CE%B4%CE%AF%CE%B1 -> cora%C3%A7%C3%A3o
#' }" )
#' names( g )
names.dagitty <- function( x ){
	ct <- .getJSContext()
	xv <- .getJSVar()
	tryCatch({
		.jsassign( xv, as.character(x) )
		.jsassign( xv, .jsp("GraphParser.parseGuess(global.",xv,")") )
		.jsassign( xv, .jsp("global.",xv,".vertices.keys()") )
		r <- .jsget(xv)},
	finally={.deleteJSVar(xv)})
	r
}

#' Plot Coordinates of Variables in Graph
#'
#' The DAGitty syntax allows specification of plot coordinates for each variable in a 
#' graph. This function extracts these plot coordinates from the graph description in a
#' \code{dagitty} object. Note that the coordinate system is undefined, typically one 
#' needs to compute the bounding box before plotting the graph.
#'
#' @param x the input graph.
#' @return A list with components \code{x} and \code{y}, giving relative coordinates 
#' for each variable.
#' 
#' @examples
#' ## Plot localization of each node in the Shier example
#' plot( coordinates( getExample("Shrier") ) )
#' @seealso
#' Function \link{layout} for automtically generating layout coordinates, and function
#' \link{plot.dagitty} for plotting graphs.
#'
#' @export
coordinates <- function( x ){
	ct <- .getJSContext()
	xv <- .getJSVar()
	yv <- .getJSVar()
	tryCatch({
		.jsassign( xv, as.character(x) )
		.jsassign( xv, .jsp("GraphParser.parseGuess(global.",xv,").vertices.values()") )
		.jsassign( yv, .jsp("_.pluck(global.",xv,",'id')") )
		labels <- .jsget(yv)
		.jsassign( yv, .jsp("_.pluck(global.",xv,",'layout_pos_x')") )
		rx <- .jsget(yv)
		.jsassign( yv, .jsp("_.pluck(global.",xv,",'layout_pos_y')") )
		ry <- .jsget(yv)},
	finally={
		.deleteJSVar(xv)
		.deleteJSVar(yv)
	})
	names(rx) <- labels
	names(ry) <- labels
	list( x=rx, y=ry )
}

#' Graph Edges
#' 
#' Extracts edge information from the input graph. 
#'
#' @param x the input graph.
#' @return a data frame with the following variables:
#' \itemize{
#'  \item{v}{Name of the start node.}
#'  \item{w}{Name of the end node. For symmetric edges (bidirected and undirected), the
#'  order of start and end node is arbitrary.}
#'  \item{e}{Type of edge. Can be one of \code{"->"}, \code{"<->"} and \code{"--"}}
#'  \item{x}{X coordinate for a control point. If this is not \code{NA}, then the edge
#'  is drawn as an \link{xspline} through the start point, this control point, and the 
#'  end point. This is especially important for cases where there is more than one edge
#'  between two variables (for instance, both a directed and a bidirected edge).}
#'  \item{y}{Y coordinate for a control point.}
#' }
#'
#' @examples
#' ## Which kinds of edges are used in the Shrier example?
#' levels( edges( getExample("Shrier") )$e )
#' @export 
edges <- function( x ){
	xv <- .getJSVar()
	tryCatch({
		.jsassign( xv, as.character(x) )
		.jsassign( xv, .jsp("GraphParser.parseGuess(global.",xv,")") )
		.jsassign( xv, .jsp("DagittyR.edge2r(global.",xv,")") )
		r <- .jsget(xv)
	}, finally={.deleteJSVar(xv)})
	as.data.frame(r)
}

#' Test for Graph Class
#' 
#' A function to check whether an object has class \code{dagitty}.
#'
#' @param x object to be tested.
#' 
#' @export
is.dagitty <- function(x) inherits(x,"dagitty")

#' Generate Graph Layout
#'
#' This function generates plot coordinates for each variable in a graph that does not
#' have them already. To this end, the well-known \dQuote{Spring} layout algorithm is
#' used. Note that this is a stochastic algorithm, so the generated layout will be 
#' different every time (which also means that you can try several times until you find
#' a decent layout).
#' 
#' @param x the input graph.
#' @param method the layout method; currently, only \code{"spring"} is supported.
#' @return the same graph as \code{x} but with layout coordinates added. 
#' 
#' @examples
#' ## Generate a layout for the M-bias graph and plot it
#' plot( layout( dagitty( "graph { X <- U1 -> M <- U2 -> Y } " ) ) )
#'
#' @export
layout <- function( x, method="spring" ){
	if( !(method %in% c("spring")) ){
		stop("Layout method ",method," not supported!")
	}
	if( !is.dagitty( x ) ){
		stop("Expecting a graph as input x!")
	}
	xv <- .getJSVar()
	tryCatch({
		.jsassign( xv, as.character(x) )
		.jsassign( xv, .jsp("GraphParser.parseGuess(global.",xv,")") )
		.jseval( paste0("(new GraphLayouter.Spring(global.",xv,")).layout()") )
		.jsassign( xv, .jsp("global.",xv,".toString()") )
		r <- .jsget(xv)
	}, finally={.deleteJSVar(xv)})
	structure( r, class="dagitty" )
}

#' Plotting Graph
#'
#' A simple plot method to quickly visualize a graph. This is intended mainly for 
#' validation purposes and is not yet meant to become a full-fledged graph drawing
#' function.
#'
#' @param x the input graph.
#' @param ... not used.
#'
#' @export
plot.dagitty <- function( x, ... ){	
	if( !is.dagitty( x ) ){
		stop("Expecting a graph as input x!")
	}
	coords <- coordinates( x )
	labels <- names(coords$x)
	plot.new()
	par(new=TRUE,mar=rep(0,4))
	wx <- sapply( paste0("mm",labels), 
		function(s) strwidth(s,units="inches") )
	wy <- sapply( paste0("\n",labels), 
		function(s) strheight(s,units="inches") )
	ppi.x <- dev.size("in")[1] / (max(coords$x)-min(coords$x))
	ppi.y <- dev.size("in")[2] / (max(coords$y)-min(coords$y))
	wx <- wx/ppi.x
	wy <- wy/ppi.y
	xlim <- c(min(coords$x-wx/2),max(coords$x+wx/2))
	ylim <- c(-max(coords$y+wy/2),-min(coords$y-wy/2))
	plot( NA, xlim=xlim, ylim=ylim, xlab="", ylab="", bty="n",
		xaxt="n", yaxt="n" )
	wx <- sapply( labels, 
		function(s) strwidth(paste0("xx",s)) )
	wy <- sapply( labels,
		function(s) strheight(paste0("\n",s)) )
	asp <- par("pin")[1]/diff(par("usr")[1:2]) /
		(par("pin")[2]/diff(par("usr")[3:4]))
	ex <- edges(x)
	ax1 <- rep(0,nrow(ex))
	ax2 <- rep(0,nrow(ex))
	ay1 <- rep(0,nrow(ex))
	ay2 <- rep(0,nrow(ex))
	axc <- rep(0,nrow(ex))
	ayc <- rep(0,nrow(ex))
	acode <- rep(2,nrow(ex))
	has.control.point <- rep(FALSE,nrow(ex))
	for( i in seq_len(nrow(ex)) ){
		if( ex[i,3] == "<->" ){
			acode[i] <- 3
			has.control.point[i] <- TRUE
		}
		if( ex[i,3] == "--" ){
			acode[i] <- 0
		}
		l1 <- as.character(ex[i,1]); l2 <- as.character(ex[i,2])
		x1 <- coords$x[l1]; y1 <- coords$y[l1]
		x2 <- coords$x[l2]; y2 <- coords$y[l2]
		if( is.na( ex[i,4] ) || is.na( ex[i,5] ) ){
			cp <- .autoControlPoint( x1, y1, x2, y2, asp,
				.2*as.integer( acode[i]==3 ) )
		} else {
			cp <- list(x=ex[i,4],y=ex[i,5])
			has.control.point[i] <- TRUE
		}
		bi1 <- .lineSegBoxIntersect( x1-wx[l1]/2,y1-wy[l1]/2,
			x1+wx[l1]/2,y1+wy[l1]/2, x1, y1, cp$x, cp$y )
		bi2 <- .lineSegBoxIntersect( x2-wx[l2]/2,y2-wy[l2]/2,
			x2+wx[l2]/2,y2+wy[l2]/2, cp$x, cp$y, x2, y2 )
		if( length(bi1) == 2 ){
			x1 <- bi1$x; y1 <- bi1$y
		}
		if( length(bi2) == 2 ){
			x2 <- bi2$x; y2 <- bi2$y
		}
		ax1[i] <- x1; ax2[i] <- x2
		ay1[i] <- y1; ay2[i] <- y2
		axc[i] <- cp$x; ayc[i] <- cp$y
	}
	directed <- acode==2 & !has.control.point
	undirected <- acode==0 & !has.control.point
	arrows( ax1[directed], -ay1[directed], 
		ax2[directed], -ay2[directed], length=0.1, col="gray" )
	segments( ax1[undirected], -ay1[undirected], 
		ax2[undirected], -ay2[undirected], col="gray" )
	for( i in which( has.control.point ) ){
		.arc( ax1[i], -ay1[i], 
			ax2[i], -ay2[i], axc[i], -ayc[i], col="gray", 
			code=acode[i], length=0.1 )
	}
	text( coords$x, -coords$y[labels], labels )
}

#' Find Adjustment Sets to Estimate Causal Effects
#' @param x the input graph.
#' @param exposure name(s) of the exposure variable(s). If not given (default), then the 
#'  exposure variables are supposed to be defined in the graph itself.
#' @param outcome name(s) of the outcome variable(s), also taken from the graph if 
#' not given.
#' @param effect which effect is to be identified. If \code{effect="total"}, then the
#' total effect is to be identified, and the adjustment criterion by van der Zander et 
#' al (2014), an extension of Pearl's back-door criterion, is used. Otherwise, if 
#' \code{effect="direct"}, then the average direct effect is to be identified, and Pearl's
#' single-door criterion is used (Pearl, 2009). In a structural equation model (Gaussian
#' graphical model), direct effects are simply the path coefficients.
#' @export
adjustmentSets <- function( x, exposure=NULL, outcome=NULL, effect="total" ){

	if( !is.null( exposure ) ){
		x <- setVariableStatus( x, exposure, "exposure" )
	}
	if( !is.null( outcome ) ){
		x <- setVariableStatus( x, outcome, "outcome" )
	}

	xv <- .getJSVar()
	tryCatch({
		.jsassign( xv, as.character(x) )
		.jsassign( xv, .jsp("GraphParser.parseGuess(global.",xv,")") )
	
		if( effect=="direct" ){	
			.jsassign( xv, .jsp("GraphAnalyzer.listMsasDirectEffect(global.",xv,")") )
		} else {
			.jsassign( xv, .jsp("GraphAnalyzer.listMsasTotalEffect(global.",xv,")") )
		}
		.jsassign( xv, .jsp("DagittyR.adj2r(global.",xv,")"))
		r <- structure( .jsget(xv), class="dagitty.sets" )
	},finally={.deleteJSVar(xv)})

	r
}

#' List Conditional Indpendencies Implied by Graphical Model
#' @param x the input graph.
#' @param max.results integer. The listing of conditional independencies is stopped once
#' this many results have been found. Use \code{Inf} to generate them all.
#' @export
impliedConditionalIndependencies <- function( x, max.results=100 ){
	xv <- .getJSVar()
	tryCatch({
		.jsassign( xv, as.character(x) )
		.jsassign( xv, .jsp("GraphParser.parseGuess(global.",xv,")") )
		if( is.finite( max.results ) ){
			.jsassign( xv, .jsp("GraphAnalyzer.listMinimalImplications(global.",xv,",",
				as.numeric(max.results),")"))
		} else {
			.jsassign( xv, .jsp("GraphAnalyzer.listMinimalImplications(global.",xv,")"))
		}
		.jsassign( xv, .jsp("DagittyR.imp2r(global.",xv,")") )
		r <- structure( .jsget(xv), class="dagitty.cis" )
	},finally={.deleteJSVar(xv)})
	
	r
}

#' Find (Conditional) Instrumental Variables to Estimate Causal Effects
#' @param x the input graph.
#' @param exposure name of the exposure variable. If not given (default), then the 
#'  exposure variable is supposed to be defined in the graph itself. Only a single exposure
#' variable is supported.
#' @param outcome name of the outcome variable, also taken from the graph if not given.
#' Only a single outcome variable is supported.
#' @export
instrumentalVariables <- function( x, exposure=NULL, outcome=NULL ){
	if( !is.null( exposure ) ){
		if( length( exposure ) > 1 ){
			stop("IV identification only supported for single exposures!")
		}
		x <- setVariableStatus( x, exposure, "exposure" )
	}
	if( !is.null( outcome ) ){
		if( length( outcome ) > 1 ){
			stop("IV identification only supported for single outcomes!")
		}
		x <- setVariableStatus( x, outcome, "outcome" )
	}

	xv <- .getJSVar()
	tryCatch({
		.jsassign( xv, as.character(x) )
		.jsassign( xv, .jsp("GraphParser.parseGuess(global.",xv,")") )
		.jsassign( xv, .jsp("GraphAnalyzer.conditionalInstruments(global.",xv,")") )
		.jsassign( xv, .jsp("DagittyR.iv2r(global.",xv,")") )
		r <- structure( .jsget(xv), class="dagitty.ivs" )
	}, finally={.deleteJSVar(xv)})
	r
}

#' List Vanishing Tetrads Implied by Gaussian Graphical Model
#' @param x the input graph.
#' @export
vanishingTetrads <- function( x ){
	xv <- .getJSVar()
	tryCatch({
		.jsassign( xv, as.character(x) )
		.jsassign( xv, .jsp("GraphParser.parseGuess(global.",xv,")") )
		.jsassign( xv, .jsp("GraphAnalyzer.vanishingTetrads(global.",xv,")") )
		r <- .jsget(xv)
	}, finally={.deleteJSVar(xv)})
	r
}

#' Convert Lavaan Model to Dagitty Graph
#' @param x data frame, lavaan parameter table such as returned by 
#' \code{\link[lavaan]{lavaanify}}.
#' @export
as.dagitty.data.frame <- function( x ){
	latents <- c()
	arrows <- c()
	for( i in seq_len( nrow(x) ) ){
		if( x$op[i] == "=~" ){
			latents <- union(latents,x$lhs[i])
			arrows <- c(arrows,paste(x$lhs[i]," -> ",x$rhs[i]))
		}
		if( x$op[i] == "~" ){
			arrows <- c(arrows,paste(x$lhs[i]," <- ",x$rhs[i]))
		}
		if( x$op[i] == "~~" && (x$lhs[i] != x$rhs[i]) ){
			arrows <- c(arrows,paste(x$lhs[i]," <-> ",x$rhs[i]))
		}
	}
	dagitty( paste("graph { ",
		paste(latents,' [latent]',collapse="\n"),"\n",
		paste(arrows,collapse="\n")," } ",collapse="\n") )
}

#' Parse Dagitty Graph
#' @param x character, string describing a graphical model in dagitty syntax.
#' @export
dagitty <- function(x){
	if(!is.character(x)){
		stop("Expecting a string to build dagitty graph!")
	}
	xv <- .getJSVar()
	tryCatch({
		.jsassign( xv, as.character(x) )
		.jsassign( xv, .jsp("GraphParser.parseGuess(global.",xv,").toString()") )
		r <- structure( .jsget(xv), class="dagitty" )
	}, finally={.deleteJSVar(xv)})
	structure( r, class="dagitty" )
}

#' Load Graph drom dagitty.net
#'
#' @param x dagitty model URL
#' @export
downloadGraph <- function(x="dagitty.net/mz-Tuw9"){
	if( !requireNamespace( "base64enc", quietly=TRUE ) ){
		stop("This function requires the package 'base64enc'!")
	}
	id <- gsub( "dagitty\\.net\\/m(.*)$", "\\1", x )
	r <- base64enc::base64decode(scan(paste0("http://dagitty.net/dags/load.php?id=",id),"character"))
	if( base64enc::checkUTF8(r) ){
		dagitty( rawToChar( r ) )
	} else {
		NULL
	}
}

#' @export
print.dagitty.ivs <- function( x, prefix="", ... ){
	for( i in x ){
		cat( prefix, i$I )
		if( length( i$Z > 0 ) ){
			cat( " | ", paste(i$Z,collapse=", ") )
		}
		cat( "\n" )
	}
}

#' @export
print.dagitty.sets <- function( x, prefix="", ... ){
	for( i in x ){
		if( length(i) == 0 ){
			cat( prefix, "{}\n")
		} else {
			cat( prefix, "{",paste(i,collapse=", "), "}\n" )
		}
	}
}

#' @export
print.dagitty.cis <- function( x, ... ){
	for( i in x ){
		cat( i$X, " _||_ ", i$Y )
		if( length( i$Z > 0 ) ){
			cat( " | ", paste(i$Z,collapse=", ") )
		}
		cat("\n")
	}
}
