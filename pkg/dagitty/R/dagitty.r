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
	structure(r,class="dagitty")
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
		r <- .jsget( paste0( xv,".toDot()" ) )
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
		r$g <- structure( .jsget( paste0(xv,".g.toDot()") ), class="dagitty" )
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

#' Get Names of All Variables
#' @param x the input graph.
#' @export
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
		.jsassign( xv, .jsp("GraphParser.parseGuess(global.",xv,").toDot()") )
		r <- structure( .jsget(xv), class="dagitty" )
	}, finally={.deleteJSVar(xv)})
	structure( r, class="dagitty" )
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
