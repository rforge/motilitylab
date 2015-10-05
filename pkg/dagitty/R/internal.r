
.dagitty.cache <- new.env()

.getJSContext <- function(){
	if( !exists("ct",.dagitty.cache) ){
		requireNamespace("V8",quietly=TRUE)
		ct <- V8::new_context()
		ct$source(system.file("js/underscore-min.js",package="dagitty"))
		ct$source(system.file("js/Class.js",package="dagitty"))
		ct$source(system.file("js/Hash.js",package="dagitty"))
		ct$source(system.file("js/Graph.js",package="dagitty"))
		ct$source(system.file("js/GraphParser.js",package="dagitty"))
		ct$source(system.file("js/GraphTransformer.js",package="dagitty"))
		ct$source(system.file("js/GraphAnalyzer.js",package="dagitty"))
		ct$source(system.file("js/RUtil.js",package="dagitty"))
		ct$source(system.file("js/example-dags.js",package="dagitty"))
		assign("ct",ct,.dagitty.cache)
	}
	get("ct",.dagitty.cache)
}

.getJSVar <- function(){
	if( !exists("nvars",.dagitty.cache) ){
		assign("nvars",0,.dagitty.cache)
		with( .dagitty.cache, vars.available <- list() )
	}
	va <- get("vars.available",.dagitty.cache)
	if( length(va)>0 ){
		r <- names(tail(va,1))
		.dagitty.cache$vars.available[[r]] <- NULL
		return(r)
	} else {
		.dagitty.cache$nvars <- .dagitty.cache$nvars+1
		return(paste0("y",get("nvars",.dagitty.cache)))
	}
}

.deleteJSVar <- function(x){
	if( !exists("nvars",.dagitty.cache) ){
		assign("nvars",0,.dagitty.cache)
		with( .dagitty.cache, vars.available <- list() )
	}
	.getJSContext()$eval( paste0("delete global.",x) )
	if( is.null(.dagitty.cache$vars.available[[x]]) ){
		.dagitty.cache$vars.available[[x]] <- 1
	}
}

.jsassign <- function(name, value, auto_unbox = TRUE, ...){
    stopifnot(is.character(name))
    requireNamespace("jsonlite",quietly=TRUE)
    ct <- .getJSContext()
    obj <- if (any(is(value, "JS_EVAL"), is(value, "AsIs"))) {
        invisible(ct$eval(paste("global.", name, "=", value)))
    }
    else {
        invisible(ct$eval(paste("global.", name, "=", jsonlite::toJSON(value, 
            auto_unbox = auto_unbox, ...))))
    }
}

.jsget <- function (name, ...) 
{
    stopifnot(is.character(name))
    requireNamespace("jsonlite",quietly=TRUE)
    ct <- .getJSContext()
    jsonlite::fromJSON(ct$eval(c("JSON.stringify(global.", name, ")")), 
        ...)
}

.jsglobals <-function (){
	setdiff( .getJSContext()$get("Object.keys(global)"),
		c("console","print","global","ArrayBuffer","Int8Array","Uint8Array","Int16Array",
		"Uint16Array","Int32Array","Uint32Array","Float32Array","Float64Array","DataView"
		) )
}

.jseval <-function (...){
	.getJSContext()$eval(...)
}

.jsp <- function (...){
	requireNamespace("V8",quietly=TRUE)
	V8::JS( paste0( ... ) )
}

.kins <- function( x, v, type="descendants" ){
	supported <- c("descendants","ancestors","neighbours","posteriors","anteriors",
		"children","parents")
	if( ! type %in% supported ){
		stop("Supported kinship types : ",paste(supported,collapse=", ") )
	}
	r <- c()
	xv <- .getJSVar()
	vv <- .getJSVar()
	tryCatch({
		.jsassign( xv, as.character(x) )
		.jsassign( xv, .jsp("GraphParser.parseGuess(global.",xv,")") )
		
		for( w in v ){
			.jsassign( vv, as.character(w) )
			.jsassign( vv, .jsp("global.",xv,".getVertex(global.",vv,")") )
			.jsassign( vv, .jsp("global.",xv,".",type,"Of([global.",vv,"])") )
			r <- union(r, .jsget( paste0("_.pluck(global.",vv,",'id')") ))
		}
	},finally={
		.deleteJSVar(xv)
		.deleteJSVar(vv)
	})
	r
}

.capitalizeFirst <- function(s){
	paste0(toupper(substring(s, 1, 1)), substring(s, 2))
}

.nodesWithProperty <- function( x, type="source" ){
	supported <- c("source","target")
	if( ! type %in% supported ){
		stop("Supported properties : ",paste(supported,collapse=", ") )
	}
	xv <- .getJSVar()
	tryCatch({
		.jsassign( xv, as.character(x) )
		.jsassign( xv, .jsp("GraphParser.parseGuess(global.",xv,")") )
		.jsassign( xv, .jsp("global.",xv,".get",.capitalizeFirst(type),"s()") )
		r <- .jsget( paste0("_.pluck(global.",xv,",'id')") )
	},finally={
		.deleteJSVar(xv)
	})
	r
}
