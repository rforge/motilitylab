/** This file contains some utility function that convert the results of important 
	DAGitty functions into a format that is more easily mapped to an R list-of-lists. */

/* globals _, examples, GraphParser  */
/* exported DagittyR */

var DagittyR = {
	adj2r : function( ss ){
		'use strict'
		var r = {}, i=1
		if( ss.length === 0 )
			return r
		_.each( ss, function(s){
			var rs = _.pluck( s, 'id').sort()
			r[i++]=rs
		})
		return r
	},

	imp2r : function( imp ){
		'use strict'
		var r = {}, k=1, rk
		_.each( imp, function( i ){
			_.each( i[2], function( z ){
				rk = { X: i[0], Y: i[1] }
				rk.Z = _.pluck( z, 'id').sort()
				r[k++]=rk
			})
		})
		return r
	},

	iv2r : function( ivs ){
		'use strict'
		var r = {}, j=1, rk
		_.each(ivs, function( i ){
			rk = { "I": i[0].id }
			if( i[1].length > 0 ){
				rk.Z = _.pluck(i[1],'id').sort()
			} 
			r[j++]=rk
		} )
		return r
	},
	
	edge2r : function( g ){
		'use strict'
		var r = { v : [], w : [], e : [], x : [], y : [] }, j=1, rk
		_.each(g.edges, function( e ){
			r.v.push( e.v1.id )
			r.w.push( e.v2.id )
			r.x.push( e.layout_pos_x )
			r.y.push( e.layout_pos_y )
			switch( e.directed ){
			case Graph.Edgetype.Directed :
				r.e.push("->"); break
			case Graph.Edgetype.Undirected :
				r.e.push("--"); break
			case Graph.Edgetype.Bidirected :
				r.e.push("<->"); break
			}
		} )
		return r
	},
	
	paths2r : function( ga, Z ){
		var i = 0, r = { path : [], open : [] }
		for( ; i < ga.length ; i ++ ){
			r.path.push( GraphSerializer.pathToDot( ga[i] ) )
			r.open.push( !GraphAnalyzer.dSeparated( ga[i], 
				ga[i].getSources(), ga[i].getTargets(), [] ) )
		}
		return r
	},
	
	getVertices : function( g, a ){
		var r = [], v, ak = Object.keys( a )
		for( var i = 0 ; i < a.length ; i ++ ){
			v = g.getVertex( a[ak[i]] )
			if( v ){
				r.push(v)
			}
		}
		return r
	},
	
	findExample : function( s ){
		'use strict'
		for( var i = 0 ; i < examples.length ; i++ ){
			if( examples[i].l.toLowerCase().indexOf(s.toLowerCase()) >= 0 ){
				return GraphParser.parseGuess(examples[i].e,examples[i].v).toString()
			}
		}
	}
}
