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
	
	findExample : function( s ){
		'use strict'
		for( var i = 0 ; i < examples.length ; i++ ){
			if( examples[i].l.toLowerCase().indexOf(s.toLowerCase()) >= 0 ){
				return GraphParser.parseGuess(examples[i].e,examples[i].v).toDot()
			}
		}
	}
}
