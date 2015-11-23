small1 <- dagitty('graph {A -> S -> T}')
exposures(small1) <- 'S'
outcomes(small1) <- 'T'

# covers old unittests 1 t/m 3

test_that("parsing / serializing", {
	expect_equal(toString(small1,"dagitty.old"), "A 1\nS E\nT O\n\nA S\nS T")	
})

test_that("ancestor moral graph", {
	expect_equal(toString(ancestorGraph(small1),"dagitty.old"), 
		"A 1\nS E\nT O\n\nA S\nS T")
	expect_equal(toString(moralize(ancestorGraph(small1)),"dagitty.old"), 
		"A 1\nS E\nT O\n\nA S\nS A T\nT S")
	expect_equal(nrow(edges(moralize(dagitty("graph{ a -> x <-> z <- b}")))),6)
})
