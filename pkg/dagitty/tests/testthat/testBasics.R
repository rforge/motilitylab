small1 <- dagitty('graph {A -> S -> T}')
exposures(small1) <- 'S'
outcomes(small1) <- 'T'

mixed <- dagitty("graph{ a -- x -> b \n c <-> x <- d \n f }")

# covers old unittests 1 t/m 3

test_that("parsing / serializing", {
	expect_equal(toString(small1,"dagitty.old"), "A 1\nS E\nT O\n\nA S\nS T")	
})

test_that("relationships", {
	expect_equal(spouses(mixed,"x"),("c"))
	expect_equal(neighbours(mixed,"x"),("a"))
	expect_equal(parents(mixed,"x"),("d"))
	expect_equal(children(mixed,"x"),("b"))
	expect_equal(sort(adjacentNodes(mixed,"x")),c("a","b","c","d"))
	expect_equal(adjacentNodes(mixed,"f"),list())
})

test_that("ancestor moral graph", {
	expect_equal(toString(ancestorGraph(small1),"dagitty.old"), 
		"A 1\nS E\nT O\n\nA S\nS T")
	expect_equal(toString(moralize(ancestorGraph(small1)),"dagitty.old"), 
		"A 1\nS E\nT O\n\nA S\nS A T\nT S")
	expect_equal(nrow(edges(moralize(dagitty("graph{ a -> x <-> z <- b}")))),6)
})
