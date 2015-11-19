test_that("maxTrackLength works",{
	expect_equal( maxTrackLength(TCells), 39 )
	expect_equal( maxTrackLength(BCells), 39 )
	expect_equal( maxTrackLength(Neutrophils), 55 )
})