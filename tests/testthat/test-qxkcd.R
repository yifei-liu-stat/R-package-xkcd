test_that("qxkcd gives the right value for in-domain p", {
    expect_equal(qxkcd(0), 0)
    expect_equal(qxkcd(1), maxy(1))
})

test_that("qxkcd broadcasts argument dimensions in the expected way", {
    expect_equal(length(pxkcd(0:10/10)), 11)
    expect_equal(length(pxkcd(1:2/10, 1:10)), 10)
    expect_equal(length(pxkcd(1:10/10, 1:10)), 10)
})

test_that("qxkcd is the inverse of pxkcd up to some numerical error", {
    expect_equal(pxkcd(qxkcd(0.5)), 0.5, tolerance = 1e-6)
    expect_equal(pxkcd(qxkcd(0)), 0, tolerance = 1e-6)
    expect_equal(pxkcd(qxkcd(1)), 1, tolerance = 1e-6)
})

test_that("qxkcd deals with log.p and swap arguments as expected", {
    testp = c(0, 1, 2, 3, 4)/3
    testlogp = log(testp)
    expect_equal(qxkcd(testp, nsd = 2),
                 qxkcd(log(testp), nsd = 2, log.p = TRUE))
    expect_equal(maxy(2) - qxkcd(testp, nsd = 2, swap.end.points = TRUE),
                 qxkcd(1-testp, nsd = 2))
})