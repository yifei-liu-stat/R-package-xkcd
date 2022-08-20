test_that("pxkcd gives proper probability for in-domain y", {
    testy <- maxy(1)/2
    expect_true(pxkcd(testy)>0 && pxkcd(testy)<1)
    expect_equal(pxkcd(0), 0)
    expect_equal(pxkcd(maxy(1)), 1)
})

test_that("pxkcd gives proper probability for out-of-domain y", {
    expect_equal(pxkcd(-100), 0)
    expect_equal(pxkcd(100, 1), 1)
})

test_that("pxkcd broadcasts argument dimensions in the expected way", {
    expect_equal(length(pxkcd(1:10)), 10)
    expect_equal(length(pxkcd(1:2, 1:10)), 10)
    expect_equal(length(pxkcd(1:10, 1:10)), 10)
})

test_that("pxkcd can be obtained by integraing dxkcd", {
    expect_equal(integrate(dxkcd, 0.1, 0.2)$value, pxkcd(0.2)-pxkcd(0.1))
    expect_equal(integrate(dxkcd, 0.07, 0.13)$value, pxkcd(0.13)-pxkcd(0.07))
})

test_that("pxkcd deals with log.p and swap arguments as expected", {
    testy <- c(-Inf, -1, 0, 1, 2, 3, 4, Inf) / 3 * maxy(2)
    expect_equal(log(pxkcd(testy, nsd = 2)), pxkcd(testy, nsd = 2, log.p = TRUE))
    expect_equal(1 - pxkcd(maxy(2) - testy, nsd = 2), pxkcd(testy, nsd = 2, swap.end.points = TRUE))
    expect_equal(log(1 - pxkcd(maxy(2) - testy, nsd = 2)), pxkcd(testy, nsd = 2, log.p = TRUE, swap.end.points = TRUE))
})