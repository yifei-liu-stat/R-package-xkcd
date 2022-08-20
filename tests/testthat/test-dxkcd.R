test_that("dxkcd gives the right density for in-domain y", {
    testy <- maxy(2)/2
    expect_equal(dnorm(dxkcd(testy, nsd = 2)/2, sd = 2), testy)
    expect_equal(dxkcd(0, 2), Inf)
    expect_equal(dxkcd(maxy(2), nsd = 2), 0)
})

test_that("dxkcd gives the right density for out-of-domain y", {
    expect_equal(dxkcd(-10), 0)
    expect_equal(dxkcd(100), 0)
})

test_that("dxkcd broadcasts argument dimensions in the expected way", {
    expect_equal(length(dxkcd(1:10)), 10)
    expect_equal(length(dxkcd(1:2, 1:10)), 10)
    expect_equal(length(dxkcd(1:10, 1:10)), 10)
})

test_that("dxkcd deals with log and swap arguments as expepcted", {
    testy <- c(-Inf, -1, 0, 1, 2, 3, 4, Inf) / 3 * maxy(2)
    expect_equal(log(dxkcd(testy, nsd = 2)), dxkcd(testy, nsd = 2, log = TRUE))
    expect_equal(log(dxkcd(maxy(2), nsd = 2)), dxkcd(maxy(2), nsd = 2, log = TRUE))
    expect_equal(dxkcd(maxy(2)-testy, nsd = 2), dxkcd(testy, nsd = 2, swap.end.points = TRUE))
    expect_equal(log(dxkcd(maxy(2)-testy, nsd = 2)), dxkcd(testy, nsd = 2, log = TRUE, swap.end.points = TRUE))
})