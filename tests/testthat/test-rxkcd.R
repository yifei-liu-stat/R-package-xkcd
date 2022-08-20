test_that("rxkcd generates random numbers of expected length and gives suitable warning when n is a vector", {
    expect_equal(length(rxkcd(5)), 5)
    expect_equal(length(rxkcd(1:6)), 6)
    expect_warning(length(rxkcd(1:6)))
})

test_that("rxkcd generates random numbers of suitable range", {
    expect_true(max(rxkcd(100, 3)) <= maxy(3))
    expect_true(min(rxkcd(100, 0.5)) >= 0)
})

test_that("rxkcd 's swap.end.points argument works as expected", {
    set.seed(1234)
    x <- rxkcd(1, 2)
    set.seed(1234)
    y <- rxkcd(1, 2, swap.end.points = TRUE)
    expect_equal(x, maxy(2) - y)
}
)