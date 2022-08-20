# upper bound of y: \frac{1}{\sqrt{2\pi} \sigma}
maxy <- function(nsd = 1) 1 / (sqrt(2 * pi) * nsd)

# upper inverse of pdf of N(0, \sigma^2)
invpdf <- function(y, sd = 1) sqrt(-2 * sd^2 * log(y * sqrt(2 * pi) * sd))