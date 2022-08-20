#' The xkcd Distribution
#'
#' Density, distribution function, quantile function and random number 
#' generation for the xkcd distribution with the underlying normal standard 
#' deviation equal to \code{nsd}.
#'
#' The default value of some arguments can be seen in the Usage section
#' above.
#'
#' The xkcd distribution is derived from a normal distribution
#' \eqn{X ~ N(0, \sigma^2)} (with density \eqn{f(x)}). Given \eqn{X=x}, we
#' generate \eqn{Y|X=x ~ Uniform(0, f(x))} accordingly. The resulting marginal
#' distribution of \eqn{Y} is called xkcd distribution and we denote it as
#' \eqn{Y ~ xkcd(\sigma)}.
#'
#' \eqn{Y} is in the range \eqn{(0, 1/[\sqrt{2\pi}\sigma])}, and it has 
#' distribution function
#' \deqn{2P[X\le -h(y)]+2yh(y), Y is in (0, 1/[\sqrt{2\pi}\sigma])}
#' and density function
#' \deqn{2h(y), Y is in (0, 1/[\sqrt{2\pi}\sigma])}
#' where \eqn{h(y)} is the nonnegative root of \eqn{f(x)=y}.
#'
#' @param y Vector of quantiles.
#' @param p Vector of probabilities.
#' @param n Number of observations. If \code{length(n) > 1}, the length is taken
#' to be the nubmer required and a warning message will be produced.
#' @param nsd Vector of standard deviations for underlying normal distribution.
#' @param swap.end.points Logical; if TRUE, \code{rxkcd, dxkcd, pxkcd, qxkcd}
#' are for distribution of \eqn{1/[\sqrt{2\pi}\sigma] - Y} rather than
#'  \eqn{Y} itself.
#' @param log,log.p Logical; if TRUE, \code{dxkcd} and \code{pxkcd} will give 
#' densities and probabilities in log scale; \code{qxkcd} will regard its first
#' argument as log probabilities.
#' @param flag Indicator of using R or C for performing core functions. It can
#' only take value from \code{c("R", "C")}.
#' @param tol The desired accuracy for \code{qxkcd} when numerically inversing
#' \code{pxkcd}
#'
#' @return   \code{dxkcd} gives the density,
#' \code{pxkcd} gives the distribution function,
#' \code{qxkcd} gives the quantile function, and
#' \code{rxkcd} generates random numbers for xkcd distribution.
#'
#' Logical arguments \code{log}, \code{log.p} and \code{swap.end.points}
#' can't accept a logical vector of length greater than 1. Other than that, all
#' the other arugments (first argument of each function and \code{nsd}) will
#' broadcast according to each other.
#'
#' \code{nsd = 0} will produce error since it makes no sense to generate a
#' nonnegative number that is uniformly distributed on \eqn{(0, Inf)}.
#'
#' \code{qxkcd} is calculated by numerically inversing \code{pxkcd} via
#' \code{stats::uniroot}. The calculation accuracy can be set up via \code{tol}
#' argument. The default accuracy is \code{.Machine$double.eps}.
#'
#' @source \code{rxkcd} is based on \href{https://xkcd.com/2118/}{xkcd comic 2118},
#' and a more detailed version (including pdf and cdf) is given by Charles Geyer
#' in \href{http://www.stat.umn.edu/geyer/8054/hw/}{Doing Assignment 2},
#' though the implementation here is a little bit different (simplified some
#' expressions).
#'
#' @references xkcd comic 2118, \emph{How To Annoy A Statistician},
#' \url{https://xkcd.com/2118/}
#'
#' Charles Geyer, \emph{Doing Assignment 2},
#' \url{http://www.stat.umn.edu/geyer/8054/hw/}
#'
#' @seealso \href{https://stat.ethz.ch/R-manual/R-devel/library/stats/html/Distributions.html}{Distributions} for other standard distributions.
#'
#' @examples
#' # some tests, see formal tests in unit test
#' testy <- dnorm(0)/2
#' dnorm(dxkcd(testy)/2) == testy
#' dxkcd(0) == Inf
#' dxkcd(dnorm(0)) == 0
#'
#' # issue of numerically inversing pxkcd()
#' pxkcd(qxkcd(0.5)) == 0.5
#'
#' # check the density
#' Y <- rxkcd(1000)
#' W <- rxkcd(1000, swap.end.points = TRUE)
#' dxkcd_w <- function(w) dxkcd(w, swap.end.points = TRUE)
#' par(mfrow = c(1, 2))
#' hist(Y, freq = FALSE, breaks = 20)
#' curve(dxkcd, 0, dnorm(0), add = TRUE, lwd = 2)
#' hist(W, freq = FALSE, breaks = 20)
#' curve(dxkcd_w, 0, dnorm(0), add = TRUE, lwd = 2)
#' par(mfrow = c(1, 1))
#'
#' @importFrom stats dnorm rnorm runif uniroot
#'
#'
#' @name xkcd
#' @aliases rxkcd dxkcd pxkcd qxkcd
#' @useDynLib xkcd
NULL

#' @rdname xkcd
#' @export

# random number generator for xkcd distribution
rxkcd <- function(n, nsd = 1, swap.end.points = FALSE, flag = "R") {
    stopifnot(is.numeric(n))
    stopifnot(is.finite(n))
    if (any(n != round(n))) stop("n must be integer")
    if (any(n <= 0)) stop("n must be positive")
    if (length(n) > 1){
        warning("n is a vector: take its length as the first argument of rxkcd")
        n = length(n)
    }

    stopifnot(is.numeric(nsd))
    stopifnot(is.finite(nsd))
    if (nsd <= 0)
      stop("variance for the normal distribution must be positive")

    stopifnot(is.logical(swap.end.points))
    if (length(swap.end.points) > 1)
      stop("swap.end.points shouldn't be a vector")

    stopifnot(is.character(flag))
    if (!(flag %in% c("R", "C")))
      stop("flag can only be either \"R\" or \"C\"")

    if (flag == "R") {
      x <- rnorm(n, sd = nsd)
      y <- runif(n, max = dnorm(x, sd = nsd))
      if (swap.end.points) {
         return(maxy(nsd) - y)
      } else {
         return(y)
      }
    }
    if (flag == "C") {
       out <- .C("crxkcd", n = as.integer(n), nsd = as.double(nsd),
       swap.end.points = as.logical(swap.end.points), x = double(n))
       return(out$x)
    }
}

#' @rdname xkcd
#' @export
dxkcd <- function(y, nsd = 1, log = FALSE, swap.end.points = FALSE,
                  flag = "R") {
    stopifnot(is.numeric(y))

    stopifnot(is.numeric(nsd))
    stopifnot(is.finite(nsd))
    if (nsd <= 0) 
      stop("variance for the normal distribution must be positive")

    stopifnot(is.logical(log))
    if (length(log) > 1) stop("log shouldn't be a vector")
    stopifnot(is.logical(swap.end.points))
    if (length(swap.end.points) > 1)
      stop("swap.end.points shouldn't be a vector")

    stopifnot(is.character(flag))
    if (!(flag %in% c("R", "C")))
      stop("flag can only be either \"R\" or \"C\"")

    if (y < 0 || y > maxy(nsd))
      return(ifelse(log, log(0), 0))

    if (flag == "R") {
      if (swap.end.points == FALSE) {
         # out-of-domain cases
         # y \in [0, maxy(nsd)]
         y <- ifelse(swap.end.points, maxy(nsd) - y, y)
         x <- invpdf(y, sd = nsd)
         return(ifelse(log, log(2 * x), 2 * x))
      } else {
         w <- maxy(nsd) - y
         dw <- dxkcd(w, nsd)
         return(ifelse(log, log(dw), dw))
      }
    }

    if (flag == "C"){
       out <- .C("cdxkcd", y = as.double(y),
                 nsd = as.double(nsd), log = as.logical(log),
                 swap.end.points = as.logical(swap.end.points), d = double(1))
      return(out$d)
    }
}
dxkcd <- Vectorize(dxkcd, vectorize.args = c("y", "nsd"))

#' @rdname xkcd
#' @export
# cdf of xkcd distribution
pxkcd <- function(y, nsd = 1, log.p = FALSE, swap.end.points = FALSE,
                  flag = "R") {
    stopifnot(is.numeric(y))

    stopifnot(is.numeric(nsd))
    stopifnot(is.finite(nsd))
    if (nsd <= 0)
      stop("variance for the normal distribution must be positive")

    stopifnot(is.logical(log.p))
    if (length(log.p) > 1) stop("log.p shouldn't be a vector")
    stopifnot(is.logical(swap.end.points))
    if (length(swap.end.points) > 1)
      stop("swap.end.points shouldn't be a vector")

    stopifnot(is.character(flag))
    if (!(flag %in% c("R", "C")))
      stop("flag can only be either \"R\" or \"C\"")

    if (y <= 0)
      return(ifelse(log.p, log(0), 0))
    if (y >= maxy(nsd))
      return(ifelse(log.p, log(1), 1))

    if (flag == "R") {
      if (swap.end.points == FALSE) {
         if (y <= 0) {
            Fy <- 0
         } else if (y >= maxy(nsd)) {
            Fy <- 1
         } else {
            # y \in (0, maxy(nsd))
            x <- invpdf(y, sd = nsd) # f(x)=y
            Fy <- 2 * pnorm(-x, sd = nsd) + 2 * y * x
         }
         return(ifelse(log.p, log(Fy), Fy))
      } else {
         w <- maxy(nsd) - y
         Fw <- 1 - pxkcd(w, nsd)
         return(ifelse(log.p, log(Fw), Fw))
      }
    }

    if (flag == "C") {
       out <- .C("cpxkcd", y = as.double(y),
                 nsd = as.double(nsd), log.p = as.logical(log.p),
                 swap.end.points = as.logical(swap.end.points), p = double(1))
      return(out$p)
    }
}
pxkcd <- Vectorize(pxkcd, vectorize.args = c("y", "nsd"))

#' @rdname xkcd
#' @export
# quantile of xkcd distribution
# numerically inverse pxkcd function via uniroot()
qxkcd <- function(p, nsd = 1, log.p = FALSE, swap.end.points = FALSE,
                  tol = .Machine$double.eps, flag = "R") {
    stopifnot(is.logical(log.p))
    if (length(log.p) > 1) stop("log.p shouldn't be a vector")
    stopifnot(is.logical(swap.end.points))
    if (length(swap.end.points) > 1)
      stop("swap.end.points shouldn't be a vector")
    if (log.p) p <- exp(p)

    stopifnot(is.numeric(p))
    stopifnot(is.finite(p))
    if (p<0 || p>1) {
       warning("p (or exp(p)) not in [0, 1]: NaNs produced")
       return(NaN)
    }

    if (p == 0) return(0)
    if (p == 1) return(maxy(nsd))

    stopifnot(is.numeric(nsd))
    stopifnot(is.finite(nsd))
    if (nsd <= 0)
      stop("variance for the normal distribution must be positive")

    stopifnot(is.character(flag))
    if (!(flag %in% c("R", "C")))
      stop("flag can only be either \"R\" or \"C\"")

    # inverse of cdf of xkcd distribution with nsd and default setup
    if (flag == "R"){
      if (swap.end.points) {
         pxkcd.fun <- function(y) pxkcd(y, nsd = nsd) - (1 - p)
         pxkcd_inv <- uniroot(pxkcd.fun, lower = 0, upper = maxy(nsd), tol = tol)$root
         return(maxy(nsd) - pxkcd_inv)
      } else {
         pxkcd.fun <- function(y) pxkcd(y, nsd = nsd) - p
         pxkcd_inv <- uniroot(pxkcd.fun, lower = 0, upper = maxy(nsd), tol = tol)$root
         return(pxkcd_inv)
      }
    }

    if (flag == "C") {
       out <- .C("cqxkcd", p = as.double(p),
                 nsd = as.double(nsd), log.p = as.logical(log.p),
                 swap.end.points = as.logical(swap.end.points), 
                 tol = as.double(tol), q = double(1))
      return(out$q)
    }
}
qxkcd <- Vectorize(qxkcd, vectorize.args = c("p", "nsd"))
