\name{gsp.v01}
\Rdversion{1.1}
\alias{gsp.v01}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Gsp -- older version
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gsp.v01(x, k, degree = rep(3, length(k) + 1), smooth = rep(2, length(k)), intercept = FALSE, debug = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{k}{
%%     ~~Describe \code{k} here~~
}
  \item{degree}{
%%     ~~Describe \code{degree} here~~
}
  \item{smooth}{
%%     ~~Describe \code{smooth} here~~
}
  \item{intercept}{
%%     ~~Describe \code{intercept} here~~
}
  \item{debug}{
%%     ~~Describe \code{debug} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, k, degree = rep(3, length(k) + 1), smooth = rep(2, 
    length(k)), intercept = FALSE, debug = FALSE) 
{
    if (!debug) 
        disp = function(x) invisible(0)
    kernel = function(L) {
        QR = qr(t(L))
        qr.Q(QR, complete = T)[, -seq_len(QR$rank)]
    }
    Xfull = function(x, k, degree) {
        if (length(degree) != length(k) + 1) 
            stop("length( degree ) must = length( k ) + 1")
        Xmat = function(x, degree) {
            do.call("cbind", lapply(0:degree, function(i) if (i == 
                0) 
                rep(1, length(x))
            else x^i))
        }
        k = sort(k)
        g = cut(x, c(-Inf, k, Inf))
        Xraw = Xmat(x, max(degree))
        do.call("cbind", lapply(seq_along(degree), function(iint) (g == 
            levels(g)[iint]) * Xraw[, 1:(degree[iint] + 1), drop = F]))
    }
    X0 = Xfull(0, k, degree)
    if (FALSE) {
        nk = length(k)
        extend = (max(k) - min(k) + 1)/nk
        ints = c(min(k) - extend, sort(k), max(k) + extend)
        xs = approx(ints, n = (max(degree) + 1) * (nk + 1))$y
        xs = c(-xs, xs)
        X.full = Xfull(xs, k, degree)
        disp(X.full)
    }
    fs = list(f = function(x) {
        c(1, x, x^2, x^3, x^4)
    }, f1 = function(x) {
        c(0, 1, 2 * x, 3 * x^2, 4 * x^3)
    }, f2 = function(x) {
        c(0, 0, 2, 6 * x, 12 * x^2)
    }, f3 = function(x) {
        c(0, 0, 0, 6, 24 * x)
    }, f4 = function(x) {
        c(0, 0, 0, 0, 24)
    })
    nints = length(k) + 1
    endcols = cumsum(degree + 1)
    disp(endcols)
    startcols = c(1, endcols[-length(endcols)] + 1)
    disp(startcols)
    ncols = degree + 1
    disp(ncols)
    Lsmooth = matrix(0, ncol = ncol(X0), nrow = sum(smooth + 
        1))
    irow = 0
    for (i in seq_along(k)) {
        disp(i)
        knot = k[i]
        iplus = startcols[i]:endcols[i]
        disp(iplus)
        iminus = startcols[i + 1]:endcols[i + 1]
        disp(iminus)
        for (j in 0:smooth[i]) {
            irow = irow + 1
            Lsmooth[irow, iplus] = (fs[[j + 1]](knot)[1:ncols[i]])
            Lsmooth[irow, iminus] = (-fs[[j + 1]](knot)[1:ncols[i + 
                1]])
        }
    }
    disp(Lsmooth)
    disp(Lsmooth \%*\% kernel(Lsmooth))
    if (intercept == FALSE) 
        Lsmooth = rbind(Lsmooth, X0)
    H = kernel(Lsmooth)
    Xret = Xfull(x, k, degree) \%*\% H
    Xret
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
