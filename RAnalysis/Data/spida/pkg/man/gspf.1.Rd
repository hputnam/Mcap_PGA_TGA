\name{gspf.1}
\Rdversion{1.1}
\alias{gspf.1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Marked for deletion
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gspf.1(x, knots, degree = 3, smooth = pmax(pmin(degree[-1], degree[-length(degree)]) - 1, 0), intercept = 0, signif = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{knots}{
%%     ~~Describe \code{knots} here~~
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
  \item{signif}{
%%     ~~Describe \code{signif} here~~
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
function (x, knots, degree = 3, smooth = pmax(pmin(degree[-1], 
    degree[-length(degree)]) - 1, 0), intercept = 0, signif = 3) 
{
    tmat = spline.T(knots, degree, smooth, intercept = intercept, 
        signif = signif)
    nC = attr(tmat, "nC")
    nE = ncol(tmat) - nC
    tmat = tmat[, c(seq(nC + 1, ncol(tmat)), 2:nC, 1)]
    X = Xf(x, knots, max(degree), signif = signif) \%*\% tmat
    qqr = qr(X)
    Q = qr.Q(qqr)
    R = qr.R(qqr)
    gsp(seq(-2, 10), 0, c(2, 2), 0)
    ret
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
