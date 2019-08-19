\name{spline.E}
\Rdversion{1.1}
\alias{spline.E}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Spline estimation function
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
spline.E(knots, degree, smooth, intercept = 0, signif = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
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
function (knots, degree, smooth, intercept = 0, signif = 3) 
{
    cmat = Cmat(knots, degree, smooth, intercept, signif)
    emat = Emat(knots, degree, smooth, !is.null(intercept), signif)
    nc = nrow(cmat)
    ne = nrow(emat)
    basisT = t(basis(cbind(t(cmat), t(emat))))
    cols = attr(basisT, "cols")
    ncc = sum(cols <= nc)
    Tmat = solve(basisT)
    Tmat[, (ncc + 1):ncol(Tmat)]
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
