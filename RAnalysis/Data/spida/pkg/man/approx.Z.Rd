\name{approx.Z}
\Rdversion{1.1}
\alias{approx.Z}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Approximation for semi-parametric splines
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
approx.Z(Z, oldtimes, newtimes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z}{
%%     ~~Describe \code{Z} here~~
}
  \item{oldtimes}{
%%     ~~Describe \code{oldtimes} here~~
}
  \item{newtimes}{
%%     ~~Describe \code{newtimes} here~~
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
function (Z, oldtimes, newtimes) 
{
    oldt.u <- sort(unique(oldtimes))
    if (length(oldt.u) != length(oldtimes) || any(oldt.u != oldtimes)) {
        Z <- Z[match(oldt.u, oldtimes), ]
        oldtimes <- oldt.u
    }
    apply(Z, 2, function(u, oldt, newt) {
        approx(oldt, u, xout = newt)$y
    }, oldt = oldtimes, newt = newtimes)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
