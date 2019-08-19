\name{ginv.1}
\Rdversion{1.1}
\alias{ginv.1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generalized inverse -- test version
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ginv.1(x, tol = .Machine$double.eps * d[1] * max(nrow(x), ncol(x)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{tol}{
%%     ~~Describe \code{tol} here~~
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
function (x, tol = .Machine$double.eps * d[1] * max(nrow(x), 
    ncol(x))) 
{
    r = svd(x)
    d = r$d
    dinv = ifelse(d < tol, 0, 1/d)
    ret = r$v \%*\% (dinv * t(r$u))
    dimnames(ret) = dimnames(x)[2:1]
    ret
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
