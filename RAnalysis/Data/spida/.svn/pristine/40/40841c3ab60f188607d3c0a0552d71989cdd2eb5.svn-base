\name{ConjComp}
\Rdversion{1.1}
\alias{ConjComp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Conjugate complement of span(X) in span(Z) with respect to inner product ip
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ConjComp(X, Z = diag(nrow(X)), ip = diag(nrow(X)), tol = 1e-07)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{Z}{
%%     ~~Describe \code{Z} here~~
}
  \item{ip}{
%%     ~~Describe \code{ip} here~~
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
function (X, Z = diag(nrow(X)), ip = diag(nrow(X)), tol = 1e-07) 
{
    help <- "\n    ConjComp returns a basis for the conjugate complement of the\n    conjugate projection of X into span(Z) with respect to inner product with\n    matrix ip.\n    Note: Z is assumed to be of full column rank but not necessarily X.\n    "
    xq <- qr(t(Z) \%*\% ip \%*\% X, tol = tol)
    if (xq$rank == 0) 
        return(Z)
    a <- qr.Q(xq, complete = T)[, -(1:xq$rank)]
    Z \%*\% a
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
