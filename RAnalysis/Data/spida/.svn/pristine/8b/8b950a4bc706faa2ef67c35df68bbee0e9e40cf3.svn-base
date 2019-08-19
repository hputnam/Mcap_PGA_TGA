\name{Vcov.rdc}
\Rdversion{1.1}
\alias{Vcov.rdc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Version of Vcov used in RDC
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Vcov.rdc(fit, L = Lmat(fit, ""))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
%%     ~~Describe \code{fit} here~~
}
  \item{L}{
%%     ~~Describe \code{L} here~~
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
function (fit, L = Lmat(fit, "")) 
{
    vc <- getFix(fit)$vcov
    L \%*\% vc \%*\% t(L)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
