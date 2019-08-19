\name{Proj.test}
\Rdversion{1.1}
\alias{Proj.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Test accuracy of projection method.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Proj.test(x, fun = Proj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{fun}{
%%     ~~Describe \code{fun} here~~
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
function (x, fun = Proj) 
{
    help = "\n     # sample test:\n     zz <- matrix( rnorm(120), ncol = 3) \%*\% diag(c(10000,1,.000001))\n     Proj.test( zz, fun = Proj)\n     Proj.test( zz, fun = Proj.1)\n"
    pp <- fun(x)
    ret <- pp \%*\% pp - pp
    attr(ret, "maxabs") <- max(abs(ret))
    attr(ret, "maxabs")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
