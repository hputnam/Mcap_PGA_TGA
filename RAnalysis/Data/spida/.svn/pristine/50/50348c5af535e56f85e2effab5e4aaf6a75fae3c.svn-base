\name{formatCoefmat}
\Rdversion{1.1}
\alias{formatCoefmat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Format a coefficient table
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
formatCoefmat(x, digits = 6, pdigits = digits - 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{digits}{
%%     ~~Describe \code{digits} here~~
}
  \item{pdigits}{
%%     ~~Describe \code{pdigits} here~~
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
function (x, digits = 6, pdigits = digits - 1) 
{
    pformat <- function(x, digits) {
        x <- format(xx <- round(x, digits))
        x[as.double(xx) == 0] <- paste(c("<.", rep("0", digits - 
            1), "1"), collapse = "")
        x
    }
    xx <- array("", dim = dim(x), dimnames = dimnames(x))
    for (i in 1:ncol(xx)) {
        xx[, i] <- format(round(x[, i], digits), digits = digits)
    }
    if (length(isp <- grep("^Pr\\(", colnames(x))) > 0) {
        xx[, isp[1]] <- pformat(x[, isp[1]], digits = pdigits)
    }
    xx
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
