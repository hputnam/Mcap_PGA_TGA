\name{apct}
\Rdversion{1.1}
\alias{apct}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Transform a frequency table to percentages
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
apct(x, MARGIN = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{MARGIN}{
%%     ~~Describe \code{MARGIN} here~~
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
function (x, MARGIN = 1) 
{
    if (length(dim(x)) == 1) {
        ret <- cbind(N = x, pct = 100 * x/sum(x, na.rm = T))
        ret <- rbind(ret, Total = apply(ret, 2, sum, na.rm = T))
        print(round(ret, 1))
        return(invisible(ret))
    }
    ret <- list(N = atotal(x), pct = 100 * acond(x, MARGIN))
    cat("\nN:\n")
    print(ret[[1]])
    cat("\nPercentage:\n")
    print(round(ret[[2]], 1))
    invisible(ret)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
