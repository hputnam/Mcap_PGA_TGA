\name{acond}
\Rdversion{1.1}
\alias{acond}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Transform a frequency table into relative frequencies relative to a margin.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
acond(x, MARGIN = NULL)
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
function (x, MARGIN = NULL) 
{
    debug <- F
    x <- atotal(x)
    if (is.null(MARGIN)) 
        return(x/x[length(x)])
    d <- dim(x)
    if (debug) {
        cat("dim(x)\n")
        print(d)
    }
    dn <- dimnames(x)
    if (debug) {
        cat("dimnames(x)\n")
        print(dn)
    }
    n <- length(d)
    m <- length(MARGIN)
    perm <- c(MARGIN, setdiff(1:n, MARGIN))
    perm.inverse <- order(perm)
    if (debug) 
        print(perm)
    if (debug) 
        print(perm.inverse)
    x <- aperm(x, perm)
    d <- dim(x)
    zl <- list(x)
    for (ii in 1:length(d)) zl[[ii + 1]] <- seq(if (ii <= m) 
        1
    else d[ii], d[ii])
    tots <- do.call("[", zl)
    ret <- array(c(x)/c(tots), dim = d)
    ret <- aperm(ret, perm.inverse)
    if (debug) 
        print(dim(ret))
    if (debug) 
        print(dn)
    dimnames(ret) <- dn
    ret
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
