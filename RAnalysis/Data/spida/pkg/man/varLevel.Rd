\name{varLevel}
\Rdversion{1.1}
\alias{varLevel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identify level of aggregation at which a variable in invariant
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
varLevel(x, form, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{form}{
%%     ~~Describe \code{form} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (x, form, ...) 
{
    sel <- model.frame(form, x)
    z <- list()
    idx <- rep("", nrow(x))
    z[[1]] <- constant(x, ...)
    for (ii in 1:ncol(sel)) {
        idx <- paste(idx, as.character(sel[[ii]]), sep = ";")
        z[[ii + 1]] <- constant(x, idx, all = T, ...)
    }
    ret <- do.call("rbind", z)
    ret <- length(z) - apply(ret * 1, 2, sum)
    ret
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
