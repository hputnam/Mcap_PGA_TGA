\name{Rbind.data.frame}
\Rdversion{1.1}
\alias{Rbind.data.frame}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Rbind applied to data frames
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Rbind.data.frame(..., vname = ".which", oname = ".order")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
  \item{vname}{
%%     ~~Describe \code{vname} here~~
}
  \item{oname}{
%%     ~~Describe \code{oname} here~~
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
function (..., vname = ".which", oname = ".order") 
{
    z <- list(...)
    for (ii in 1:length(z)) {
        z[[ii]] <- as.data.frame(z[[ii]], stringsAsFactors = FALSE)
        z[[ii]][, vname] <- rep(ii, nr <- nrow(z[[ii]]))
        z[[ii]][, oname] <- 1:nr
    }
    ret <- z[[1]]
    for (ii in 2:length(z)) ret <- merge(ret, z[[ii]], all = TRUE, 
        sort = FALSE)
    ret[order(ret[[vname]], ret[[oname]]), ]
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
