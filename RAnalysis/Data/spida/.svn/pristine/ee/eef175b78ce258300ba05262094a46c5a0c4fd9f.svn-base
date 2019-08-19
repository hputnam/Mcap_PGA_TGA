\name{oplot}
\Rdversion{1.1}
\alias{oplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Indicate number of points overplotted
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
oplot(x, y, ..., verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
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
function (x, y, ..., verbose = T) 
{
    pat <- paste(x, y, sep = ",")
    keep <- !duplicated(pat)
    ns <- table(pat)
    ns <- ns[pat[keep]]
    nps <- as.character(ns)
    x <- x[keep]
    y <- y[keep]
    if (verbose) {
        print(pat)
        print(table(pat))
        print(keep)
        print(ns)
    }
    plot(x, y, pch = "o", cex = 5, ...)
    text(x, y, nps)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
