\name{Lag.0}
\Rdversion{1.1}
\alias{Lag.0}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Lag within subject: older less efficient version of Lag
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Lag.0(x, id, idx, lag = 1, at = NULL, check = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{id}{
%%     ~~Describe \code{id} here~~
}
  \item{idx}{
%%     ~~Describe \code{idx} here~~
}
  \item{lag}{
%%     ~~Describe \code{lag} here~~
}
  \item{at}{
%%     ~~Describe \code{at} here~~
}
  \item{check}{
%%     ~~Describe \code{check} here~~
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
function (x, id, idx, lag = 1, at = NULL, check = T) 
{
    if (check) {
        comb <- paste(id, idx)
        if (any(duplicated(comb))) 
            stop("id not unique in each level of idx")
    }
    if (any(is.na(idx))) 
        stop("NAs in idx")
    if (any(round(idx) != idx)) 
        stop("Non integers in idx")
    ret <- x
    id <- as.character(id)
    names(x) <- id
    for (i in max(idx):min(idx)) {
        to.pos <- idx == i
        if (is.null(at)) 
            from <- x[idx == (i - lag)]
        else from <- x[idx == at]
        ids <- names(x[to.pos])
        ret[to.pos] <- from[ids]
    }
    ret
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
