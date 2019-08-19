\name{Lag}
\Rdversion{1.1}
\concept{longitudinal}
\concept{hierarchical}
\alias{cLag}
\alias{Lag}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
'Contextual' Lag with respect to time within id
}
\description{
Lag a vector with respect to time order
}
\usage{
cLag(x, id = rep(1, length(x)), time = 1:length(x), lag = 1, at = NULL, check = TRUE, idx)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
vector of values to be lagged
}
  \item{id}{
identifies clusters within which lagging takes place
}
  \item{time}{
time variable for lagging, generally integer valued
}
  \item{lag}{
period for lagging: lag = 1 reaches one unit of \code{time} into the past. Negative
numbers reach into the future.
}
  \item{at}{
uses the value of \code{x} at a particular value of \code{time}. e.g. \code{at = 1}
returns the value of \code{x} corresponding to \code{time == 1}.
}
  \item{check}{
that \code{id/time} combinations are unique.
}
  \item{idx}{
for compatibility with a previous version.
}
}
\details{
The function can also be called as \code{Lag} which is now deprecated
to avoid conflicts with \code{Hmisc::Lag}. Use \code{cLagI} for intra(extra)polative lagging.
}
\value{
A vector of values of \code{x} lagged within levels of \code{id}. The value returned at \code{time == t} is the value
of \code{x} at \code{time == t-1}.  If there is no observation for which \code{time == t-1} then the returned value is
\code{NA}.  To lag to the previous value of \code{time}, one can use \code{rank}. Consider, also, \link{\code{cLagI}} that intrapolates backwards one
unit of \code{time}.
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
\code{\link{cLagI}},
\code{\link{cDiff}},
\code{\link{cDiffI}},
\code{\link{capply}},
\code{\link{up}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, id = rep(1, length(x)), time = 1:length(x), lag = 1, 
    at = NULL, check = T, idx) 
{
    if (check) {
        comb <- paste(id, time)
        if (any(duplicated(comb))) 
            stop("id not unique in each level of time")
    }
    if (!missing(idx)) 
        time = idx
    ret <- x
    names(ret) <- paste(id, time, sep = "::")
    retnn <- if (is.null(at)) 
        paste(id, time - lag, sep = "::")
    else paste(id, at, sep = "::")
    ret[retnn]
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
