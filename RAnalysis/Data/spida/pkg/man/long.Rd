\name{long}
\Rdversion{1.1}
\alias{long}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a longitudinal file from a wide file
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
long(data, varying = NULL, sep.varying = "\\.", v.names = names(Varying), timevar = "time", idvar = "id", ids = 1:NROW(data), times = seq(length = length(Varying[[1]])), drop = NULL, new.row.names = NULL, split = list(regexp = "\\.", include = FALSE), debug = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{varying}{
%%     ~~Describe \code{varying} here~~
}
  \item{sep.varying}{
%%     ~~Describe \code{sep.varying} here~~
}
  \item{v.names}{
%%     ~~Describe \code{v.names} here~~
}
  \item{timevar}{
%%     ~~Describe \code{timevar} here~~
}
  \item{idvar}{
%%     ~~Describe \code{idvar} here~~
}
  \item{ids}{
%%     ~~Describe \code{ids} here~~
}
  \item{times}{
%%     ~~Describe \code{times} here~~
}
  \item{drop}{
%%     ~~Describe \code{drop} here~~
}
  \item{new.row.names}{
%%     ~~Describe \code{new.row.names} here~~
}
  \item{split}{
%%     ~~Describe \code{split} here~~
}
  \item{debug}{
%%     ~~Describe \code{debug} here~~
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
function (data, varying = NULL, sep.varying = "\\.", v.names = names(Varying), 
    timevar = "time", idvar = "id", ids = 1:NROW(data), times = seq(length = length(Varying[[1]])), 
    drop = NULL, new.row.names = NULL, split = list(regexp = "\\.", 
        include = FALSE), debug = FALSE) 
{
    help <- "\n    Example:   long(dwl, varying = c('wmrss', 'lm1tot','lm2tot'), sep='')\n    "
    nn <- names(data)
    if (debug) 
        disp(nn)
    Varying <- list()
    for (ii in 1:length(varying)) {
        if (length(prefix <- varying[[ii]]) == 1) {
            ns <- grep(paste("^", prefix, sep.varying, sep = ""), 
                nn, value = T)
            if (debug) 
                disp(ns)
            Varying[[ii]] <- ns
            names(Varying)[ii] <- prefix
        }
        else {
            Varying[[ii]] <- varying[[ii]]
            names(Varying)[ii] <- names(varying)[ii]
        }
    }
    if (debug) 
        disp(Varying)
    if (debug) 
        disp(times)
    ret <- stats::reshape(data, Varying, v.names, timevar, idvar, 
        ids, times, drop, direction = "long", new.row.names, 
        split)
    ret[order(ret[[idvar]], ret[[timevar]]), ]
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
