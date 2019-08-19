\name{atotal.rdc}
\Rdversion{1.1}
\alias{atotal.rdc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Version of atotal in RDC
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
atotal.rdc(arr, FUN = sum, name = "Total", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{arr}{
%%     ~~Describe \code{arr} here~~
}
  \item{FUN}{
%%     ~~Describe \code{FUN} here~~
}
  \item{name}{
%%     ~~Describe \code{name} here~~
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
function (arr, FUN = sum, name = "Total", ...) 
{
    d <- dim(arr)
    if (length(d) == 1) {
        arr <- c(arr)
        d <- dim(arr)
    }
    if (is.character(FUN)) 
        FUN <- get(FUN, mode = "function")
    else if (mode(FUN) != "function") {
        farg <- substitute(FUN)
        if (mode(farg) == "name") 
            FUN <- get(farg, mode = "function")
        else stop(paste("\"", farg, "\" is not a function", sep = ""))
    }
    if (is.null(d)) {
        ret <- c(arr, FUN(arr))
        names(ret)[length(ret)] = name
        return(ret)
    }
    n <- length(d)
    name <- rep(name, length = n)
    ret <- arr
    ind <- 1:n
    for (i in n:1) {
        new <- apply(ret, ind[-i], FUN, ...)
        ret <- abind(ret, new, i, name[i])
    }
    ret
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
