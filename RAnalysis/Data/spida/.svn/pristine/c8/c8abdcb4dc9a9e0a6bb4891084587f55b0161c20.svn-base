\name{atotal}
\Rdversion{1.1}
\alias{atotal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Add all marginal sums to an array
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
atotal(arr, FUN = sum, label = "Total", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{arr}{
%%     ~~Describe \code{arr} here~~
}
  \item{FUN}{
%%     ~~Describe \code{FUN} here~~
}
  \item{label}{
%%     ~~Describe \code{label} here~~
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
\code{\link{tab}},
\code{\link{acond}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (arr, FUN = sum, label = "Total", ...) 
{
    help <- "\natotal                coursefun.R     for PSYC 6140/MATH 6630 05/06\n\nAdds border of sums to an array\n\nDescription:\n\n     'atotal' adds by default a border of sums to an array.\n     The function FUN may be used instead of 'sum'. Additional\n     arguments to FUN can also be given.\n\nUsage:\n\n     atotal( arr , FUN = sum, label = 'Total', ...)\n\nArguments:\n\n     arr: array, matrix or vector\n\n     FUN: function to be applied to cross sections of arr\n\n     ...: additional arguments to FUN\n\nDetails:\n\nValue:\n\n     An array with dimension dim(arr) + 1\n\nReferences:\n\nContributed by:  G. Monette  2005-10-10\n\nModifications:\n     2007-12-17: Fixed bug so dimnames is preserved for one-dimensional tables\n\n"
    d <- dim(arr)
    cls <- class(arr)
    dim1 <- FALSE
    if (length(d) == 1) {
        dim1 <- TRUE
        dn <- dimnames(arr)
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
        ret <- structure(c(arr, FUN(arr, ...)), names = c(names(arr), 
            label), class = cls)
        if (dim1) {
            dn[[1]] <- c(dn[[1]], label)
            ret <- structure(ret, dim = length(ret), dimnames = dn)
        }
        return(ret)
    }
    n <- length(d)
    ret <- arr
    ind <- 1:n
    for (i in n:1) {
        new <- apply(ret, ind[-i], FUN, ...)
        ret <- abind(ret, new, i, label)
    }
    class(ret) <- cls
    ret
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
