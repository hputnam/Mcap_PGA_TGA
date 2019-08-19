\name{print.data.frame.lab}
\Rdversion{1.1}
\alias{print.data.frame.lab}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Print a labeled data frame
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
print.data.frame.lab(x, ..., digits = NULL, quote = FALSE, right = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
  \item{digits}{
%%     ~~Describe \code{digits} here~~
}
  \item{quote}{
%%     ~~Describe \code{quote} here~~
}
  \item{right}{
%%     ~~Describe \code{right} here~~
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
function (x, ..., digits = NULL, quote = FALSE, right = TRUE) 
{
    labs <- attributes(x)$labs
    if (length(x) == 0) {
        cat("NULL data frame with", length(row.names(x)), "rows\n")
    }
    else if (length(row.names(x)) == 0) {
        print.default(names(x), quote = FALSE)
        cat("<0 rows> (or 0-length row.names)\n")
    }
    else {
        mat <- as.matrix(format.data.frame(x, digits = digits, 
            na.encode = FALSE))
        labs <- c(labs, "", "")
        labs <- labs[1:2]
        names(dimnames(mat)) <- labs
        print(mat, ..., quote = quote, right = right)
    }
    invisible(x)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
