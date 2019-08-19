\name{print.glh.rdc}
\Rdversion{1.1}
\alias{print.glh.rdc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Print a 'glh' object tested in RDC
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
print.glh.rdc(x, round = 6, pround = round - 1, L = T, cov = T, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{round}{
%%     ~~Describe \code{round} here~~
}
  \item{pround}{
%%     ~~Describe \code{pround} here~~
}
  \item{L}{
%%     ~~Describe \code{L} here~~
}
  \item{cov}{
%%     ~~Describe \code{cov} here~~
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
function (x, round = 6, pround = round - 1, L = T, cov = T, ...) 
{
    rnd <- function(x, digits) {
        if (is.numeric(x)) 
            x <- round(x, digits = digits)
        format(x)
    }
    for (ii in 1:length(x)) {
        nn <- names(x)[[ii]]
        tt <- x[[ii]]
        ta <- tt$anova
        tap <- array("", dim = dim(ta), dimnames = dimnames(ta))
        cat("\n", nn, "\n", sep = "")
        print(formatCoefmat(ta, digits = round, pdigits = pround), 
            quote = F, right = T)
        cat("\n")
        te <- tt$estimate
        if (!is.null(zhead <- attr(tt, "heading"))) 
            cat(zhead, "\n")
        print(formatCoefmat(te, digits = round, pdigits = pround), 
            quote = F, right = T)
        if (L == T) {
            cat("\nL:\n")
            print(tt$L)
            if (dim(tt$L.full)[1] < dim(tt$L)[1]) {
                cat("\nL (full rank):\n")
                print(tt$L.full)
            }
        }
        if (cov == T) {
            cat("\nVar-Cov of estimates:\n")
            print(tt$vcov)
            cat("\nCorrelations:\n")
            print(tt$vcor)
        }
    }
    invisible(x)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
