\name{Lmat}
\Rdversion{1.1}
\alias{Lmat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate a hypothesis matrix 
}
\description{
Generates a hypothesis matrix to 
test whether a group of coefficients in a linear model are
jointly zero, selecting the coefficients that match a regular expression.
}
\usage{
Lmat(fit, pattern)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
%%     ~~Describe \code{fit} here~~
}
  \item{pattern}{
%%     ~~Describe \code{pattern} here~~
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
function (fit, pattern) 
{
    umatch <- function(pat, x) {
        ret <- rep(0, length(pat))
        for (ii in 1:length(pat)) {
            imatch <- grep(pat[ii], x)
            if (length(imatch) != 1) {
                cat("\nBad match of:\n")
                print(pat)
                cat("in:\n")
                print(x)
                stop("Bad match")
            }
            ret[ii] <- imatch
        }
        ret
    }
    if (is.character(fit)) {
        x <- pattern
        pattern <- fit
        fit <- x
    }
    fe <- getFix(fit)$fixed
    ne <- names(fe)
    if (is.character(pattern)) {
        L.indices <- grep(pattern, names(fe))
        ret <- diag(length(fe))[L.indices, , drop = F]
        rownames(ret) <- names(fe)[L.indices]
        labs(ret) <- "Coefficients"
    }
    else if (is.list(pattern)) {
        ret <- matrix(0, nrow = length(pattern), ncol = length(fe))
        colnames(ret) <- ne
        for (ii in 1:length(pattern)) {
            Lcoefs <- pattern[[ii]]
            pos <- umatch(names(Lcoefs), ne)
            if (any(is.na(pos))) 
                stop("Names of L coefs not matched in fixed effects")
            ret[ii, pos] <- Lcoefs
        }
        rownames(ret) <- names(pattern)
    }
    labs(ret) <- "Coefficients"
    ret
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
