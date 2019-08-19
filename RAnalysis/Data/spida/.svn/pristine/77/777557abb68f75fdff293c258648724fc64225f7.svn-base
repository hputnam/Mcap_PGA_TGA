\name{smspline}
\Rdversion{1.1}
\alias{smspline}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate a smoothing spline
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
smspline(formula, data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
%%     ~~Describe \code{formula} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
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
function (formula, data) 
{
    if (is.vector(formula)) {
        x <- formula
    }
    else {
        if (missing(data)) {
            mf <- model.frame(formula)
        }
        else {
            mf <- model.frame(formula, data)
        }
        if (ncol(mf) != 1) 
            stop("formula can have only one variable")
        x <- mf[, 1]
    }
    x.u <- sort(unique(x))
    Zx <- smspline.v(x.u)$Zs
    Zx[match(x, x.u), ]
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
