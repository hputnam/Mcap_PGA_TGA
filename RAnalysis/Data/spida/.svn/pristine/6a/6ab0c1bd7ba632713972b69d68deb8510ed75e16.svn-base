\name{constant.data.frame}
\Rdversion{1.1}
\alias{constant.data.frame}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identify variables that are constant within levels of a factor
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
constant.data.frame(x, id, all = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{id}{
%%     ~~Describe \code{id} here~~
}
  \item{all}{
%%     ~~Describe \code{all} here~~
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
function (x, id, all = FALSE, ...) 
{
    if (missing(id)) 
        ret <- sapply(x, constant, ...)
    else {
        id <- eval(substitute(id), x, parent.frame())
        if (inherits(id, "formula")) 
            id <- c(model.frame(id, x))
        ret <- sapply(x, function(xx) tapply(xx, id, constant, 
            ...))
        if (all) 
            ret <- apply(ret, 2, all)
    }
    ret
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
