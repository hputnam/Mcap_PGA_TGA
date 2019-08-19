\name{eg}
\Rdversion{1.1}
\alias{eg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Attempt at streamlining expand.grid for prediction
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
eg(df, by, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{
%%     ~~Describe \code{df} here~~
}
  \item{by}{
%%     ~~Describe \code{by} here~~
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
function (df, by, ...) 
{
    dots = list(...)
    by = model.frame(by, df, na.action = na.include)
    byn = names(by)
    names(byn) = byn
    args = lapply(byn, function(x) {
        vv = df[[x]]
        if (is.factor(vv)) 
            levels(vv)
        else unique(vv)
    })
    args = c(args, dots)
    do.call("expand.grid", args)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
