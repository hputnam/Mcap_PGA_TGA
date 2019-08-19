\name{na.pad.lme}
\Rdversion{1.1}
\alias{na.pad.lme}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Add NAs -- lme
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
na.pad.lme(x, obj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{obj}{
%%     ~~Describe \code{obj} here~~
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
function (x, obj) 
{
    ind <- rep(NA, nrow(x$data))
    names(ind) <- rn <- rownames(x$data)
    ind.part <- 1:nrow(x$resid)
    names(ind.part) <- rownames(x$resid)
    ind[names(ind.part)] <- ind.part
    if (!is.null(dim(obj))) 
        ret <- obj[ind, ]
    else ret <- obj[ind]
    names(ret) <- rn
    ret
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
