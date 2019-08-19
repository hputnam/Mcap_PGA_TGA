\name{na.include}
\Rdversion{1.1}
\alias{na.include}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Include NAs
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
na.include(obj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
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
function (obj) 
{
    if (inherits(obj, "data.frame")) 
        for (i in seq(along = obj)) obj[[i]] <- na.include(obj[[i]])
    else {
        if (length(levels(obj)) && any(is.na(obj))) 
            obj <- factor(obj, exclude = NULL)
    }
    obj
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
