\name{select.first}
\Rdversion{1.1}
\alias{select.first}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find the first non-missing element of a vector and check whether there are other inconsistent non-missing values
}
\description{
Find the first non-missing element of a vector and check whether there are other inconsistent non-missing values
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
select.first(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
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
function (x, ...) 
{
    ret <- unique(x[!is.na(x)])
    if (length(ret) > 1 && is.character(ret) && ("" \%in\% ret)) 
        ret <- ret[ret != ""]
    if (length(ret) > 1) {
        cat("\n\n============= Multiple values in select.first ==========\n")
        for (i in 1:length(ret)) cat("\nValue", i, ": <<", ret[i], 
            ">>\n")
    }
    ret[1]
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
