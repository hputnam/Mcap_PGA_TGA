\name{as.data.frame.wald}
\Rdversion{1.1}
\alias{as.data.frame.wald}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Transform output of a Wald test into a data frame
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
as.data.frame.wald(obj, se = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
%%     ~~Describe \code{obj} here~~
}
  \item{se}{
%%     ~~Describe \code{se} here~~
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
function (obj, se = FALSE) 
{
    if (length(obj) == 1) {
        if (is.null(dd <- obj[[1]]$data)) 
            return(NULL)
        cf <- obj[[1]]$coef
        ret <- if (is.logical(se)) {
            if (se) 
                data.frame(dd, coef = cf, se = obj[[1]]$se)
            else data.frame(dd, coef = cf)
        }
        else if (se > 0) {
            data.frame(dd, coef = cf, coefp = cf + se * obj[[1]]$se, 
                coefm = cf - se * obj[[1]]$se, se = obj[[1]]$se)
        }
    }
    else ret <- sapply(obj, as.data.frame.wald)
    ret
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
