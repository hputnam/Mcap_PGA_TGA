\name{com}
\Rdversion{1.1}
\alias{com}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prediction with 'lmer' objects
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
com(fit, varpattern = "", form = terms(fit), data = getData(fit), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
%%     ~~Describe \code{fit} here~~
}
  \item{varpattern}{
%%     ~~Describe \code{varpattern} here~~
}
  \item{form}{
%%     ~~Describe \code{form} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
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
function (fit, varpattern = "", form = terms(fit), data = getData(fit), 
    ...) 
{
    getData <- function(x) UseMethod("getData")
    getData.lme <- function(x) x$data
    getData.lmer <- function(x) x@frame
    model.mat <- model.matrix(form, data, ...)
    fe <- fixef(fit)
    effect.names <- grep(varpattern, names(fe), value = T)
    ret <- c(model.mat[, effect.names] \%*\% fe[effect.names])
    names(ret) <- rownames(model.mat)
    retnames <- rownames(fit$resid)
    ret <- ret[retnames]
    attr(ret, "predictors") <- effect.names
    ret
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
