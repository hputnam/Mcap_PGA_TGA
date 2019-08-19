\name{av.frame}
\Rdversion{1.1}
\alias{av.frame}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create data frame for added variable plot
}
\description{
av.frame( model, variable)
returns a data frame with model.frame(model) augmented
by y.res and x.res, the residuals for an added variable
plot

The purpose of this function is to facilitate OLS av.plots
for mixed models.
}
\usage{
av.frame(model, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
%%     ~~Describe \code{model} here~~
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
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
       library(nlme)
       library(lattice)
       hs <- read.csv( 'http://www.math.yorku.ca/~georges/Data/hs.csv')

       # Mixed model where ses and Sex are Level 1 and Sector is Level 2

       fit.mm <- lme( mathach ~ ses * Sex * Sector, hs, random = ~ 1+ses| school)

       # for diagnostics fit an OLS model using only level 1 variables interacting
       # with the id variable

       fit.ols <- lm( mathach ~ (ses * Sex ) * factor(school), hs)
       xyplot( y.res ~ x.res | factor(school), cbind(av.frame(fit.ols, 'ses:Sex'),hs), sub = 'ses:Sex')
       xyplot( y.res ~ x.res | factor(school), cbind(av.frame(fit.ols, '^Sex'),hs), sub = 'Sex')
       xyplot( y.res ~ x.res | factor(school), cbind(av.frame(fit.ols, '^ses$|^ses:f'),hs), sub = 'ses')


       Note : y.res is the residual from fitting the response on
              the model matrix for fit.ols omitting any column
              whose names is matched (as a regular expression)
              by 'effect'
              x.res is the residual of the first column of the
              model matrix that is matched by 'effect' on the
              same matrix used for y.res.
       Caution: To make sure that the correct columns were
              matched, the list of matched columns that are omitted
              is printed.

## The function is currently defined as
function (model, ...) 
{
    UseMethod("av.frame")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
