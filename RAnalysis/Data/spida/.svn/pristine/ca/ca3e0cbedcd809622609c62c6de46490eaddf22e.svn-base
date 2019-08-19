\name{Read.dta}
\Rdversion{1.1}
\alias{Read.dta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
read a STATA .dta file
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Read.dta(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
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
function (...) 
{
    help <- "\n  Read.dta reads Stata files using 'read.dta' in 'library(foreign)'\n  This appears to be an ideal way of importing spss files in order\n  to keep full variable names. Direct use of 'read.spss' on a SPSS\n  '.sav' file abbreviates variable names to 8 characters.\n  Note: missing.type = T produces warnings.\n"
    require("foreign")
    dd <- read.dta(...)
    cls <- sapply(dd, class)
    ch.nams <- names(dd)[cls == "character"]
    for (nn in ch.nams) dd[[nn]] <- factor(trim(dd[[nn]]))
    dd
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
