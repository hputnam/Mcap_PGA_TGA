\name{cell}
\Rdversion{1.1}
\alias{cell}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Confidence ellipse
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cell(...)
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
    UseMethod("cell")
    help <- "\n        See help for car::confidence.ellipse.lm\n        except that 'cell' returns the points to form the ellipse\n        which must be plotted with plot(...,type='l') or lines(...)\n        -- Use dfn to determine Sheffe dimension, i.e. dfn = 1 to generate ordinary CIs, dfn = 2 for 2-dim CE, etc.\n        -- TODO: extend to 3 dimensions if which.coef has length 3\n        "
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
