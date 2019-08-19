\name{tab}
\Rdversion{1.1}
\alias{tab}
\alias{tab.formula}
\alias{tab.data.frame}
\alias{tab.default}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Table of frequencies or relative frequencies bordered with totals and including NAs
}
\description{
Generates a table of frequencies or relative frequencies or relative percentages
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{

tab( df, ~ a + b )

tab( df, ~ a + b , pct = 1)

tab( ~ a + b + c , df, pr = c(1,2))

tab( a, b, c )

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{
as with \code{table}, one or more objects which can be interpreted as factors (including character strings),
or a list (or data frame) whose components can be so interpreted.
}
\item{useNA}{
 whether to include NA levels. The default is "ifany". Can also be set to "no" or "always".
}
\item{pct}{
margins to be scaled to sum to 100
}
\item{pct}{
margins to be scaled to sum to 1
}
\item{total.margins}{
if FALSE, generate table without margins
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A matrix or an array of dimension equal to the number of variable with margins showing totals
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
    UseMethod("tab")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
