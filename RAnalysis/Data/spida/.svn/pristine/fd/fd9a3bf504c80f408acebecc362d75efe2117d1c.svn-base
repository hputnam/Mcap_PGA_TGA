\name{up}
\Rdversion{1.1}
\alias{up}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a data frame at a higher level of aggregation
}
\description{
Produce a higher level data set with one row per cluster. The data set can contain only
variables that are invariant in each cluster or it can also include summaries (mean or modes)
of variables that vary by cluster. Adapted from \code{gsummary} in the \code{nlme} package.
}
\usage{
up(object, form = formula(object), all = FALSE, FUN = function(x) mean(x, na.rm = TRUE), omitGroupingFactor = FALSE, groups, invariantsOnly = !all, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
a data frame to be summarized.
}
  \item{form}{
a one-sided formula identifying the variable(s) in \code{object} that identifies clusters. e.g. ~ school/Sex to
get a summary within each Sex of each school.
}
  \item{all}{
if TRUE, include summaries of variables that vary within clusters, otherwise keep
only cluster-invariant variables.
}
  \item{FUN}{
function to be used for summaries.
}
  \item{omitGroupingFactor}{
kept for compatibility with \code{gsummary}
}
  \item{groups}{
kept for compatibility with \code{gsummary}
}
  \item{invariantsOnly}{
kept for compatibility with \code{gsummary}
}
  \item{\dots}{
additional arguments to \code{tapply} when summarizing numerical variables. e.g. \code{na.rm = TRUE}
}
}
\details{
\code{up} was created from \code{nlme::gsummary} and modified to make it easier to
use and to make an equivalent of \code{gsummary} available when using \code{lme4}.
}
\value{
a data frame with one row per value of the variable in \code{form}
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

    data(hs)
    dim( hs )
    hsu <- up( hs, ~ school )
    dim( hsu )
    
    # to also get cluster means of cluster-varying numeric variables and modes of factors:

    hsa <- up( hs, ~ school , all = TRUE )

    # to get summary proportions of cluster varying factors:

    up( cbind( hs, model.matrix( ~ Sex -1 , hs)), ~ school, all = T)


    ## To plot a summary between-cluster panel along with within-cluster panels:

    hsu <- up( hs, ~ school, all = TRUE)
    hsu$school <- ' between'  # space to make it come lexicographically before cluster names

    require( lattice )
    xyplot( mathach ~ ses | school, rbind(hs,hsu),
        panel = function( x, y, ...) {
            panel.xyplot( x, y, ...)
            panel.lmline( x, y, ...)
        } )

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
