\name{gicc}
\alias{gicc}
\alias{gicc.default}
\alias{gicc.data.frame}
\alias{gicc.factor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generalized intraclass correlation coefficient for multilevel data
%%  ~~function to do ... ~~
}
\description{
\code{gicc} facilitates the identification of between-subject variables versus balanced within-subject variables with hierarchical data.
Methods handle factors, character variables and numerical variables. For categorical variables, Goodman-Kruskal's tau is returned
and for numerical variables, the simple ICC with (variance between) / (variance between + variance within) applied to the rank of variable.
In either case, a value of 1
signifies a variable that is constant within clusters and a value of 0, a variable that is perfectly balanced within clusters.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gicc(x, by, method = 'raw', ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    a data frame, factor, character variable or a numerical variable
%%     ~~Describe \code{x} here~~
}
  \item{by}{
    if \code{x} is a data frame, \code{by} can be a formula (e.g. \code{~id}) evaluated in \code{x}. Otherwise,   \code{by} is a variable
    defining clusters.
%%     ~~Describe \code{by} here~~
}
  \item{method}{
    a character string indicating whether to work with the raw variable, \code{x}, stripped of missing values, or with \code{is.na(x)}. Can be one of "raw" or "na".

%%     ~~Describe \code{method} here~~
}
  \item{\dots}{
  -- not used
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a measure of relative variability within clusters so that 1 represents no variability and 0 perfect balance.
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
G. Monette <georges@yorku.ca>
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
function(x, by, method = NULL, ...) {

  # Computes a 'generalized intra-class correlation coefficient' using
  # for factors: Goodman and Kruskal's tau based on proportional reduction in
  #    prediction error probability from a proportional prediction rule, and
  # for numerical variables: intraclass correlation coefficient(1,1) based
  #    on ranks of x
  #
  # Args:
  #   x: a variable or a data frame
  #   by: a variable to form cluster, or if x is a data frame, a formula
  #   method: not yet implemented
  # Returns:
  #   The gicc for a variable or for each variable in a data frame
  # TODO: implement different variation indices to see which discriminate
  #       most usefully (e.g. see p. 25 of Agresti (1990) Categorical Data Analysis, Wiley.
    UseMethod("gicc")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
