\name{cell.glh}
\Rdversion{1.1}
\alias{cell.glh}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Confidence ellipse for estimated parameters in Wald test
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cell.glh(obj, which.coef = 1:2, levels = 0.95, Scheffe = FALSE, dfn = 2, center.pch = 19, center.cex = 1.5, segments = 51, xlab, ylab, las = par("las"), col = palette()[2], lwd = 2, lty = 1, add = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
%%     ~~Describe \code{obj} here~~
}
  \item{which.coef}{
%%     ~~Describe \code{which.coef} here~~
}
  \item{levels}{
%%     ~~Describe \code{levels} here~~
}
  \item{Scheffe}{
%%     ~~Describe \code{Scheffe} here~~
}
  \item{dfn}{
%%     ~~Describe \code{dfn} here~~
}
  \item{center.pch}{
%%     ~~Describe \code{center.pch} here~~
}
  \item{center.cex}{
%%     ~~Describe \code{center.cex} here~~
}
  \item{segments}{
%%     ~~Describe \code{segments} here~~
}
  \item{xlab}{
%%     ~~Describe \code{xlab} here~~
}
  \item{ylab}{
%%     ~~Describe \code{ylab} here~~
}
  \item{las}{
%%     ~~Describe \code{las} here~~
}
  \item{col}{
%%     ~~Describe \code{col} here~~
}
  \item{lwd}{
%%     ~~Describe \code{lwd} here~~
}
  \item{lty}{
%%     ~~Describe \code{lty} here~~
}
  \item{add}{
%%     ~~Describe \code{add} here~~
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
function (obj, which.coef = 1:2, levels = 0.95, Scheffe = FALSE, 
    dfn = 2, center.pch = 19, center.cex = 1.5, segments = 51, 
    xlab, ylab, las = par("las"), col = palette()[2], lwd = 2, 
    lty = 1, add = FALSE, ...) 
{
    obj <- obj[[1]]
    coef <- obj$coef[which.coef]
    xlab <- if (missing(xlab)) 
        paste(names(coef)[1], "coefficient")
    ylab <- if (missing(ylab)) 
        paste(names(coef)[2], "coefficient")
    dfd <- obj$anova$denDF
    shape <- obj$vcov[which.coef, which.coef]
    ret <- ell(coef, shape, sqrt(dfn * qf(levels, dfn, dfd)))
    colnames(ret) <- c(xlab, ylab)
    ret
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
