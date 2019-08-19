\name{dell}
\Rdversion{1.1}
\alias{dell}
\alias{dellplus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate coordinates of a data ellipse
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
\code{dell} to calculates the coordinates of a 2D  data ellipse
(concentration ellipse) from (X, Y) variables.

\code{dellplus} can produce, in addition to the points of an ellipse, the
conjugate axes corresponding to a \code{chol} or other decomposition
and the surrounding parallelogram defined by these axes.

}
\usage{
dell(x, y, radius = 1, ...)

dellplus(x, y, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, y}{
Either a two-column matrix or
numeric vectors of the same length
}

  \item{radius}{
Radius of the ellipse-generating unit circle.  The default, \code{radius=1}
corresponds to a "standard" ellipse.
}
  \item{\dots}{
Other arguments passed down to \code{ell} or \code{ellplus}.
}

}

\details{
These functions simply calculate the mean vector and
covariance matrix and call \code{ell} or \code{ellplus}.
}

\value{
Returns a 2-column matrix of (X,Y) coordinates suitable for drawing with
\code{lines()}.  

For \code{dellplus}, when more than one of the options \code{ellipse}, \code{diameters},
and \code{box} is \code{TRUE}, the different parts are separated by a row of \code{NA}.
}
\references{
Monette, G. (1990).
Geometry of Multiple Regression and Interactive 3-D Graphics.
In Fox, J. & Long, S. (ed.) 
\emph{Modern Methods of Data Analysis}, Sage Publications, 209-256.

}
\author{
Georges Monette
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{cell}}, \code{\link{ell}}, \code{\link{ellplus}},
}
\examples{
data(Prestige)   # from car
attach(Prestige)
fit.simple <- lm( prestige ~ education, Prestige)

plot(prestige ~ education, type='p')
lines(dell(education, prestige), col="blue", lwd=3)
lines(bbox <- dellplus(education, prestige, box=TRUE))
lines(dellplus(education, prestige, diameter=TRUE, radius=2), col="gray")
detach(Prestige)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{dplot}
\keyword{aplot}% __ONLY ONE__ keyword per line
