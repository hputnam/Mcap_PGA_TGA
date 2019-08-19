\name{ell}
\Rdversion{1.1}
\alias{ell}
\alias{ellplus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate coordinates of an ellipse
}
\description{
\code{ell} is a utility function used to calculate the (X, Y) coordinates of a 2D ellipse
for the purpose of drawing statistical diagrams and plots.

\code{ellplus} can produce, in addition to the points of an ellipse, the
conjugate axes corresponding to a \code{chol} or other decomposition
and the surrounding parallelogram defined by these axes.


}
\usage{
ell(center, shape, radius = 1, n = 100)

ellplus(center = rep(0, 2), shape = diag(2), radius = 1, n = 100, angles = (0:n) * 2 * pi/n, fac = chol, ellipse = all, diameters = all, box = all, all = FALSE)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{center}{
X,Y location of the center of the ellipse
}
  \item{shape}{
A 2x2 matrix, typically a covariance matrix of data (for a data ellipse), or a
covariance matrix of estimated parameters in a model (for a confidence ellipse).
}
  \item{radius}{
Radius of the ellipse-generating unit circle.  The default, \code{radius=1}
corresponds to a "standard" ellipse.
}
  \item{n}{
Number of points on the unit circle used to calculate the ellipse
}
  \item{angles}{
Angles around the unit circle used to calculate the ellipse
}
  \item{fac}{
A function defining the conjugate axes used to transform the unit circle
into an ellipse.  The default, \code{chol}, uses the right Cholesky factor
of \code{shape}.
}
  \item{ellipse}{
Logical to indicate if the points on the ellipse should be returned
}
  \item{diameters}{
Logical to indicate if the points defining the ends of the
conjugate axes of the ellipse should be returned
}
  \item{box}{
Logical to indicate if the points on the conjugate-axes bounding box should be returned
}
  \item{all}{
Logical to request all of \code{ellipse}, \code{diameters} and \code{box}.
If \code{FALSE}, only the components specified separately by 
\code{ellipse}, \code{diameters} and \code{box} are returned.
}

}
%\details{
%}
\value{
Returns a 2-column matrix of (X,Y) coordinates suitable for drawing with
\code{lines()}.  

For \code{ellplus}, when more than one of the options \code{ellipse}, \code{diameters},
and \code{box} is \code{TRUE}, the different parts are separated by a row of \code{NA}.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Georges Monette
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{cell}}, \code{\link{dell}}, \code{\link{dellplus}},
}
\examples{
plot( x=0,y=0, xlim = c(-3,3), ylim = c(-3,3),
      xlab = '', ylab = '', type = 'n', asp=1)
abline( v=0, col="gray")
abline( h=0, col="gray")
A <- cbind( c(1,2), c(1.5,1))
W <- A \%*\% t(A)    
    
lines( ell(center=c(0,0), shape = W ), col = 'blue', lwd=3)
lines( ellplus(center=c(0,0), shape = W, box=TRUE, diameters=TRUE ), col = 'red')

# show conjugate axes for PCA factorization
pca.fac <- function(x) {
    xx <- svd(x)
    ret <- t(xx$v) * sqrt(pmax( xx$d,0))
    ret 
}

plot( x=0,y=0, xlim = c(-3,3), ylim = c(-3,3),
      xlab = '', ylab = '', type = 'n', asp=1)
abline( v=0, col="gray")
abline( h=0, col="gray")
lines( ell(center=c(0,0), shape = W ), col = 'blue', lwd=3)
lines( ellplus(center=c(0,0), shape = W, box=TRUE, diameters=TRUE, fac=pca.fac ), col = 'red')

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{dplot}
\keyword{aplot}% __ONLY ONE__ keyword per line
