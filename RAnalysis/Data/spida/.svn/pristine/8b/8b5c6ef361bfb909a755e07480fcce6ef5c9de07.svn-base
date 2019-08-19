\name{brace}
\Rdversion{1.1}
\alias{brace}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a brace for a graph
}
\description{
The function \code{brace()} calculates the (X,Y) coordinates to draw a brace 
between two points in a graph.  
}
\usage{
brace(x1 = 0, y1 = 0, x2 = 0, y2 = 1, right = TRUE, rad = 0.2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x1, y1}{
Coordinate of the initial point for the brace
}
  \item{x2, y2}{
Coordinate of the final point for the brace
}
  \item{right}{
Logical to indicate if a vertical brace should have its' central point facing right.
If the brace is horizontal (\code{y1 == y2}), "right" means "down"
}
  \item{rad}{
Radius of the quarter circles used to draw the brace
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A 2-column matrix of (X,Y) points, suitable for use with \code{lines}
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Georges Monette, modified for R by Michael Friendly
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[graphics]{lines}}, \code{\link[graphics]{arrows}}, ~~~
}
\examples{
plot(c(-1,1), c(-1,1), type="n", xlab="x", ylab="y")
abline(h=0, col="gray")
abline(v=0, col="gray")
b <- 0.6
abline(0, b, col="blue")
lines(brace(0, 0, 0, b, right=FALSE, rad=0.1))
lines(brace(0, 0, 1, 0, rad=0.1))
text(0.5, -.2, '1', cex=2)
text(-.2, b/2, 'b', cex=2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{dplot}
\keyword{aplot}% __ONLY ONE__ keyword per line
