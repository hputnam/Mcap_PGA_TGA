\name{p}
\Rdversion{1.1}
\alias{p}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
p - quick paste with sep = ''
}
\description{
Works like \code{paste}, using an empty separator string.
}
\usage{
p(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{
one or more R objects, to be converted to character vectors.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A character vector of the concatenated values. 
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
\code{\link[base]{paste}}
}
\examples{
p(letters[1:5], 1:5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
