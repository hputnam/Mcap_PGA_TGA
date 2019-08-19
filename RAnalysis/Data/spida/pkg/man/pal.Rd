\name{pal}
\Rdversion{1.1}
\alias{pal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate a palette of colours -- possibly superseded
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pal(col = c("blue", "pink"), border = "light gray", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{col}{
%%     ~~Describe \code{col} here~~
}
  \item{border}{
%%     ~~Describe \code{border} here~~
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
function (col = c("blue", "pink"), border = "light gray", ...) 
{
    n <- length(col)
    plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, 
        xlab = "", ylab = "", ...)
    rect(0, 0:(n - 1)/n, 0.6, 1:n/n, col = col, border = border)
    ret <- col2rgb(col)
    dimnames(ret)[[2]] <- col
    ret <- t(ret)
    txt <- paste(as.character(col), "(", apply(ret, 1, paste, 
        collapse = " "), ")")
    text(rep(0.6, n), (0:(n - 1) + 0.5)/n, txt, adj = 0)
    ret <- col2rgb(col)
    dimnames(ret)[[2]] <- col
    t(ret)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
