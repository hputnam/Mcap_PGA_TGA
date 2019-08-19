\name{tran}
\Rdversion{1.1}
\alias{tran}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Translate elements of a vector
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
tran(from, to, x, tofactor = is.factor(x))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{from}{
%%     ~~Describe \code{from} here~~
}
  \item{to}{
%%     ~~Describe \code{to} here~~
}
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{tofactor}{
%%     ~~Describe \code{tofactor} here~~
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
function (from, to, x, tofactor = is.factor(x)) 
{
    if (is.factor(from)) 
        from <- as.character(from)
    if (is.factor(to)) 
        to <- as.character(to)
    to <- rep(to, length = length(from))
    ret <- x
    if (is.factor(x)) {
        ret <- as.character(ret)
        levs <- levels(x)
    }
    ret <- c(to, unique(ret))[match(ret, c(from, unique(ret)))]
    if (tofactor) {
        if (is.factor(x)) 
            tolevs <- tran(from, to, levs)
        else tolevs <- to
        tolevs <- c(tolevs, unique(ret))
        tolevs <- intersect(tolevs, unique(ret))
        ret <- factor(ret, levels = tolevs)
    }
    ret
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
