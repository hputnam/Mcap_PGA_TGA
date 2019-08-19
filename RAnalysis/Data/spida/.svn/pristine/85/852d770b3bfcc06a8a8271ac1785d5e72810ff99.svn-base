\name{Q}
\Rdversion{1.1}
\alias{Q}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Q matrix of QR decomposion with missing data
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Q(x, verbose = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
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
function (x, verbose = 0) 
{
    miss <- apply(x, 1, function(xx) any(is.na(xx)))
    xc <- x[!miss, ]
    qf <- qr.Q(qqr <- qr(xc))
    if (verbose > 0) {
        cat("xc:", dim(xc), "\n")
        cat("qf:", dim(qf), "\n")
        print(qf)
    }
    if (ncol(xc) > ncol(qf)) 
        xc <- xc[, 1:ncol(qf)]
    ip <- sign(apply(qf * xc, 2, sum))
    qf <- qf * rep(ip, rep(nrow(qf), length(ip)))
    ret <- array(NA, dim = c(nrow(x), ncol(qf)))
    rownames(ret) <- rownames(x)
    colnames(ret) <- colnames(xc)
    ret[!miss, ] <- qf
    attr(ret, "rank") <- qqr$rank
    attr(ret, "miss") <- miss
    ret
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
