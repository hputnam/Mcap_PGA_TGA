\name{smspline.v}
\Rdversion{1.1}
\alias{smspline.v}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
smspline.v from library splines with better ginverse
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
smspline.v(time)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{time}{
%%     ~~Describe \code{time} here~~
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
function (time) 
{
    t1 <- sort(unique(time))
    p <- length(t1)
    h <- diff(t1)
    h1 <- h[1:(p - 2)]
    h2 <- h[2:(p - 1)]
    Q <- matrix(0, nr = p, nc = p - 2)
    Q[cbind(1:(p - 2), 1:(p - 2))] <- 1/h1
    Q[cbind(1 + 1:(p - 2), 1:(p - 2))] <- -1/h1 - 1/h2
    Q[cbind(2 + 1:(p - 2), 1:(p - 2))] <- 1/h2
    Gs <- matrix(0, nr = p - 2, nc = p - 2)
    Gs[cbind(1:(p - 2), 1:(p - 2))] <- 1/3 * (h1 + h2)
    Gs[cbind(1 + 1:(p - 3), 1:(p - 3))] <- 1/6 * h2[1:(p - 3)]
    Gs[cbind(1:(p - 3), 1 + 1:(p - 3))] <- 1/6 * h2[1:(p - 3)]
    Gs
    Zus <- t(ginv(Q))
    R <- chol(Gs, pivot = FALSE)
    tol <- max(1e-12, 1e-08 * mean(diag(R)))
    if (sum(abs(diag(R))) < tol) 
        stop("singular G matrix")
    Zvs <- Zus \%*\% t(R)
    list(Xs = cbind(rep(1, p), t1), Zs = Zvs, Q = Q, Gs = Gs, 
        R = R)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
