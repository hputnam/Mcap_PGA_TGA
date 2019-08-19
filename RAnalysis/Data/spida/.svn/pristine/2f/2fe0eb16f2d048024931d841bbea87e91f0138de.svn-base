\name{Lequal}
\Rdversion{1.1}
\alias{Lequal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hypothesis matrix for equality of factor level effects
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Lequal(fit, pat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
%%     ~~Describe \code{fit} here~~
}
  \item{pat}{
%%     ~~Describe \code{pat} here~~
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
function (fit, pat) 
{
    L <- Lmat(fit, pat)
    nam <- rownames(L)
    n <- nrow(L)
    if (n < 2) 
        return(L)
    plus <- unlist(apply(rbind(2:n), 2, seq, n))
    minus <- rep(1:(n - 1), (n - 1):1)
    Lp <- L[plus, ]
    Lm <- L[minus, ]
    Lret <- rbind(Lp - Lm)
    rn <- paste(nam[plus], nam[minus], sep = " - ")
    rownames(Lret) <- rn
    Lret
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
