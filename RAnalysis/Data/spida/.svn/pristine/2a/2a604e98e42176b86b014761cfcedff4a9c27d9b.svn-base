\name{Ldiff.rdc}
\Rdversion{1.1}
\alias{Ldiff.rdc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Version of Ldiff used in RDC
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Ldiff.rdc(fit, nam, ref = "no longer used")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
%%     ~~Describe \code{fit} here~~
}
  \item{nam}{
%%     ~~Describe \code{nam} here~~
}
  \item{ref}{
%%     ~~Describe \code{ref} here~~
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
function (fit, nam, ref = "no longer used") 
{
    Lm <- Lmu(fit, nam)
    levs <- rownames(Lm)
    n <- nrow(Lm)
    if (n < 2) 
        return(Lm)
    plus <- unlist(apply(rbind(2:n), 2, seq, n))
    minus <- rep(1:(n - 1), (n - 1):1)
    Lret <- Lm[plus, ] - Lm[minus, ]
    rn <- paste(levs[plus], levs[minus], sep = " - ")
    rownames(Lret) <- rn
    Lret
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
