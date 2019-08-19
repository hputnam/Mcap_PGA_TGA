\name{Ldiff}
\Rdversion{1.1}
\alias{Ldiff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hypothesis matrix for ...
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Ldiff(fit, pat, levnames = c(reflevel, substring(rownames(L), cut + 1)), reflevel = "<ref>", cut = nchar(pat), verbose = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
%%     ~~Describe \code{fit} here~~
}
  \item{pat}{
%%     ~~Describe \code{pat} here~~
}
  \item{levnames}{
%%     ~~Describe \code{levnames} here~~
}
  \item{reflevel}{
%%     ~~Describe \code{reflevel} here~~
}
  \item{cut}{
%%     ~~Describe \code{cut} here~~
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
function (fit, pat, levnames = c(reflevel, substring(rownames(L), 
    cut + 1)), reflevel = "<ref>", cut = nchar(pat), verbose = F) 
{
    L <- Lmat(fit, pat)
    nam <- rownames(L)
    n <- nrow(L)
    zm <- matrix(1:n, nrow = n, ncol = n)
    plus <- zm[col(zm) < row(zm)]
    minus <- rep(1:(n - 1), (n - 1):1)
    Lp <- L[plus, ]
    Lm <- L[minus, ]
    Lret <- rbind(L, Lp - Lm)
    pnames <- levnames[c(1:n, plus) + 1]
    mnames <- levnames[c(rep(0, n), minus) + 1]
    if (verbose) {
        print(levnames)
        print(plus)
        print(minus)
        print(Lp)
        print(Lm)
        print(L)
        print(Lret)
        print(pnames)
        print(mnames)
    }
    rn <- paste(levnames[c(1:n, plus) + 1], levnames[c(rep(0, 
        n), minus) + 1], sep = " - ")
    rownames(Lret) <- rn
    Lret
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
