\name{Lc}
\Rdversion{1.1}
\alias{Lc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hypothesis matrix for lmer objects: comparisons with reference level
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Lc(fit, nam, ref = 1, verbose = 0)
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
function (fit, nam, ref = 1, verbose = 0) 
{
    if (class(fit) != "lmer") 
        stop("only implemented for lmer")
    L <- Lmu(fit, nam)
    Lref <- L[ref, , drop = F]
    index <- 1:nrow(L)
    names(index) <- rownames(L)
    refind <- index[ref]
    if (length(refind) != 1) 
        stop(paste(ref, "does not refer to a single level"))
    Lret <- L[-refind, ]
    Lret <- Lret - cbind(rep(1, nrow(Lret))) \%*\% Lref
    attr(Lret, "heading") <- paste("Comparisons with reference level:", 
        rownames(L)[refind])
    Lret
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
