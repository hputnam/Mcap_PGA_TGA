\name{Lall}
\Rdversion{1.1}
\alias{Lall}
\concept{wald}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hypothesis matrix for lmer objects
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Lall(fit, nam)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
%%     ~~Describe \code{fit} here~~
}
  \item{nam}{
%%     ~~Describe \code{nam} here~~
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
function (fit, nam) 
{
    if (class(fit) != "lmer") 
        stop("only implemented for lmer")
    v <- fit@frame[[nam]]
    if (!is.factor(v)) 
        stop("nam needs to specify the name of a factor")
    lev0 <- levels(v)[1]
    ret <- list()
    namf <- nam
    if (substring(namf, 1, 1) != "^") 
        namf <- paste("^", namf, sep = "")
    ret[[nam]] <- Lmat(fit, namf)
    ret[[paste(nam, "mu", sep = ".")]] <- Lmu(fit, nam)
    ret[[paste(nam, "diff", sep = ".")]] <- Ldiff(fit, nam)
    ret
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
