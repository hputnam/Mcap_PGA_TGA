\name{Lmu}
\Rdversion{1.1}
\alias{Lmu}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate predicted response for a factor level.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Lmu(fit, nam, verbose = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
%%     ~~Describe \code{fit} here~~
}
  \item{nam}{
%%     ~~Describe \code{nam} here~~
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
function (fit, nam, verbose = 0) 
{
    if (class(fit) != "lmer") 
        stop("only implemented for lmer")
    v <- fit@frame[[nam]]
    if (!is.factor(v)) 
        stop("nam needs to specify the name of a factor")
    levs <- levels(v)
    if (verbose > 0) 
        print(levs)
    cmat <- contrasts(v)
    if (verbose > 0) 
        print(cmat)
    fe <- getFix(fit)$fixed
    if (verbose > 0) 
        print(fe)
    if (substring(nam, 1, 1) != "^") 
        nam <- paste("^", nam, sep = "")
    L.indices <- grep(nam, names(fe))
    if (verbose > 0) 
        print(L.indices)
    L <- matrix(0, nrow = length(levs), ncol = length(fe))
    colnames(L) <- names(fe)
    if (verbose > 0) 
        print(L)
    rownames(L) <- levs
    L[, L.indices] <- cmat
    if ("(Intercept)" \%in\% colnames(L)) 
        L[, "(Intercept)"] <- 1
    L
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
