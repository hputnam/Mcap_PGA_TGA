\name{vif.lme}
\Rdversion{1.1}
\alias{vif.lme}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Variance Inflation Factors for Mixed Models
}
\description{
Calculates versions of the variance-inflation and generalized variance-inflation factors for mixed models.
}
\usage{
vif.lme(mod)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mod}{
%%     ~~Describe \code{mod} here~~
}
}
\details{
The concept of Variance Inflation in linear models can be applied to mixed models in a
number of ways since the variance-covariance matrix of the estimated fixed coefficients is not
simply proportional to the inverse of the cross-product matrix for the data. This
method for the generic function \code{vif} in the \code{car} package, implements the version
based on the variance-covariance funtion of the estimated fixed coefficients.  Since the \code{vcov}
method is available for \code{lm} objects and \code{lme} objects, uses the code for
the \code{lm} method for \code{lme} objects.
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
function (mod) 
{
    if (any(is.na(fixef(mod)))) 
        stop("there are aliased coefficients in the model")
    v <- vcov(mod)
    mm <- model.matrix(formula(mod), mod$data)
    assign <- attributes(mm)$assign
    if (names(fixef(mod)[1]) == "(Intercept)") {
        v <- v[-1, -1]
        assign <- assign[-1]
    }
    else warning("No intercept: vifs may not be sensible.")
    terms <- labels(terms(mod))
    n.terms <- length(terms)
    if (n.terms < 2) 
        stop("model contains fewer than 2 terms")
    R <- cov2cor(v)
    detR <- det(R)
    result <- matrix(0, n.terms, 3)
    rownames(result) <- terms
    colnames(result) <- c("GVIF", "Df", "GVIF^(1/2Df)")
    for (term in 1:n.terms) {
        subs <- which(assign == term)
        result[term, 1] <- det(as.matrix(R[subs, subs])) * det(as.matrix(R[-subs, 
            -subs]))/detR
        result[term, 2] <- length(subs)
    }
    if (all(result[, 2] == 1)) 
        result <- result[, 1]
    else result[, 3] <- result[, 1]^(1/(2 * result[, 2]))
    result
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
