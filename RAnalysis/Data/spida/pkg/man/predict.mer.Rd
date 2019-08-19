\name{predict.mer}
\Rdversion{1.1}
\alias{predict.mer}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A tentative version of predict for mer objects
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
predict.mer(model, data = model.matrix(model), form, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{form}{
%%     ~~Describe \code{form} here~~
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
function (model, data = model.matrix(model), form, verbose = FALSE) 
{
    help = "\n    This is a crude predict method for class 'mer'\n    Limitations:\n    1) When invoked on the model it only returns the linear predictor for the\n       complete data. Note that a complete data frame is easily obtained with\n       model.frame( model )\n    2) When the 'data' argument is provided for 'new' data, it is also\n       necessary to provide the correct rhs of the fixed effects formula.\n\n    "
    if (missing(form) & !missing(data)) {
        cat(help)
        stop("Need 'form' if 'data' given")
    }
    if (!missing(data)) {
        data = model.matrix(form, data)
        cnames = colnames(data)
        if (verbose) 
            print(cnames)
        fnames = names(fixef(model))
        if (verbose) 
            print(fnames)
        if (any(cnames != fnames)) {
            cat("\nMatrix names:\n")
            print(cnames)
            cat("\nCoeff names:\n")
            print(fnames)
            warning("matrix and coeff names not the same")
        }
    }
    data \%*\% fixef(model)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
