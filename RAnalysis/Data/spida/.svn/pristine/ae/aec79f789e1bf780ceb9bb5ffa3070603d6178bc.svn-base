\name{av.frame.lm}
\Rdversion{1.1}
\alias{av.frame.lm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
lm method for av.frame
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
av.frame.lm(model, variable, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{variable}{
%%     ~~Describe \code{variable} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (model, variable, ...) 
{
    mod.mat <- model.matrix(model)
    var.names <- colnames(mod.mat)
    omit <- grep(variable, var.names)
    if (0 == length(omit)) 
        stop(paste(variable, "is not matched among columns of the model matrix."))
    cat("x.var =", var.names[omit[1]], "\n", "omitted vars =", 
        var.names[omit[-1]], "\n")
    response <- response(model)
    x.var <- mod.mat[, omit[1]]
    Xpred <- mod.mat[, -omit]
    preds <- predict(update(model, na.action = na.exclude))
    responseName <- responseName(model)
    if (is.null(weights(model))) 
        wt <- rep(1, length(response))
    else wt <- weights(model)
    res <- lsfit(mod.mat[, -omit], cbind(mod.mat[, omit[1]], 
        response), wt = wt, intercept = FALSE)$residuals
    ret <- matrix(NA, nrow = length(preds), ncol = 2)
    ret[!is.na(preds), ] <- res
    data.frame(x.res = ret[, 1], y.res = ret[, 2])
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
