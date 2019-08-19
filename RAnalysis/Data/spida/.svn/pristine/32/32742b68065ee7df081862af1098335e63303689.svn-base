\name{diags.lm}
\Rdversion{1.1}
\alias{diags.lm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Standard diagnostics for lm objects
}
\description{
Standard diagnostics for lm objects
}
\usage{
diags.lm(x, y, ..., ask, labels = names(residuals(x)), showlabs = text)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
  \item{ask}{
%%     ~~Describe \code{ask} here~~
}
  \item{labels}{
%%     ~~Describe \code{labels} here~~
}
  \item{showlabs}{
%%     ~~Describe \code{showlabs} here~~
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
function (x, y, ..., ask, labels = names(residuals(x)), showlabs = text) 
{
    if (!missing(ask)) {
        op <- par(ask = ask)
        on.exit(par(op))
    }
    form <- formula(x)
    f <- predict(x)
    r <- residuals(x)
    nams <- names(r)
    if (!missing(labels)) {
        nams <- names(residuals(x))
        if (length(nams) != length(labels)) 
            labels <- labels[nams]
    }
    ret <- NULL
    if (missing(y)) {
        y <- f + r
        yname <- deparse(form[[2]])
    }
    else yname <- deparse(substitute(y))
    fname <- paste("Fitted:", deparse(form[[3]]), collapse = " ")
    plot(f, y, xlab = fname, ylab = yname, main = "Dependent var. vs. Predicted", 
        ...)
    abline(0, 1, lty = 1)
    lines(supsmu(f, y))
    showlabs(f, y, labels, ...)
    lmi <- lm.influence(x)
    hat <- lmi$hat
    sigma <- lmi$sigma
    mm <- scale(model.matrix(x), scale = F)
    mp <- predict(x, type = "terms")
    comp.res <- mp + r
    plot(f, abs(r), xlab = fname, ylab = deparse(substitute(abs(resid(x)))), 
        main = "Absolute Residual vs. Predicted", ...)
    showlabs(f, abs(r), labels, ...)
    zq <- qqnorm(r, main = "Normal Quantile Plot", ylab = "Residual", 
        sub = fname)
    qqline(r)
    showlabs(zq, labels, ...)
    n <- length(r)
    r.o <- sort(r)
    half <- (n + 1)/2
    if (n\%\%2 == 1) {
        med <- r.o[half]
        below <- med - r.o[half:1]
        above <- r.o[half:n] - med
    }
    else {
        med <- sum(r.o[c(half, half + 1)])/2
        below <- med - r.o[(n/2):1]
        above <- r.o[(n/2 + 1):n] - med
    }
    opt <- par(pty = "s")
    ran <- range(c(below, above))
    plot(below, above, main = "Symmetry plot of residuals", xlab = "Distance below median", 
        ylab = "Distance above median", xlim = ran, ylim = ran)
    abline(0, 1, lty = 2)
    par(opt)
    std.r <- r/(sigma * sqrt(1 - hat))
    plot(hat, std.r, xlab = "Leverage (hat)", ylab = yname, sub = fname, 
        main = "Studentized residual vs. Leverage", ...)
    showlabs(hat, std.r, labels, ...)
    nams <- dimnames(lmi$coefficients)[[1]]
    pairs(lmi$coefficients)
    pairs(lmi$coefficients, panel = function(x, y, nams) {
        points(x, y)
        text(x, y, nams)
    }, nams = nams)
    invisible(0)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
