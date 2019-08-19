\name{print.wald}
\Rdversion{1.1}
\alias{print.wald}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Print method for wald objects
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
print.wald(x, round = 6, pround = 5, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{round}{
%%     ~~Describe \code{round} here~~
}
  \item{pround}{
%%     ~~Describe \code{pround} here~~
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
function (x, round = 6, pround = 5, ...) 
{
    pformat <- function(x, digits = pround) {
        x <- format(xx <- round(x, digits))
        x[as.double(xx) == 0] <- paste(c("<.", rep("0", digits - 
            1), "1"), collapse = "")
        x
    }
    rnd <- function(x, digits) {
        if (is.numeric(x)) 
            x <- round(x, digits = digits)
        format(x)
    }
    for (ii in 1:length(x)) {
        nn <- names(x)[ii]
        tt <- x[[ii]]
        ta <- tt$anova
        ta[["p-value"]] <- pformat(ta[["p-value"]])
        print(as.data.frame(ta, row.names = nn))
        te <- tt$estimate
        rowlab <- attr(te, "labs")
        te[, "p-value"] <- pformat(te[, "p-value"])
        if (!is.null(round)) {
            for (ii in 1:length(te)) {
                te[[ii]] <- rnd(te[[ii]], digits = round)
            }
        }
        labs(te) <- rowlab
        print(te, digits = round, ...)
        cat("\n")
    }
    invisible(x)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
