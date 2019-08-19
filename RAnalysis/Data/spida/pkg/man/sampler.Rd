\name{sampler}
\Rdversion{1.1}
\alias{sampler}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Show characters, colours, etc.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sampler(n = 24)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
%%     ~~Describe \code{n} here~~
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
function (n = 24) 
{
    old.par <- par(ask = T)
    on.exit(par(old.par))
    require(lattice)
    y <- 0:n
    x <- 0:n
    print(xyplot(y ~ x, type = "n", xlab = "lty", ylab = "col", 
        panel = function(x, y, ...) {
            for (i in x) {
                panel.xyplot(c(i, i), range(y), type = "l", lty = i, 
                  col = 1, lwd = 3)
            }
            for (i in y) {
                for (j in seq(0, 0.9, by = 0.1)) {
                  panel.xyplot(c(min(x) + j * (max(x) - min(x)), 
                    min(x) + (j + 0.1) * (max(x) - min(x))), 
                    c(i, i), type = "l", lty = 1, col = i, lwd = 3)
                }
            }
        }))
    spl <- trellis.par.get("superpose.line")
    z <- expand.grid(y = 1:length(spl$lty), x = 0:2)
    print(xyplot(y ~ x, z, ylim = c(0, length(spl$lty)), groups = y, 
        type = "b", main = "superpose.line and .symbol"))
    y <- 10 * (0:25)
    x <- 0:9
    print(xyplot(y ~ x, type = "n", main = "pch", xlab = expression(~alpha + 
        beta + gamma + delta[epsilon] + zeta^eta + theta + iota + 
        kappa), ylab = expression(~lambda + mu + nu + xi + omicron + 
        pi + rho + sigma + tau + upsilon + phi + chi + psi + 
        omega), panel = function(x, y, ...) {
        for (i in x) {
            for (j in y) {
                panel.xyplot(i, j, pch = i + j, cex = 2)
            }
        }
    }))
    invisible(0)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
