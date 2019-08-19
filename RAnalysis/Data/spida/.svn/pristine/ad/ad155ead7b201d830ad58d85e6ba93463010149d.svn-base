\name{xanova.lmer}
\Rdversion{1.1}
\alias{xanova.lmer}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Modified Anova for lmer objects
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xanova.lmer(fit, Llist, df = NULL, clevel = 0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
%%     ~~Describe \code{fit} here~~
}
  \item{Llist}{
%%     ~~Describe \code{Llist} here~~
}
  \item{df}{
%%     ~~Describe \code{df} here~~
}
  \item{clevel}{
%%     ~~Describe \code{clevel} here~~
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
function (fit, Llist, df = NULL, clevel = 0.95) 
{
    warning("xanova.lmer uses Chi-Square tests")
    ret <- list()
    for (ii in 1:length(Llist)) {
        L <- rbind(Llist[[ii]])
        QR <- qr(L)
        R <- qr.R(QR)
        dfH <- QR$rank
        eta <- R \%*\% fixef(fit)
        vv <- R \%*\% vcov(fit) \%*\% t(R)
        chisq <- t(eta) \%*\% qr.solve(vv, eta)
        test <- list(ChiSquare = chisq, DF = dfH, `p-value` = 1 - 
            pchisq(chisq, dfH))
        ret[[ii]]$anova <- test
        eta <- L \%*\% fixef(fit)
        vv <- diag(L \%*\% vcov(fit) \%*\% t(L))
        etasd <- sqrt(vv)
        zval <- c(eta/etasd)
        aod <- cbind(Estimate = c(eta), Std.Error = etasd, `z-value` = zval, 
            `p-value` = 2 * pnorm(-abs(zval)))
        if (!is.null(clevel)) {
            hw <- qnorm(1 - (1 - clevel)/2) * etasd
            aod <- cbind(aod, LL = eta - hw, UL = eta + hw)
            labs <- paste(c("Lower", "Upper", format(clevel)))
            colnames(aod)[ncol(aod) + c(-1, 0)] <- labs
        }
        aod <- as.data.frame(aod)
        class(aod) <- c("estimate.lme", "data.frame")
        ret[[ii]]$estimate <- aod
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
