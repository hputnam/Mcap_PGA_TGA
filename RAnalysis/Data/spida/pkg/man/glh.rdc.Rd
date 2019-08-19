\name{glh.rdc}
\Rdversion{1.1}
\alias{glh.rdc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Test version of glh in RDC
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
glh.rdc(fit, Llist, help = F, clevel = 0.95, debug = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
%%     ~~Describe \code{fit} here~~
}
  \item{Llist}{
%%     ~~Describe \code{Llist} here~~
}
  \item{help}{
%%     ~~Describe \code{help} here~~
}
  \item{clevel}{
%%     ~~Describe \code{clevel} here~~
}
  \item{debug}{
%%     ~~Describe \code{debug} here~~
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
function (fit, Llist, help = F, clevel = 0.95, debug = F) 
{
    if (help) {
        cat("help!!!")
        return(0)
    }
    if (!is.list(Llist)) 
        Llist <- list(Llist)
    ret <- list()
    fix <- getFix(fit)
    beta <- fix$fixed
    vc <- fix$vcov
    dfs <- fix$df
    for (ii in 1:length(Llist)) {
        ret[[ii]] <- list()
        L <- rbind(zz <- Llist[[ii]])
        heading <- attr(zz, "heading")
        nam <- names(Llist)[ii]
        qqr <- qr(t(L))
        L.rank <- qqr$rank
        L.full <- t(qr.Q(qqr))[1:L.rank, , drop = F]
        if (F) {
            cat("\n+++++++++++++++++++\n")
            print(nam)
            print(L)
            print(svd(L)$d)
            print(L.full)
            print(svd(L.full)$d)
        }
        if (debug) {
            print(L)
            print(dim(L.full))
            print(dim(vc))
        }
        vv <- L.full \%*\% vc \%*\% t(L.full)
        eta.hat <- L.full \%*\% beta
        Fstat <- (t(eta.hat) \%*\% qr.solve(vv, eta.hat))/L.rank
        Fstat2 <- (t(eta.hat) \%*\% solve(vv) \%*\% eta.hat)/L.rank
        included.effects <- apply(L, 2, function(x) sum(abs(x))) != 
            0
        denDF <- min(dfs[included.effects])
        numDF <- L.rank
        ret.anova <- rbind(c(numDF, denDF, Fstat, Fstat2, pf(Fstat, 
            numDF, denDF, lower.tail = F)))
        colnames(ret.anova) <- c("numDF", "denDF", "F value", 
            "F2", "Pr(>F)")
        rownames(ret.anova) <- nam
        ret[[ii]]$anova <- ret.anova
        etahat <- L \%*\% beta
        etavar <- L \%*\% vc \%*\% t(L)
        etasd <- sqrt(diag(etavar))
        denDF <- apply(L, 1, function(x, dfs) min(dfs[x != 0]), 
            dfs = dfs)
        aod <- cbind(c(etahat), etasd, denDF, c(etahat/etasd), 
            2 * pt(-abs(etahat/etasd), denDF))
        colnames(aod) <- c("Estimate", "Std.Error", "DF", "t value", 
            "Pr(>|t|)")
        if (!is.null(clevel)) {
            hw <- qt(1 - (1 - clevel)/2, denDF) * etasd
            aod <- cbind(aod, LL = etahat - hw, UL = etahat + 
                hw)
            labs <- paste(c("Lower", "Upper"), format(clevel))
            colnames(aod)[ncol(aod) + c(-1, 0)] <- labs
        }
        rownames(aod) <- rownames(L)
        ret[[ii]]$estimate <- aod
        ret[[ii]]$vcov <- Vcov(fit, L)
        ret[[ii]]$vcor <- Vcor(fit, L)
        ret[[ii]]$L <- L
        ret[[ii]]$L.full <- L.full
        if (!is.null(heading)) 
            attr(ret[[ii]], "heading") <- heading
    }
    names(ret) <- names(Llist)
    attr(ret, "class") <- "glh"
    ret
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
