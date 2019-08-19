\name{xqplot}
\Rdversion{1.1}
\alias{xqplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extended Quantile Plots
}
\description{
An easy way to see a dataset's variables at a glance. Show 'uniform quantile plot' for
numeical varibles and barcharts for factors. Quantile plots also show a horizontal line at the position of the mean
and at mean plus or minus one standard deviation.
}
\usage{
xqplot(x, ptype = "quantile", labels = dimnames(x)[[2]], ..., mfrow = findmfrow(ncol(x)), ask = prod(mfrow) < ncol(x) && dev.interactive(), mcex = 0.8, maxlab = 12, debug = F, mar = c(5, 2, 3, 1), text.cex.factor = 1, left.labs = F, maxvarnamelength = 20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{ptype}{
%%     ~~Describe \code{ptype} here~~
}
  \item{labels}{
%%     ~~Describe \code{labels} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
  \item{mfrow}{
%%     ~~Describe \code{mfrow} here~~
}
  \item{ask}{
%%     ~~Describe \code{ask} here~~
}
  \item{mcex}{
%%     ~~Describe \code{mcex} here~~
}
  \item{maxlab}{
%%     ~~Describe \code{maxlab} here~~
}
  \item{debug}{
%%     ~~Describe \code{debug} here~~
}
  \item{mar}{
%%     ~~Describe \code{mar} here~~
}
  \item{text.cex.factor}{
%%     ~~Describe \code{text.cex.factor} here~~
}
  \item{left.labs}{
%%     ~~Describe \code{left.labs} here~~
}
  \item{maxvarnamelength}{
%%     ~~Describe \code{maxvarnamelength} here~~
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
function (x, ptype = "quantile", labels = dimnames(x)[[2]], ..., 
    mfrow = findmfrow(ncol(x)), ask = prod(mfrow) < ncol(x) && 
        dev.interactive(), mcex = 0.8, maxlab = 12, debug = F, 
    mar = c(5, 2, 3, 1), text.cex.factor = 1, left.labs = F, 
    maxvarnamelength = 20) 
{
    help <- "\nxqplot                   coursefun.R     for PSYC 6140/MATH 6630 05/06\n\nExtended Quantile Plots and Barcharts\n\nDescription:\n\n     'xqplot' produces uniform quantile plots of numeric variables and\n     barcharts of factor variables.  The display is tuned to provide a\n     quick view of a data frame at a glance.\n\nUsage:\n     xqplot( x, mcex = 0.8, mar = c(5,2,3,1), text.cex.factor = 1,\n             mfrow,\n             maxlabs = 12,\n             left.labs = F)\n\nArguments:\n\n     x    : a data frame or list of variables to plot\n\n     mcex : character expansion factor for marginal text\n\n     mar  : size of margins\n\n     text.cex.factor : character expansion factor for barchart labels\n\n     left.labs : determines placement of barchart labels\n\n     maxlab : maximum number of categories to label in barcharts\n\n     mfrow :  number of rows and columns per page. If missing, an attempt\n            is made to choose a reasonable number.\n\n     maxvarnamelength : maximum length of variable name without splitting\n            on two lines.\nDetails:\n\nValue:\n\nBugs:\n\n     'mfrow' should take the total number of variables into account if they will\n     fill more than one page so the last page is close to being full.\n\n     The current version of the function could be made much simpler and\n     more transparent. Some code is redundant.\n\nReferences:\n\nContributed by:  G. Monette  2005-10-10\n\nModifications:\n\n\n"
    if (!is.list(x)) 
        x <- as.data.frame(x)
    left.labs <- rep(left.labs, length = length(x))
    findmfrow <- function(x) {
        if (x > 9) 
            c(3, 4)
        else cbind(`1` = c(1, 1), `2` = c(1, 2), `3` = c(2, 2), 
            `4` = c(2, 2), `5` = c(2, 3), `6` = c(2, 3), `7` = c(3, 
                3), `8` = c(3, 3), `9` = c(3, 3))[, x]
    }
    opt <- par(mfrow = mfrow, ask = ask, mar = mar + 0.1)
    on.exit(par(opt))
    if (debug) {
        cat("opt:\n")
        print(opt)
    }
    iscat <- function(x) is.factor(x) || is.character(x)
    Levels <- function(x) {
        if (is.factor(x)) 
            levels(x)
        else unique(x)
    }
    compute.cex <- function(x) {
        ll <- length(x)
        cex <- 2 * ifelse(ll < 5, 2, ifelse(ll < 10, 1, ifelse(ll < 
            20, 0.7, 0.5)))/mfrow[1]
    }
    for (ii in 1:dim(x)[2]) {
        vv <- x[[ii]]
        nam <- labels[[ii]]
        Nmiss <- sum(is.na(vv))
        N <- length(vv)
        if (iscat(vv)) {
            tt <- table(vv)
            xlab <- paste("N =", N)
            if (Nmiss > 0) {
                tt <- c(`<NA>` = sum(is.na(vv)), tt)
                xlab <- paste(xlab, "  Nmiss =", Nmiss)
            }
            ll <- names(tt)
            nn <- length(ll)
            if (left.labs[ii]) 
                barplot(tt, horiz = T, xlab = xlab, cex.names = text.cex.factor * 
                  compute.cex(nn))
            else {
                zm <- barplot(tt, names = rep("", nn), horiz = T, 
                  xlab = xlab)
                sel <- rep(T, length(tt))
                tt.sorted <- rev(sort(tt))
                if (nn > maxlab) 
                  sel <- tt > tt.sorted[maxlab]
                if (debug) {
                  disp(sel)
                  disp(nam)
                  disp(tt)
                  disp(tt.sorted)
                  disp(maxlab)
                  disp(tt.sorted[maxlab])
                  disp(sel)
                  disp(zm[sel])
                  disp(rep(max(tt), nn)[sel])
                  disp(ll[sel])
                }
                if (any(sel)) 
                  text(rep(max(tt), nn)[sel], zm[sel], ll[sel], 
                    adj = 1, cex = text.cex.factor * compute.cex(nn))
            }
        }
        else {
            sublab <- ""
            N <- length(vv)
            Ninfinite <- 0
            if (any(is.infinite(vv))) {
                n.pi <- sum(vv == Inf, na.rm = T)
                n.ni <- sum(vv == -Inf, na.rm = T)
                Ninfinite <- n.pi + n.ni
                vv <- vv[!is.infinite(vv)]
                sublab <- paste(sublab, "-Inf:", n.ni, "+Inf:", 
                  n.pi)
            }
            Nmiss <- 0
            if (any(is.na(vv))) {
                Nmiss <- sum(is.na(vv))
                vv <- vv[!is.na(vv)]
                sublab <- paste(sublab, "NA:", Nmiss)
            }
            Nok <- N - Nmiss - Ninfinite
            if (pmatch(ptype, "normal", nomatch = 0) == 1) {
                xxvar <- qnorm(ppoints(length(vv)))
                xlab <- paste("Normal quantile for", Nok, "obs.")
            }
            else {
                xxvar <- ppoints(length(vv))
                xlab <- paste("Fraction of", Nok, "obs.")
            }
            if (Nok == 0) {
                xxvar <- 1
                vv <- 1
                if (sublab == "") {
                  plot(xxvar, vv, xlab = xlab, ylab = "", type = "n")
                }
                else {
                  plot(xxvar, vv, xlab = xlab, ylab = "", type = "n", 
                    sub = sublab)
                }
                text(1, 1, "NA")
            }
            else {
                if (sublab == "") {
                  plot(xxvar, sort(vv), xlab = xlab, ylab = "Data", 
                    ...)
                }
                else {
                  plot(xxvar, sort(vv), xlab = xlab, ylab = "Data", 
                    ..., sub = sublab)
                }
                xm <- mean(vv)
                xs <- sqrt(var(vv))
                abline(h = xm, lty = 1)
                abline(h = c(xm - xs, xm + xs), lty = 2)
            }
        }
        vlab <- labels[ii]
        line.offset <- 1
        if (nchar(vlab) > maxvarnamelength) {
            vlab <- paste(substring(vlab, 1, maxvarnamelength), 
                "\n", substring(vlab, maxvarnamelength + 1))
            line.offset <- 0.2
        }
        mtext(vlab, 3, line.offset, cex = mcex)
    }
    if (debug) {
        disp(par())
    }
    invisible(0)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
