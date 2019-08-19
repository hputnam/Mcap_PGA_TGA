\name{xmerge}
\Rdversion{1.1}
\alias{xmerge}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extended merge with diagnostics
}
\description{
Extended merge with diagnostics.
This is a modification of \code{merge} that combines consistent variables
even if not specified in 'by' to keep a common name.
}
\usage{
xmerge(x, y, by, all = TRUE, dropdots = FALSE, verbose = FALSE, debug = TRUE, from = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, y}{
data frames, or objects to be coerced to one
}
  \item{by}{
%%     ~~Describe \code{by} here~~
}
  \item{all}{
%%     ~~Describe \code{all} here~~
}
  \item{dropdots}{
%%     ~~Describe \code{dropdots} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
  \item{debug}{
%%     ~~Describe \code{debug} here~~
}
  \item{from}{
%%     ~~Describe \code{from} here~~
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
%\references{
%%% ~put references to the literature/web site here ~
%}
\author{
Georges Monette
}
%\note{
%%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[base]{merge}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, y, by, all = T, dropdots = F, verbose = F, debug = T, 
    from = F, ...) 
{
    help <- "This is a modification of merge that combines consistent variables\neven if not specified in 'by' to keep a common name.\n-- Some errors fixed Apr 24, 2007"
    xm <- function(a, b, tofac = is.factor(a) || is.factor(b)) {
        if (tofac) {
            levs <- union(levels(b), levels(a))
            a <- as.character(a)
            b <- as.character(b)
        }
        b[is.na(b)] <- a[is.na(b)]
        if (tofac) {
            levs <- union(levs, unique(b))
            b <- factor(b, levels = levs)
        }
        b
    }
    na2f <- function(x) {
        x[is.na(x)] <- F
        x
    }
    consistent <- function(a, b) {
        if (is.factor(a)) 
            a <- as.character(a)
        if (is.factor(b)) 
            b <- as.character(b)
        !na2f(a != b)
    }
    if (from) {
        xname <- deparse(substitute(x))
        yname <- deparse(substitute(y))
        x[[".F"]] <- rep("x", nrow(x))
        y[[".F"]] <- rep("y", nrow(y))
    }
    xby <- x[, by, drop = F]
    yby <- y[, by, drop = F]
    xby$.file <- rep("x", nrow(xby))
    yby$.file <- rep("y", nrow(yby))
    by2 <- rbind(xby, yby)
    if (verbose) 
        cat("\nby in x and y:\n")
    if (verbose) 
        print(atotal(do.call("tab", by2), sum, "Total"))
    nams <- union(names(x), names(y))
    if (verbose) 
        print(c(DimX = dim(x), DimY = dim(y)))
    if (verbose) 
        cat("\nVariables in both:\n")
    if (verbose) 
        print(intersect(names(x), names(y)))
    if (verbose) 
        cat("\nVariables in X only:\n")
    if (verbose) 
        print(setdiff(names(x), names(y)))
    if (verbose) 
        cat("\nVariables in Y only:\n")
    if (verbose) 
        print(setdiff(names(y), names(x)))
    x$FromX <- 1:nrow(x)
    y$FromY <- 1:nrow(y)
    mm <- merge(x, y, by, all = T, ...)
    newroots <- setdiff(intersect(names(x), names(y)), by)
    if (verbose) 
        cat("\nDimension of merged data frames:\n")
    if (verbose) 
        print(c(DimMerge = dim(mm)))
    if (verbose) 
        cat("\nNames of variables in merged data frame:\n")
    if (verbose) 
        print(names(mm))
    if (F) {
        dotx <- grep("\\.x", names(mm), value = T)
        if (verbose) 
            print(c(dotx = dotx))
        doty <- grep("\\.y", names(mm), value = T)
        if (verbose) 
            print(c(doty = doty))
        rootx <- substring(dotx, 1, nchar(dotx) - 2)
        rooty <- substring(doty, 1, nchar(doty) - 2)
        newroots <- intersect(rootx, rooty)
    }
    FromBoth <- !is.na(mm$FromX) & !is.na(mm$FromY)
    Xonly <- !is.na(mm$FromX) & is.na(mm$FromY)
    Yonly <- is.na(mm$FromX) & !is.na(mm$FromY)
    if (verbose) 
        cat("\nRows in:\n")
    if (verbose) 
        print(c(Both = sum(FromBoth), Xonly = sum(Xonly), Yonly = sum(Yonly)))
    if (verbose) 
        cat("\nThe following variables occur in both data frames:\n")
    if (verbose) 
        print(newroots)
    drop.list <- character(0)
    for (nn in newroots) {
        nn.x <- paste(nn, ".x", sep = "")
        nn.y <- paste(nn, ".y", sep = "")
        mm[[nn]] <- xm(mm[[nn.x]], mm[[nn.y]])
        if (all(same <- consistent(mm[[nn.x]], mm[[nn.y]]))) {
            if (verbose) 
                cat("Variable ", nn, " is consistent\n")
            drop.list <- c(drop.list, nn)
        }
        else {
            if (verbose) 
                cat("Variable ", nn, " is inconsistent in the following rows:\n")
            if (verbose) 
                print(mm[same, c(by, nn.x, nn.y, nn)])
        }
    }
    if (dropdots) 
        drop.list <- newroots
    drop <- if (length(drop.list) > 0) {
        c(paste(drop.list, "x", sep = "."), paste(drop.list, 
            "y", sep = "."))
    }
    else character(0)
    if (verbose) 
        cat("\nDrop list:\n")
    if (verbose) 
        print(drop)
    if (length(drop) > 0) {
        if (verbose) 
            print(c(drop = drop))
        mm <- mm[, -match(drop, names(mm))]
    }
    onams <- 1:length(nams)
    onams <- c(onams, onams + 0.1, onams + 0.2)
    names(onams) <- c(nams, paste(nams, ".x", sep = ""), paste(nams, 
        ".y", sep = ""))
    keep <- intersect(names(sort(onams)), names(mm))
    mm[, keep]
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
