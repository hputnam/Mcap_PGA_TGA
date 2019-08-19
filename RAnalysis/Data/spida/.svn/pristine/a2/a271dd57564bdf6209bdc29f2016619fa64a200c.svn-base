\name{describe.vector}
\Rdversion{1.1}
\alias{describe.vector}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Describe a vector
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
describe.vector(x, descript, exclude.missing = TRUE, digits = 4, weights = NULL, normwt = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{descript}{
%%     ~~Describe \code{descript} here~~
}
  \item{exclude.missing}{
%%     ~~Describe \code{exclude.missing} here~~
}
  \item{digits}{
%%     ~~Describe \code{digits} here~~
}
  \item{weights}{
%%     ~~Describe \code{weights} here~~
}
  \item{normwt}{
%%     ~~Describe \code{normwt} here~~
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
function (x, descript, exclude.missing = TRUE, digits = 4, weights = NULL, 
    normwt = FALSE, ...) 
{
    oldopt <- options(digits = digits)
    on.exit(options(oldopt))
    if (length(weights) == 0) 
        weights <- rep(1, length(x))
    special.codes <- attr(x, "special.miss")$codes
    labx <- attr(x, "label")
    if (missing(descript)) 
        descript <- as.character(sys.call())[2]
    if (length(labx) && labx != descript) 
        descript <- paste(descript, ":", labx)
    un <- attr(x, "units")
    if (length(un) && un == "") 
        un <- NULL
    fmt <- attr(x, "format")
    if (length(fmt) && (is.function(fmt) || fmt == "")) 
        fmt <- NULL
    if (length(fmt) > 1) 
        fmt <- paste(as.character(fmt[[1]]), as.character(fmt[[2]]))
    present <- if (all(is.na(x))) 
        rep(FALSE, length(x))
    else if (is.character(x)) 
        (if (.R.) 
            x != "" & x != " " & !is.na(x)
        else x != "" & x != " ")
    else !is.na(x)
    present <- present & !is.na(weights)
    if (length(weights) != length(x)) 
        stop("length of weights must equal length of x")
    if (normwt) {
        weights <- sum(present) * weights/sum(weights[present])
        n <- sum(present)
    }
    else n <- round(sum(weights[present]), 2)
    if (exclude.missing && n == 0) 
        return(structure(NULL, class = "describe"))
    missing <- round(sum(weights[!present], na.rm = TRUE), 2)
    atx <- attributes(x)
    atx$names <- atx$dimnames <- atx$dim <- atx$special.miss <- NULL
    atx$class <- atx$class[atx$class != "special.miss"]
    isdot <- testDateTime(x, "either")
    isdat <- testDateTime(x, "both")
    x <- x[present, drop = FALSE]
    x.unique <- sort(unique(x))
    weights <- weights[present]
    n.unique <- length(x.unique)
    attributes(x) <- attributes(x.unique) <- atx
    isnum <- (is.numeric(x) || isdat) && !is.category(x)
    timeUsed <- isdat && testDateTime(x.unique, "timeVaries")
    z <- list(descript = descript, units = un, format = fmt)
    counts <- c(n, missing)
    lab <- c("n", "missing")
    if (length(special.codes)) {
        tabsc <- table(special.codes)
        counts <- c(counts, tabsc)
        lab <- c(lab, names(tabsc))
    }
    if (length(atx$imputed)) {
        counts <- c(counts, length(atx$imputed))
        lab <- c(lab, "imputed")
    }
    if (length(pd <- atx$partial.date)) {
        if ((nn <- length(pd$month)) > 0) {
            counts <- c(counts, nn)
            lab <- c(lab, "missing month")
        }
        if ((nn <- length(pd$day)) > 0) {
            counts <- c(counts, nn)
            lab <- c(lab, "missing day")
        }
        if ((nn <- length(pd$both)) > 0) {
            counts <- c(counts, nn)
            lab <- c(lab, "missing month,day")
        }
    }
    if (length(atx$substi.source)) {
        tabss <- table(atx$substi.source)
        counts <- c(counts, tabss)
        lab <- c(lab, names(tabss))
    }
    counts <- c(counts, n.unique)
    lab <- c(lab, "unique")
    x.binary <- n.unique == 2 && isnum && x.unique[1] == 0 && 
        x.unique[2] == 1
    if (x.binary) {
        counts <- c(counts, sum(weights[x == 1]))
        lab <- c(lab, "Sum")
    }
    if (isnum) {
        xnum <- if (.SV4.) 
            as.numeric(x)
        else oldUnclass(x)
        if (isdot) {
            dd <- sum(weights * xnum)/sum(weights)
            fval <- formatDateTime(dd, atx, !timeUsed)
            counts <- c(counts, fval)
        }
        else counts <- c(counts, format(sum(weights * x)/sum(weights), 
            ...))
        lab <- c(lab, "Mean")
    }
    if (n.unique >= 10 & isnum) {
        q <- if (any(weights != 1)) 
            wtd.quantile(xnum, weights, normwt = FALSE, na.rm = FALSE, 
                probs = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95))
        else quantile(xnum, c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 
            0.95), na.rm = FALSE)
        fval <- if (isdot) 
            formatDateTime(q, atx, !timeUsed)
        else format(q, ...)
        counts <- c(counts, fval)
        lab <- c(lab, ".05", ".10", ".25", ".50", ".75", ".90", 
            ".95")
    }
    names(counts) <- lab
    z$counts <- counts
    counts <- NULL
    if (n.unique >= 20) {
        if (isnum) {
            r <- range(xnum)
            xg <- pmin(1 + floor((100 * (xnum - r[1]))/(r[2] - 
                r[1])), 100)
            z$intervalFreq <- list(range = as.single(r), count = as.integer(tabulate(xg)))
        }
        lo <- x.unique[1:5]
        hi <- x.unique[(n.unique - 4):n.unique]
        fval <- if (isdot) 
            formatDateTime(c(oldUnclass(lo), oldUnclass(hi)), 
                atx, !timeUsed)
        else format(c(format(lo), format(hi)), ...)
        counts <- fval
        names(counts) <- c("L1", "L2", "L3", "L4", "L5", "H5", 
            "H4", "H3", "H2", "H1")
    }
    if (n.unique > 1 && n.unique < 20 && !x.binary) {
        tab <- wtd.table(if (isnum) 
            format(x)
        else x, weights, normwt = FALSE, na.rm = FALSE, type = "table")
        pct <- round(100 * tab/sum(tab))
        counts <- t(as.matrix(tab))
        counts <- rbind(counts, pct)
        dimnames(counts)[[1]] <- c("Frequency", "\%")
    }
    z$values <- counts
    structure(z, class = "describe")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
