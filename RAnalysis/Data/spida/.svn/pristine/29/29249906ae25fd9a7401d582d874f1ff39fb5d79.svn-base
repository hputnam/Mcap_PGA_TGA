\name{write.sas}
\Rdversion{1.1}
\alias{write.sas}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Modified version of write.sas
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
write.sas(df, datafile = "G:/SAS/R2SAS.txt", codefile = "G:/SAS/R2SAS.sas")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{
%%     ~~Describe \code{df} here~~
}
  \item{datafile}{
%%     ~~Describe \code{datafile} here~~
}
  \item{codefile}{
%%     ~~Describe \code{codefile} here~~
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
function (df, datafile = "G:/SAS/R2SAS.txt", codefile = "G:/SAS/R2SAS.sas") 
{
    debug <- F
    pr <- function(x) {
        cat(deparse(substitute(x)), "\n")
        print(x)
        cat("==========\n")
        cat(x)
        cat("\n=============\n")
    }
    lrecl <- 256
    if (!debug) {
        write.table(df, file = datafile, row = FALSE, col = FALSE, 
            sep = ";", na = ".")
        lines <- scan(file = datafile, what = list("character"), 
            sep = "\n")
        lrecl <- max(nchar(lines))
    }
    nms <- names(df)
    nms.sas <- gsub("\\.", "_", nms)
    if (length(unique(toupper(nms.sas))) != length(nms.sas)) {
        ind <- duplicated(toupper(nms.sas))
        ind.rev <- duplicated(rev(toupper(nms.sas)))
        cat("Warning:\n")
        cat("The following R names may yield duplicate SAS names", 
            "\n", paste(nms[ind | rev(ind.rev)], collapse = " "), 
            "\n")
        warning("Possible duplicate SAS names")
    }
    factors <- sapply(df, is.factor) | sapply(df, is.character)
    classes <- sapply(df, class)
    odd.classes <- setdiff(sapply(df, class), c("numeric", "factor", 
        "character"))
    if (length(odd.classes) > 0) {
        cat("Warning:\n")
        cat("The following variables have classes that might not be handled properly by SAS\n")
        print(classes[grep(odd.classes, classes)])
        cat("\n")
    }
    factor.names <- nms[factors]
    factor.names.sas <- nms.sas[factors]
    dollarsign <- ifelse(factors, "$", "")
    factor.lengths <- sapply(df[factor.names], function(x) {
        if (is.factor(x)) 
            max(nchar(levels(x)))
        else max(nchar(x))
    })
    length.stmt <- paste(paste("   ", factor.names.sas, "$", 
        factor.lengths, "\n"), collapse = "")
    length.stmt <- paste("LENGTH\n", length.stmt, ";\n")
    if (debug) 
        pr(length.stmt)
    input.stmt <- paste(paste("    ", nms.sas, dollarsign, "\n"), 
        collapse = "")
    input.stmt <- paste("INPUT\n", input.stmt, ";\n")
    if (debug) 
        pr(input.stmt)
    code <- paste("filename r2sas '", datafile, "';\n", "libname to 'G:/SAS';\n", 
        "data to.r2sas;\n", "infile r2sas delimiter=';' dsd LRECL =", 
        lrecl + 100, ";\n", sep = "")
    code <- paste(code, length.stmt, input.stmt, "\nrun;\n")
    if (debug) 
        pr(code)
    if (!debug) 
        cat(code, file = codefile)
    invisible(0)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
