\name{Write.spss}
\Rdversion{1.1}
\alias{Write.spss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Write and SPSS file
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Write.spss(dataframe, file, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataframe}{
%%     ~~Describe \code{dataframe} here~~
}
  \item{file}{
%%     ~~Describe \code{file} here~~
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
function (dataframe, file, ...) 
{
    require(foreign)
    dname <- deparse(substitute(dataframe))
    disp(dname)
    cls <- sapply(dataframe, class)
    for (nn in names(dataframe)[cls == "Date"]) {
        dataframe[[nn]] <- as.character(dataframe[[nn]], "\%Y/\%m/\%d")
    }
    if (any(cls == "Date")) {
        cat("\nOpen .dta file in SPSS and convert following variables to dates\nwith yyyy/mm/dd format:\n")
        for (nn in names(dataframe)[cls == "Date"]) cat("       ", 
            nn, "\n")
    }
    for (nn in names(dataframe)[cls == "factor"]) {
        dataframe[[nn]] <- as.character(dataframe[[nn]])
    }
    if (missing(file)) 
        file <- paste(dname, ".dta", sep = "")
    else file <- sub("\\.dta$|\\.DTA$|$", ".dta", file)
    cat(paste("\nData saved in", file, "\n"))
    write.dta(dataframe, file, version = 7, ...)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
