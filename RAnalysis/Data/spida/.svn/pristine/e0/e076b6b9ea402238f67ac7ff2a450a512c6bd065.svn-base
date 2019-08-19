\name{sasin}
\Rdversion{1.1}
\alias{sasin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read a SAS ODS file
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sasin(file, tfile = tempfile())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{
%%     ~~Describe \code{file} here~~
}
  \item{tfile}{
%%     ~~Describe \code{tfile} here~~
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
function (file, tfile = tempfile()) 
{
    help = "\n    sasin reads a .csv file created by SAS with\n       ODS CSV FILE = 'file';\n        < SAS procedure statements >\n       ODS CSV CLOSE;\n    The tables produced by SAS are elements in the list\n    returned by sasin.\n    "
    todf <- function(ll) {
        if (length(ll) < 3) 
            return(character(0))
        if (length(ll) == 3) 
            return(ll[2])
        cat(ll[2], "\n", file = tfile)
        for (ii in 3:(length(ll) - 1)) {
            cat(ll[ii], "\n", file = tfile, append = T)
        }
        df <- read.csv(tfile, header = F)
        if (!any(sapply(df, is.numeric))) 
            df <- read.csv(tfile)
        df
    }
    readin <- scan(file, what = "", sep = "\n", blank.lines.skip = F)
    blanks <- which(readin == "")
    head.pos <- c(1, 1 + head(blanks, -1))
    heads <- gsub("\\\"|,", "", readin[head.pos])
    reps <- diff(c(head.pos, 1 + length(readin)))
    heads <- rep(heads, reps)
    readin <- split(readin, heads)
    readin <- lapply(readin, todf)
    readin
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
