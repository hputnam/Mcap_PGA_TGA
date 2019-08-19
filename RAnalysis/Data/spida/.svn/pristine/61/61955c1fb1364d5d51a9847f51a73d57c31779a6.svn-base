\name{cap1}
\Rdversion{1.1}
\alias{cap1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Capitalize first character of each word
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cap1(x, tofactor = is.factor(x), stop = c(" The", " Of", " By", " To", " And"), blanks = c(" ", "(", "\"", "/", "+"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{tofactor}{
%%     ~~Describe \code{tofactor} here~~
}
  \item{stop}{
%%     ~~Describe \code{stop} here~~
}
  \item{blanks}{
%%     ~~Describe \code{blanks} here~~
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
function (x, tofactor = is.factor(x), stop = c(" The", " Of", 
    " By", " To", " And"), blanks = c(" ", "(", "\"", "/", "+")) 
{
    under2blank <- T
    if (is.factor(x)) {
        ret <- cap1(levels(x))
        if (length(unique(ret)) != length(ret)) 
            warning("factor levels have been shortened")
        levels(x) <- ret
        return(x)
    }
    ret <- as.character(x)
    for (ii in 1:length(ret)) {
        z <- ret[[ii]]
        if (under2blank) 
            z <- gsub("_", " ", z)
        n <- nchar(z)
        z <- substring(z, 1:n, 1:n)
        zu <- toupper(z)
        zl <- tolower(z)
        zb <- c(" ", zu[-n])
        z <- paste(ifelse(zb \%in\% blanks, zu, zl), collapse = "")
        for (ss in stop) z <- gsub(ss, tolower(ss), z)
        ret[[ii]] <- z
    }
    ret
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
