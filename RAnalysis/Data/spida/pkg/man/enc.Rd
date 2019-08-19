\name{enc}
\Rdversion{1.1}
\alias{enc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read coding tables in RDC
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
enc(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
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
function (x) 
{
    tonum <- function(x) as.numeric(gsub(",", "", as.character(x)))
    ff <- function(x) format(x, big.mark = ",")
    tran.table <- scan(what = list("character", "integer", "character", 
        "character"), flush = T)
    sam <- tonum(tran.table[[3]])
    pop <- tonum(tran.table[[4]])
    samp.pct <- 100 * (sam/pop)/(sum(sam, na.rm = T)/sum(pop, 
        na.rm = T))
    print(data.frame(Code = tran.table[[2]], Content = tran.table[[1]], 
        Sample = ff(sam), Popn = ff(pop), Sampling.Pct = round(samp.pct, 
            1)))
    tran(tran.table[[2]], tran.table[[1]], x, tofactor = T)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
