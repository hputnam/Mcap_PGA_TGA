\name{td}
\Rdversion{1.1}
\alias{td}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Set lattice parameters for multiple groups
}
\description{
Designed to easily sets lattice parameters for multiple groups. Setting parameters
before calling the lattice function allows parameters to be used consistently in the group key.
}
\usage{
td(new = FALSE, col = c("#0080ff", "#ff00ff", "darkgreen", "#ff0000", "orange", "#00ff00", "brown"), lty = 1:7, lwd = 1, pch = 1:7, cex = 0.8, font = 1, fill = "transparent", col.line = col, col.symbol = col, alpha = 1, alpha.line = alpha, alpha.symbol = alpha, len = 7, long = FALSE, record = FALSE, basecol = NULL, colsets = c("plot.symbol", "plot.line", "dot.symbol", "dot.line", "cloud.3d", "box.dot"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{new}{
If new = TRUE, open a new window, otherwise modify the existing active window, if there is one.
%%     ~~Describe \code{new} here~~
}
  \item{col}{
%%     ~~Describe \code{col} here~~
}
  \item{lty}{
%%     ~~Describe \code{lty} here~~
}
  \item{lwd}{
%%     ~~Describe \code{lwd} here~~
}
  \item{pch}{
%%     ~~Describe \code{pch} here~~
}
  \item{cex}{
%%     ~~Describe \code{cex} here~~
}
  \item{font}{
%%     ~~Describe \code{font} here~~
}
  \item{fill}{
%%     ~~Describe \code{fill} here~~
}
  \item{col.line}{
%%     ~~Describe \code{col.line} here~~
}
  \item{col.symbol}{
%%     ~~Describe \code{col.symbol} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{alpha.line}{
%%     ~~Describe \code{alpha.line} here~~
}
  \item{alpha.symbol}{
graphical parameters for superpose.line and superpose.symbol
%%     ~~Describe \code{alpha.symbol} here~~
}
  \item{len}{
extend the length of parameters by recycling to lenght 'len'
%%     ~~Describe \code{len} here~~
}
  \item{long}{
  if TRUE generate a default combination of col, lty and pch with length 42.
%%     ~~Describe \code{long} here~~
}
  \item{record}{
  If TRUE, set history to 'recording'. Caution: this can use a lot of memory in .GlobalEnv.
  Consider adding selected graphs to memory with the 'Insert' key instead.
%%     ~~Describe \code{record} here~~
}
  \item{basecol}{
%%     ~~Describe \code{basecol} here~~
}
  \item{colsets}{
%%     ~~Describe \code{colsets} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
'td' calls 'trellis.device' and sets graphical parameters for
'superpose.line' and 'superpose.symbol'. 'td' also initializes
a new trellis device with a white background if new = TRUE.
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
function (new = F, col = c("#0080ff", "#ff00ff", "darkgreen", 
    "#ff0000", "orange", "#00ff00", "brown"), lty = 1:7, lwd = 1, 
    pch = 1:7, cex = 0.8, font = 1, fill = "transparent", col.line = col, 
    col.symbol = col, alpha = 1, alpha.line = alpha, alpha.symbol = alpha, 
    len = 7, long = F, record = T, basecol = NULL, colsets = c("plot.symbol", 
        "plot.line", "dot.symbol", "dot.line", "cloud.3d", "box.dot"), 
    ...) 
{
    help <- "\ntd                   coursefun.R     for PSYC 6140/MATH 6630 05/06\n\nEasier way to call 'trellis.device'\n\nDescription:\n\n     'td' calls 'trellis.device' and sets graphical parameters for\n     'superpose.line' and 'superpose.symbol'. 'td' also initializes\n     a new trellis device with a white background by default.\n\nUsage:\n\n     td( new = F,\n         col=c(1,3,5,4,6,8,2),\n         col.line = col,\n         col.symbol = col,\n         alpha = 1,\n         alpha.line = alpha,\n         alpha.symbol = alpha,\n         lty=1:7,\n         lwd=1,\n\t pch = 1:7, cex = 0.8, font = 1,\n\t len = 7,\n\t long = F,\n         record = T,\n         basecol = NULL,\n\t colsets = c('plot.symbol','plot.line','dot.symbol',\n\t\t\t'dot.line','cloud.3d','box.dot'),...)\n\nArguments:\n\n     new : open a new device with: 'td(T)'.\n\n     col, lty, lwd, pch, cex: graphical parameters for superpose.line and superpose.symbol\n\n     len : extend the length of graphical parameters by recycling\n\n     long: if TRUE generate a default combination of col, lty and pch with length 42\n\n     record : initiate the history mechanism for the graphical device\n\n     colsets: additional graphical parameter lists to be modified\n\nDetails:\n\n     By using col and lty/pch with lengths that are relatively prime and\n     by using the len argument, one can generate unique combinations,\n     e.g. for len = 42 with col of length 6 and pch of length 7\n\nValue:\n\nReferences:\n\nContributed by:  G. Monette  2005-10-10\n\nModifications:\n\n"
    require(lattice)
    if (long) {
        col <- c(3, 5, 4, 6, 8, 2)
        len <- 42
    }
    if (new) 
        trellis.device(theme = col.whitebg, record = record, 
            new = new, ...)
    len <- max(len, length(col), length(lty), length(lwd), length(pch), 
        length(cex), length(font))
    spl <- trellis.par.get("superpose.line")
    spl$alpha <- rep(alpha.line, length = len)
    spl$lty <- rep(lty, length = len)
    spl$col <- rep(col.line, length = len)
    spl$lwd <- rep(lwd, length = len)
    trellis.par.set("superpose.line", spl)
    sps <- trellis.par.get("superpose.symbol")
    sps$alpha <- rep(alpha.symbol, length = len)
    sps$pch <- rep(pch, length = len)
    sps$col <- rep(col.symbol, length = len)
    sps$cex <- rep(cex, length = len)
    sps$font <- rep(font, length = len)
    sps$fill <- rep(fill, length = len)
    trellis.par.set("superpose.symbol", sps)
    list(superpose.symbol = sps, superpose.line = spl)
    if (!is.null(basecol)) {
        for (ii in colsets) {
            tt <- trellis.par.get(ii)
            tt$col <- basecol
            trellis.par.set(ii, tt)
        }
    }
    invisible(attr(.Device, "trellis.settings"))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
