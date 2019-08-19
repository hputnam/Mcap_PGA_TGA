\name{abind}
\Rdversion{1.1}
\alias{abind}
\alias{abind.rdc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bind comformable arrays
}
\description{
Bind comformable arrays
}
\usage{
abind(arr1, arr2, d, facename = "")

bind.rdc(arr1, arr2, d, facename = "")

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{arr1}{
first array
}
  \item{arr2}{
second array
}
  \item{d}{
dimension along which arr1 and arr2 are joined, or, if the \code{length(dim(arr1)) == length(dim(arr2)) + 1}
the dimension in arr1 that is extended by arr2.
}
  \item{facename}{
Name for the new array dimension
}
}
\details{
\code{abind} binds two conformable arrays along a dimension.
\code{abind.rdc} is a version used in the RDC.

dim( arr1 ) and dim( arr2 ) must be equal except in the dth dimension.
If the length of dim( arr2 ) is 1 less than that of dim( arr1 ), then
'arr2' is treated as if it had dropped the dth dimension with size 1.
}
\value{
The returned value, ret, is an array with dimension dim(arr1)
except for the dth dimension where dim(ret)[d] == dim(arr1)[d] + dim(arr2)[d].

If length(dim(arr2)) == length(dim(arr1)) - 1, then arr2 is treated as if
it dropped a dimension of size 1 in the dth dimension. 'facename' is used
as the name of the dimnames list for this dimension.
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
\code{\link[base]{aperm}}, to permute arrays
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (arr1, arr2, d, facename = "") 
{
    d1 <- dim(arr1)
    n1 <- length(d1)
    dn1 <- dimnames(arr1)
    d2 <- dim(arr2)
    n2 <- length(d2)
    dn2 <- dimnames(arr2)
    arenull <- is.null(dn1) & is.null(dn2)
    if (is.null(dn1)) {
        dn1 <- lapply(as.list(d1), function(x) seq(1, x))
        dimnames(arr1) <- dn1
    }
    if (n1 != n2) {
        d2 <- d1
        d2[d] <- 1
        dn2 <- dn1
        dn2[[d]] <- facename
        dimnames(arr2) <- NULL
        dim(arr2) <- d2
        dimnames(arr2) <- dn2
        n2 <- n1
    }
    if (is.null(dn2)) {
        dn2 <- lapply(as.list(d2), function(x) seq(1, x))
        dimnames(arr2) <- dn2
    }
    perm <- 1:n1
    perm[c(d, n1)] <- c(n1, d)
    arr.p1 <- aperm(arr1, perm)
    arr.p2 <- aperm(arr2, perm)
    dret <- d1[perm]
    dret[n1] <- dret[n1] + d2[d]
    dnret <- dn1
    dnret[[d]] <- c(dnret[[d]], dn2[[d]])
    ret <- c(arr.p1, arr.p2)
    dim(ret) <- dret
    ret <- aperm(ret, perm)
    dimnames(ret) <- dnret
    ret
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
\keyword{array}% __ONLY ONE__ keyword per line
