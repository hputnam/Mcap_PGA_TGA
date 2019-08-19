# gicc
# Author: G. Monette  georges@yorku.ca
# Date: January 3, 2011
# Package: spida

gicc <- function(x, by, method = c('raw','na'), ...) {

  # Computes a 'generalized intra-class correlation coefficient' using
  # for factors: Goodman and Kruskal's tau based on proportional reduction in
  #    prediction error probability from a proportional prediction rule, and
  # for numerical variables: intraclass correlation coefficient(1,1) based
  #    on ranks of x
  #
  # Args:
  #   x: a variable or a data frame
  #   by: a variable to form cluster, or if x is a data frame, a formula
  #   method: not yet implemented
  # Returns:
  #   The gicc for a variable or for each variable in a data frame
  # TODO: implement different variation indices to see which discriminate
  #       most usefully (e.g. see p. 25 of Agresti (1990) Categorical Data Analysis, Wiley.
    UseMethod("gicc")
}


gicc.factor <- function( x, by, method = 'raw', ...) {
  # See generic function for description
   Verppr <- function( x ) {    # probability of error with proportional prediction rule
      px <- table( x ) / length(x)
      sum( px * (1-px))
   }
   if (method == 'na') x <- is.na(x)
   else {
     drop <- is.na(x) | is.na(by)
     x <- x[!drop]
     by <- by[!drop]
   }
   xc <- split(x,by)
   marg.perror <- Verppr(x)
   pcond <- sapply(xc,length)/length(x)
   cond.perror <- sum( pcond * sapply( xc, Verppr))
   #list( marginal = marg.perror, conditional = cond.perror, tau = (marg.perror - cond.perror)/marg.perror)
   (marg.perror - cond.perror)/marg.perror
}

gicc.default <- function(x, by, method = 'raw', ...) {
# intraclass correlation based on ranks to reduce influence of outliers
# perhaps we should use normal quantiles
   if ( method == 'na')
      return( gicc.factor( x, by, method ))
   if (is.character(x))
      return(gicc(as.factor(x), by))
   drop <- is.na(x) | is.na(by)
   x <- rank(x[!drop])
   by <- by[!drop]
   xc <- split(x,by)
   var.bet <- var( sapply( xc, mean))
   wt <- sapply(xc, length) -1
   wt <- wt / sum(wt)
   var.within <- sum(wt * sapply(xc, var), na.rm = TRUE)
   #disp(sapply(xc,var))
   #list(between = var.bet, within = var.within, icc = (var.bet)/(var.bet+var.within))
   var.bet/(var.bet+var.within)
}

gicc.data.frame <- function(x, by, method = 'raw', ...) {
# intraclass correlation based on ranks to reduce influence of outliers
# perhaps we should use normal quantiles
  if (inherits(by, "formula"))
    by <- model.frame(by, x, na.action = na.include)
    sapply( x, gicc, by, method)
}
