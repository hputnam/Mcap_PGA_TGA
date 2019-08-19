###
###
### From spida-test
### Modified by GM 2009-09-26
### tab( ..., pct = 1) standardizes rows so thay add to 100
### tab( ..., pr = 1) standardizes rows so thay add to 1
### tab( ..., pr = 0) standardizes entire table so it sums to 1
###
###    Incorporate testing from package rdc
####



.mat2arr <- function(x) {
      ret <- as.list(x)
      dim(ret) <- NULL
      nams <- expand.grid( dimnames(x))
      for ( ii in 1:length(nams)) {
          nams[[ii]] <- paste( names(nams)[ii] , nams[[ii]], sep = " = ")
      }
      nams <- c( nams, sep = ", ")
      nams <- do.call( paste, nams)
      names(ret) <- nams
      ret
}
#mat2arr( zza)

Tab <- function(...) tab(...,total.margins = FALSE)
pab <- Tab     # legacy

tab <- function(x,...) UseMethod("tab")

tab.default <- function (..., total.margins = TRUE, pct = NULL, pr = NULL, useNA = "ifany",test=FALSE)
{
    #disp(pct)
    #disp(pr)
    aa <- list(...)
    if (length(aa) == 1 && is.list(aa[[1]])) {
        aa <- c(aa[[1]],list( total.margins = total.margins,  useNA = useNA, pr = pr, pct = pct, test=test))
        disp(aa[[1]])
        return(do.call("tab", aa[[1]]))
    }
    if (is.null(names(aa))) {
        nns = names(match.call())
        names(aa) = nns[2:(1 + length(aa))]
    }
    if(useNA=="ifany") for (ii in 1:length(aa)) aa[[ii]] <- factor(aa[[ii]], exclude = NULL)  # table uses 'ifany' correctly when a number but not for factors
    aa[["useNA"]] <- useNA
    # disp(aa)
    ret <- do.call("table", aa)
    if (test) {
          if( length(dim(ret)) < 3) {
              test.out <- chisq.test(ret)
          } else if ( length(dim(ret)) == 3)  {
              test.out <- apply( ret, 3:length(dim(ret)), chisq.test)

          } else if ( length(dim(ret)) > 3)  {
              test.out <- .mat2arr(apply( ret, 3:length(dim(ret)), chisq.test))

          }
    }
    if ( !is.null(pr)) ret <- acond( ret, MARGIN = pr)
    else if ( !is.null(pct)) ret <- 100* acond( ret, MARGIN = pct)
    else if (total.margins) ret = atotal(ret)
    if( test ) attr(ret,'test')<- test.out
    unclass(ret)
}


tab.data.frame <- function (dd, fmla, total.margins = TRUE, useNA = "ifany",  pct = NULL, pr = NULL, test = FALSE)
{
    if (missing(fmla))    {
        args <- c( as.list(dd), list(total.margins = total.margins,
            useNA = useNA, pct=pct, pr = pr,test=test))
        disp(args)
        return(do.call("tab",args))
    }
    xx = model.frame(fmla, dd, na.action = na.include)
    xx = c(xx, list( total.margins = total.margins, useNA = useNA, pct=pct, pr = pr, test=test) )
    do.call("tab", xx)
}

tab.formula <- function (fmla, dd, total.margins = TRUE, useNA = "ifany", ...)
{
    tab(dd, fmla, total.margins = total.margins, useNA = useNA,
        ...)
}

acond <- function (x, MARGIN = NULL, total.margins = TRUE)
{
    debug <- F
    x <- if( total.margins ) atotal(x)
    if (is.null(MARGIN)|| max(MARGIN) < 1)
        return(x/x[length(x)])
    d <- dim(x)
    if (debug) {
        cat("dim(x)\n")
        print(d)
    }
    dn <- dimnames(x)
    if (debug) {
        cat("dimnames(x)\n")
        print(dn)
    }
    n <- length(d)
    m <- length(MARGIN)
    if ( length(setdiff( MARGIN, 1:n)) > 0) stop("MARGIN must select dimensions of x, or be equal to 0")
    perm <- c(MARGIN, setdiff(1:n, MARGIN))
    perm.inverse <- order(perm)
    if (debug)
        disp(perm)
    if (debug)
        disp(perm.inverse)
    x <- aperm(x, perm)
    d <- dim(x)
    zl <- list(x)
    for (ii in 1:length(d)) zl[[ii + 1]] <- seq(if (ii <= m)
        1
    else d[ii], d[ii])
    tots <- do.call("[", zl)
    ret <- array(c(x)/c(tots), dim = d)
    ret <- aperm(ret, perm.inverse)
    if (debug)
        disp(dim(ret))
    if (debug)
        disp(dn)
    dimnames(ret) <- dn
    ret
}

aprop <- acond    # older name

