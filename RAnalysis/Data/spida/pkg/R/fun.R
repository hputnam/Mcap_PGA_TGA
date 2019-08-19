

"
LOG OF CHANGES:
2010:
   Aug 13: removed obsolete version of Lmat in fun.R superceded by version in 
           wald.R  
           Ldiff moved to wald.R
           TODO: clean up other wald and eventually spline stuff
   Jul 21: renames Lag and LagI, cLag and cLagI respectively avoid conflicts
           with Hmisc::Lag.  Lag and and LagI are kept as aliases.
   Jun 16: Added non-ordered version of reorder.factor
   Jan 30: Added cvars

2009:
   May 8:  Added naresid to model.matrix.lme
   


'''fun.R:  A collection of utility functions and functions for multilevel modeling'''
'''        kept by Georges Monette.  These file can be downloaded or sourced in R'''
'''        from http://www.math.yorku.ca/~georges/R/fun.R'''

'''NOTE 1: Please send any suggestions or problems to georges@yorku.ca.

'''NOTE 2: There is a copy of this file at //wiki.math.yorku.ca/index.php/R:fun.R'''
'''        which can be edited. Changes will be incorporated regularly into the'''
'''        the downloadable file.'''

'''NOTE 3: THIS FILE USES QUOTES AND TAGS SO IT CAN BE RENDERED AS A WIKI FILE AS WELL AS R SOURCE'''

'''BE SURE TO LEAVE WIKI TEXT IN DOUBLE QUOTES AND AVOID DOUBLE QUOTES IN WIKI TEXT'''

Last uploaded to http://www.math.yorku.ca/~georges/R/fun.R : Auguest 23, 2006

:Decribe modifications here
:March 2009
:: panel.subgroups  allows different plotting 'types' within subgroups
:July 13, 2008
:: gsp: general spline program
:: smsp: smoothing spline using random effects model
:New: Aug 10, 2007
:: pchisq.mix: mixture of chisquares for different dfs for testing hypotheses
:::  concering null random effects in lme. Use simulate.lme to verify whether correct
:New: May 27, 2007
:: sasin      - read a SAS ODS CSV file and extract individual tables into a list
:: brace
:New: April 9, 2007
:: vif.lme    - variance inflation factors for lme fits
:: Rbind      - combine data and prediction data frame to plot
:::          results together
:New: November 9, 2006
:: oplot - plots number of observations that overplot
:New: August 23, 2006
:: cvar( x, id ) ; dvar( x, id )
::: cvar creates a contextual variable that is the group mean of 'x' within each level of the factor 'id'. If 'x' is a factor, 'cvar' returns a suitably labelled matrix that is the group mean of the coding variables for the factor 'x'.
::: dvar is x - cvar( x, id ) where x is turned into its coding matrix if x is a factor.
:New: July 10, 2006
:: xmerge( x , y , by )
::: merge that tests consistency of common names not in the by argument. Values are combined with priority given to 'y'. If variables are inconsistent, the '.x' and '.y' versions are also left in the merged data frame.
::: merge two data.frames with diagnostics
:: long(data, varying)
::: reshapes data from wide to long format using variable names given in list 'varying'. Vectors in the list are old names, names of vectors are new names. Names alone serve are roots. 
:Modified: June 10, 2006
:: added:
::: up(dd, form)  - creates a higher level data set with one row per case defined by 'form'
:Modified:  --[[User:Georges|Georges Monette]] 14:30, 28 May 2006 (EST)
:: added functions from funRDC
:Modified:  --[[User:Georges|Georges Monette]] 05:56, 22 Feb 2006 (EST)
:: added Lmat
:Modified:  --[[User:Georges|Georges Monette]] 09:45, 2 Nov 2005 (EST)
:: added class 'cat' to coursefun
:: defined print.cat to print with cat
:Modified:
:: plot3d, identify3d by John Fox so they work with matrices, data.frames or variable arguments
:: ell  - corrected error when using radius
::



== General description ==
<pre>
"
##
##
##  Some R functions for PSYC 6140 and MATH 6630
##  2005-2006
##
## Last update: October 27, 2005
## Summary:
##
##
## Tables:
##
##    atotal: border an array with sums
##    abind : glue two arrays together on a selected dimension
##
## General Linear Hypothesis
##
#     glh   : glh( fit, L )
##    Lmat  : generates L matrix for anova in library(lmer) or lht in library(car)
##
## Graphics:
##
##    td    : easy front end to trellis.device and related functions
##    xqplot: extended quantile plots
##
## 3D graphics by John Fox:
##
##    scatter3d
##    identify3d
##    ellipsoid
##    plot3d     - wrapper for scatter3d
##
## Inference
##    cell   - a modified version of car::confidence.ellipse.lm that
##             creates a confidence ellipse to plot with lines(...)
##    dell   - data ellipse
##
##            
# Note that some functions that were transferred and improved in coursefun.R
# have been 'disabled' by appending '.rdc' to function name
#
#  Splus functions written in RDC
#  2005:
#  April 19       Modified for R
#  May 3          copied atotal, abind
#                 wrote acond
#  May 10         new function: cap1 to capitalize initial letters and turn underscores to blanks
#  May 12         new function adapted from gm: td 
#  May 13         new function: write.sas to write data frame to SAS without truncating variable names
#  May 16         added Poly and centering to splines
#  June 13        added constant to check if values are constant within levels of id
#                 added varLevel( data.frame, ~lev1/lev2) to report level of each variable
#                 added Lmat to generate L matrix with 1's for effects containing a string
#                 modified anova.lme to use min DFs in 'denominator' DFs
#  August 3       modified Lag to accept 'at' argument
#  August 5       changed 'arep' to 'apct' in order to parallel atotal and acond
#                 changed acond to aprop
#  August 15      getFix, glh and and print.glh, Q, Vcov, Vcor
#  August 25      Contrasts
#  2006
#  May 19         fill and capply
#  October 2      cvar:  create contextual variable


##
##  Crude predict methods for 'mer'  Dec 6, 2008
##

predict.mer <- function( model, data = model.matrix(model), form , verbose = FALSE) {

help   = "
    This is a crude predict method for class 'mer'
    Limitations:
    1) When invoked on the model it only returns the linear predictor for the
       complete data. Note that a complete data frame is easily obtained with
       model.frame( model )
    2) When the 'data' argument is provided for 'new' data, it is also
       necessary to provide the correct rhs of the fixed effects formula.

    "
        if (missing(form) & ! missing(data)) {
            cat( help)
            stop( "Need 'form' if 'data' given")
        }
        if( ! missing( data )){
            data = model.matrix(form, data)
            cnames = colnames(data)
            if( verbose ) print( cnames )
            fnames = names( fixef( model ))
            if (verbose) print( fnames)
            if ( any( cnames != fnames)) {
                cat("\nMatrix names:\n")
                print( cnames )
                cat("\nCoeff names:\n")
                print( fnames)
                warning("matrix and coeff names not the same")
            }
        }
        data %*% fixef( model)

}




##
## Linear algebra
##


ConjComp <- function( X , Z = diag( nrow(X)) , ip = diag( nrow(X)), tol = 1e-07 ) {
    help <- "
    ConjComp returns a basis for the conjugate complement of the
    conjugate projection of X into span(Z) with respect to inner product with
    matrix ip.
    Note: Z is assumed to be of full column rank but not necessarily X.
    "
    
    xq <- qr(t(Z) %*% ip %*% X, tol = tol)
    if ( xq$rank == 0 ) return( Z )
    a <- qr.Q( xq, complete = TRUE ) [ ,-(1:xq$rank)]
    Z %*% a
}

OrthoComp <- function (X, Z , tol = 1e-07) ConjComp( X, Z, tol = tol)

oplot <- function( x, y, ..., verbose = TRUE) {
    pat <- paste( x, y, sep = ",")
    keep <- !duplicated(pat)
    
    ns <- table(pat)
    ns <- ns[pat[keep]]      # to order ns so it matches pat[keep]
    nps <- as.character(ns)
    #nps [ns>9] <- "*"
    x <- x[keep]
    y <- y[keep]
    
    
    if (verbose) {
        print(pat)
        print(table(pat))
        print(keep)
        print(ns)
    }
    plot( x, y,pch = "o",cex = 5,...)
    text( x, y, nps)
    
}
#plot( c(1,2,3,2,1,2,3), c(1,2,3,2,1,2,3), pch = 'AB')
#oplot( c(rep(1,10),2,3,2,1,2,3), c(rep(1,10),2,3,2,1,2,3), type = 'b')

  
fac <- function(x) {
   xx <- svd(x)
   ret <- t(xx$v) * sqrt(pmax( xx$d,0))
   ret  #ret [ nrow(ret):1,]
}



disp <- function( x , head = deparse(substitute(x))) {
    # for debugging
    cat("\n::: ", head , " :::\n")
    print(x)
}


"
</pre>

== General Linear Hypothesis ==
<pre>
"


eg <-
function( df, by, ...) {
    # a quicker version of expand.grid
    # should work with fits
        dots = list(...)
        by = model.frame( by, df, na.action = na.include)
        byn = names(by)  # will this stay
        names(byn) = byn
        args = lapply( byn, function(x){
                vv = df[[x]]
                if ( is.factor(vv) ) levels(vv) else unique(vv)
        })
        args = c( args,dots)
        do.call('expand.grid', args)
}

##
##
##   Functions to perform a GLH on lm, lme or lmer models
##   August 13, 2005
##
##
##
##   Lmat: generate a hypothesis matrix based on a pattern
##
##   glh
##   Lmat
##   Ldiff
##   getFix
##
##   print.glh
##


getFix <- function(fit,...) UseMethod("getFix")

getFix.lm <- function(fit,...) {
       ss <- summary(fit)
       ret <- list()
       ret$fixed <- coef(fit)
       ret$vcov <- ss$sigma^2 * ss$cov.unscaled
       ret$df <- rep(ss$df[2], length(ret$fixed))
       ret
}


getFix.glm <- function(fit,...) {
       ss <- summary(fit)
       ret <- list()
       ret$fixed <- coef(fit)
       ret$vcov <- vcov(fit)
       ret$df <- rep(ss$df.residual, length(ret$fixed))
       ret
}

getFix.lme <- function(fit,...) {
       require(nlme)
       ret <- list()
       ret$fixed <- nlme::fixef(fit)
       ret$vcov <- fit$varFix
       ret$df <- fit$fixDF$X
       ret
}

getFix.lmer <- function(fit,...) {
       ret <- list()
       ret$fixed <- fixef(fit)
       ret$vcov <- as.matrix( vcov(fit) )
       # ret$df <- Matrix:::getFixDF(fit)
       ret$df <- rep( Inf, length(ret$fixed))
       ret
}

getFix.glmer <- function(fit,...) {
       ret <- list()
       ret$fixed <- fixef(fit)
       ret$vcov <- as.matrix(vcov(fit))
       # ret$df <- Matrix:::getFixDF(fit)
       ret$df <- rep( Inf, length(ret$fixed))
       ret
}


getFix.mer <- function(fit,...) {
       ret <- list()
       ret$fixed <- fixef(fit)
       ret$vcov <- as.matrix(vcov(fit))
       # ret$df <- Matrix:::getFixDF(fit)
       ret$df <- rep( Inf, length(ret$fixed))
       ret
}

Vcov <- function(fit) {
     getFix(fit)$vcov
}

Vcor <- function(fit) {
     vc <- cov2cor(getFix(fit)$vcov)
     svds <- svd(vc)$d
     attribute(vc,'conditionNumber') <- svds[1]/svds[length(svds)]
     vc
}


na2f <- function(x)  {
    x[is.na(x)] <- FALSE
    x
}


# fit <- lmer( Yield ~ Location *   Family  + (1|Block), data = Genetics)
# getFix(fit)
# fit <- lme( Yield ~ Location *   Family  , data = Genetics, random = ~1|Block)

# L <- rbind( c(1,1,1,1), c(0,1,0,0), c(1,0,1,1))


tfun <- function( x) {
      ret <- NULL
      if ( is.character(x)) ret <- "it's character"
      else { 
         if ( is.numeric(x)) {
            if ( is.null(dim(x))) {
                if ( length(x) != 4 ) ret <- diag(4)[x,]
                else ret <- rbind( x )
            }
         }
      }
      ret 
}

##
##  Extension of avp from car
##

av.frame <- function( model, ..., help = FALSE) {
if(help) {
 cat("
       av.frame( model, variable)
       returns a data frame with model.frame(model) augmented
       by y.res and x.res, the residuals for an added variable
       plot
       
       The purpose of this function is to facilitate OLS av.plots
       for mixed models.



       Example:

       library(nlme)
       library(lattice)
       hs <- read.csv( 'http://www.math.yorku.ca/~georges/Data/hs.csv')
       
       # Mixed model where ses and Sex are Level 1 and Sector is Level 2
       
       fit.mm <- lme( mathach ~ ses * Sex * Sector, hs, random = ~ 1+ses| school)

       # for diagnostics fit an OLS model using only level 1 variables interacting
       # with the id variable

       fit.ols <- lm( mathach ~ (ses * Sex ) * factor(school), hs)
       xyplot( y.res ~ x.res | factor(school), cbind(av.frame(fit.ols, 'ses:Sex'),hs), sub = 'ses:Sex')
       xyplot( y.res ~ x.res | factor(school), cbind(av.frame(fit.ols, '^Sex'),hs), sub = 'Sex')
       xyplot( y.res ~ x.res | factor(school), cbind(av.frame(fit.ols, '^ses$|^ses:f'),hs), sub = 'ses')

       
       Note : y.res is the residual from fitting the response on
              the model matrix for fit.ols omitting any column
              whose names is matched (as a regular expression)
              by 'effect'
              x.res is the residual of the first column of the
              model matrix that is matched by 'effect' on the
              same matrix used for y.res.
       Caution: To make sure that the correct columns were
              matched, the list of matched columns that are omitted
              is printed.

")
  return( invisible(0))
}
         UseMethod("av.frame")
}

av.frame.lm <- function (model, variable,...){
# code borrowed from 'car' by J. Fox, function 'av.plot'
# labels = names(residuals(model)[!is.na(residuals(model))]),
#    identify.points = TRUE, las = par("las"), col = palette()[2],
#    pch = 1, lwd = 2, main = "Added-Variable Plot", ...)


    mod.mat <- model.matrix(model)
    var.names <- colnames(mod.mat)
    omit <- grep( variable, var.names)
    if (0 == length(omit))
        stop(paste(variable, "is not matched among columns of the model matrix."))

    cat( "x.var =", var.names[ omit[1] ], "\n",
    "omitted vars =", var.names[omit[-1]], "\n")
    
    response <- response(model)
    x.var <- mod.mat[,omit[1]]
    Xpred <- mod.mat[, - omit ]
    preds <- predict( update( model, na.action = na.exclude))

    responseName <- responseName(model)
    if (is.null(weights(model)))
        wt <- rep(1, length(response))
    else wt <- weights(model)
    res <- lsfit(mod.mat[, -omit], cbind(mod.mat[, omit[1]], response),
        wt = wt, intercept = FALSE)$residuals
    ret <- matrix(NA, nrow = length( preds), ncol = 2)
    ret[ !is.na(preds),] <- res
    data.frame( x.res = ret[,1], y.res = ret[,2])
    
}



"vif.lme" <-
function (mod)
{
    if (any(is.na(fixef(mod))))
        stop("there are aliased coefficients in the model")
    v <- vcov(mod)  # vcov.lme is in library(stats)
    mm <- model.matrix( formula(mod), mod$data)
    assign <- attributes(mm)$assign
    if (names(fixef(mod)[1]) == "(Intercept)") {
        v <- v[-1, -1]
        assign <- assign[-1]
    }
    else warning("No intercept: vifs may not be sensible.")
    terms <- labels(terms(mod))
    n.terms <- length(terms)
    if (n.terms < 2)
        stop("model contains fewer than 2 terms")
    R <- cov2cor(v)
    detR <- det(R)
    result <- matrix(0, n.terms, 3)
    rownames(result) <- terms
    colnames(result) <- c("GVIF", "Df", "GVIF^(1/2Df)")
    for (term in 1:n.terms) {
        subs <- which(assign == term)
        result[term, 1] <- det(as.matrix(R[subs, subs])) * det(as.matrix(R[-subs,
            -subs]))/detR
        result[term, 2] <- length(subs)
    }
    if (all(result[, 2] == 1))
        result <- result[, 1]
    else result[, 3] <- result[, 1]^(1/(2 * result[, 2]))
    result
}







"</pre>

== Trellis graphics ==
<pre>
"




###
###  td
###



td <- function(
    new = FALSE,
    col=c("#0080ff",   "#ff00ff",   "darkgreen", "#ff0000" ,  "orange" ,   "#00ff00",   "brown" ),
    lty=1:7, lwd=1,
	pch = 1:7, cex = 0.8, font = 1,
	fill = "transparent",
	col.line = col,
	col.symbol = col,
	alpha = 1,
	alpha.line = alpha,
	alpha.symbol = alpha,
	len = 7,
	long = FALSE,
    record = FALSE,
    basecol = NULL,
	colsets = c('plot.symbol','plot.line','dot.symbol',
			'dot.line','cloud.3d','box.dot'),...) {

help <- "
td                   coursefun.R     for PSYC 6140/MATH 6630 05/06

Easier way to call 'trellis.device'

Description:

     'td' calls 'trellis.device' and sets graphical parameters for
     'superpose.line' and 'superpose.symbol'. 'td' also initializes
     a new trellis device with a white background by default.

Usage:

     td( new = F,
         col=c(1,3,5,4,6,8,2),
         col.line = col,
         col.symbol = col,
         alpha = 1,
         alpha.line = alpha,
         alpha.symbol = alpha,
         lty=1:7,
         lwd=1,
	 pch = 1:7, cex = 0.8, font = 1,
	 len = 7,
	 long = F,
         record = FALSE,
         basecol = NULL,
	 colsets = c('plot.symbol','plot.line','dot.symbol',
			'dot.line','cloud.3d','box.dot'),...)

Arguments:

     new : open a new device with: 'td(T)'.

     col, lty, lwd, pch, cex: graphical parameters for superpose.line and superpose.symbol

     len : extend the length of graphical parameters by recycling

     long: if TRUE generate a default combination of col, lty and pch with length 42

     record : initiate the history mechanism for the graphical device

     colsets: additional graphical parameter lists to be modified

Details:

     By using col and lty/pch with lengths that are relatively prime and
     by using the len argument, one can generate unique combinations,
     e.g. for len = 42 with col of length 6 and pch of length 7

Value:

References:

Contributed by:  G. Monette  2005-10-10

Modifications:

"

        # Modified for R: Oct. 10, 2004
	#
	# reset superpose.symbol and superpose.line so they have consistent length
	# equal to the max of input parameters:
	#			sps: cex, col, font, pch
	#			spl: lty, col, lwd
	# or to len
	#  This allows distinctive line styles, for example for 12 objects by using
	#  good lty's:     1,4,8,6,3,5,7
	#  and good col's: 1,8,6,3,4,5,2
	#

  require(lattice)
  if ( long ) {
    col <- c(3,5,4,6,8,2)   # drop yellow
    len <- 42 # generates 42 unique combinations of pch/lty and col
  }
  if(new) trellis.device(theme = col.whitebg, record = record, new = new, ...)
                                        # NOTE: fixed panel.superpose so lty argument
                                        # is passed to points for type = 'b'
  len <- max(len,length(col),length(lty),length(lwd),length(pch),length(cex),length(font))
  spl <- trellis.par.get("superpose.line")
  spl$alpha <- rep(alpha.line, length = len)
  spl$lty <- rep(lty,length=len)
  spl$col <- rep(col.line,length=len)
  spl$lwd <- rep(lwd,length=len)
  trellis.par.set("superpose.line",spl)
  sps <- trellis.par.get("superpose.symbol")
  sps$alpha <- rep( alpha.symbol, length = len)
  sps$pch <- rep(pch, length=len)
  sps$col <- rep(col.symbol, length=len)
  sps$cex <- rep(cex, length=len)
  sps$font <- rep(font, length=len)
  sps$fill <- rep(fill, length=len)

  trellis.par.set("superpose.symbol",sps)
  list(superpose.symbol = sps, superpose.line = spl)
  if ( !is.null(basecol)) {
    for ( ii in colsets ) {
      tt <- trellis.par.get(ii)
      tt$col <- basecol
      trellis.par.set(ii,tt)
    }
  }
  ret <- trellis.par.get()
  invisible(ret[grep('superpose',names(ret))])
}


"</pre>

== Quantile plots ==
<pre>"



###
###  xqplot
###


xqplot <- function(x, 
                    ptype = "quantile", 
                    labels = dimnames(x)[[2]], ..., 
                    mfrow = findmfrow ( ncol(x)),
                    ask = prod(mfrow) < 
                            ncol(x) && dev.interactive(),
                   
                    mcex = 0.8, maxlab = 12 , debug = F,
                    mar = c(5,2,3,1),
                    text.cex.factor = 1 ,
                    left.labs = F,
                    maxvarnamelength = 20)
{
help <- "
xqplot                   coursefun.R     for PSYC 6140/MATH 6630 05/06

Extended Quantile Plots and Barcharts

Description:

     'xqplot' produces uniform quantile plots of numeric variables and
     barcharts of factor variables.  The display is tuned to provide a
     quick view of a data frame at a glance.

Usage:
     xqplot( x, mcex = 0.8, mar = c(5,2,3,1), text.cex.factor = 1,
             mfrow,
             maxlabs = 12,
             left.labs = F)

Arguments:

     x    : a data frame or list of variables to plot

     mcex : character expansion factor for marginal text

     mar  : size of margins

     text.cex.factor : character expansion factor for barchart labels

     left.labs : determines placement of barchart labels

     maxlab : maximum number of categories to label in barcharts

     mfrow :  number of rows and columns per page. If missing, an attempt
            is made to choose a reasonable number.

     maxvarnamelength : maximum length of variable name without splitting
            on two lines.
Details:

Value:

Bugs:

     'mfrow' should take the total number of variables into account if they will
     fill more than one page so the last page is close to being full.

     The current version of the function could be made much simpler and
     more transparent. Some code is redundant.

References:

Contributed by:  G. Monette  2005-10-10

Modifications:


"
  ## Adapted from myplot.data.frame for R by G. Monette, Oct. 25, 2004
  ##    maxlab is maximum number of labels
  # Turn matrices into variables:
    if (! is.list(x)) x <- as.data.frame(x)
  
  if ( any ( sapply( x, class) == 'matrix') ) {
       zz <- list()
       for ( ii in seq_along( x )) {
           if ( is.matrix( x[[ii]])) {
                  if ( is.null (colnames( x[[ii]]))) {
                        cnames <- paste( names(x)[ii], 1:ncol(x[[ii]]), sep ='.')
                  } else {
                         cnames <- paste( names(x)[ii], colnames(x[[ii]]), sep = '.')
                  }
                  for ( jj in seq_len( ncol ( x[[ii]]))) {
                       zz[[cnames[jj] ]] <- x[[ii]][,jj]
                  }

           } else {
               zz[[ names(x)[[ii]] ]] <- x[[ii]]
           }
          }
          x <- as.data.frame(zz)
          #disp( x )
  }
  

  left.labs <- rep( left.labs, length = length(x))
  findmfrow <- function( x ) {
	   if ( x > 9) c(3,4)
	   else cbind( '1'=c(1,1),'2'=c(1,2),'3'=c(2,2),
                   '4'=c(2,2),'5'=c(2,3),'6'=c(2,3),
                   '7'=c(3,3), '8'=c(3,3), '9'=c(3,3)) [, x]
  }

  opt <- par( mfrow = mfrow, ask = ask , mar = mar + 0.1 )
  on.exit(par(opt))
  if(debug) { cat("opt:\n");print(opt)}

  iscat <- function( x ) is.factor(x) || is.character(x)

  Levels <- function(x) {
      if ( is.factor(x)) levels(x) else unique(x)
  }


  compute.cex <- function( x ) {
    ll <- length(x)
    cex <- 2 * ifelse( ll < 5, 2,
                      ifelse( ll < 10, 1,
                             ifelse( ll < 20, .7, .5)))/mfrow[1]
  }
  for ( ii in 1: dim(x)[2]) {
    vv <- x[[ii]]
    nam <- labels[[ii]]
    Nmiss <- sum(is.na(vv))
    N <- length(vv)
    if ( iscat(vv) ){
      tt <- table(vv)

      xlab <- paste("N =", N )
      if ( Nmiss > 0 ) {
        tt <- c( "<NA>" = sum(is.na(vv)), tt)
        xlab <- paste(xlab, "  Nmiss =", Nmiss)
      }
      ll <- names(tt)
      nn <- length(ll)
      if ( left.labs[ii] ) barplot( tt, horiz = TRUE,
                                   xlab = xlab,
                                   cex.names = text.cex.factor * compute.cex(nn) )
      else {
        zm <- barplot( tt, names = rep("",nn), horiz = TRUE, xlab = xlab)
        ## If nn > maxlab drop labels for smaller frequencies
        sel <- rep( T, length(tt))
        tt.sorted <- rev(sort(tt))
        if ( nn > maxlab ) sel <- tt > tt.sorted[maxlab]
        if (debug) {
          disp(sel)
          disp(nam)
          disp(tt)
          disp(tt.sorted)
          disp(maxlab)
          disp(tt.sorted[maxlab])
          disp(sel)
          disp(zm[sel])
          disp(rep(max(tt),nn)[sel])
          disp( ll[sel])
        }
        if ( any(sel) ) text( rep( max( tt ), nn)[sel]  ,
                             zm[sel], ll[sel], adj = 1, cex = text.cex.factor * compute.cex( nn ))
      }
    } # end of iscat(vv)
    else {
      sublab <- ""
      N <- length( vv )
      Ninfinite <- 0
      if ( any( is.infinite ( vv ) ) ){
            n.pi <- sum( vv == Inf , na.rm = TRUE)
            n.ni <- sum( vv == -Inf, na.rm = TRUE )
            Ninfinite <- n.pi + n.ni
            vv <- vv[!is.infinite(vv)]
            sublab <- paste( sublab,"-Inf:",n.ni,"+Inf:",n.pi)
      }
      Nmiss <- 0
      if ( any ( is.na( vv )  )) {
            Nmiss <- sum( is.na(vv) )
            vv  <- vv[!is.na(vv)]
            sublab <- paste( sublab, "NA:", Nmiss)
      }
      Nok <- N - Nmiss - Ninfinite
      if ( pmatch( ptype, 'normal', nomatch = 0) == 1 ) {
            xxvar <- qnorm( ppoints(length(vv)) )
            xlab <- paste("Normal quantile for", Nok, "obs.")
      }
      else {
          xxvar <- ppoints( length(vv) )
          xlab <- paste("Fraction of", Nok, "obs.")
      }

      ## Plot continuous
      if ( Nok == 0 ) {
        xxvar <- 1
        vv <- 1
        if ( sublab == "") {
            plot( xxvar, vv, xlab = xlab, ylab="", type = 'n')
        } else {
            plot( xxvar, vv, xlab = xlab, ylab="", type = 'n', sub = sublab)
        }
        text( 1, 1, "NA")
      }
      else {
        if ( sublab == "") {
            plot(xxvar, sort(vv), xlab = xlab, ylab = "Data", ...)
        } else {
            plot(xxvar, sort(vv), xlab = xlab, ylab = "Data", ..., sub = sublab)
        }
        xm <- mean(vv)
        xs <- sqrt(var(vv))
        abline( h= xm,lty=1)
        abline( h= c(xm-xs,xm+xs),lty=2)
      }
    }
    ## Titles for all plots
    vlab <- labels[ii]
    line.offset <- 1.0
    if ( nchar( vlab ) > maxvarnamelength) {
        vlab <- paste( substring(vlab,1,maxvarnamelength), "\n",substring(vlab, maxvarnamelength + 1))
        line.offset <- 0.2
    }
    mtext(vlab, 3, line.offset , cex = mcex)
  }
  # par(opt)
  if(debug) { disp(par()) }
  invisible(0)
}




sampler <-
    function( n=24 ) {
    # sample of lines and symbols
     old.par <- par(ask=T)
     on.exit( par(old.par))
      require(lattice)

     y <- 0:n
     x <- 0:n
     
     print(xyplot( y ~ x, type = 'n', xlab = 'lty', ylab = 'col',
      panel = function(x,y,...) {
      for ( i in x) {
       panel.xyplot(c(i,i),range(y),type='l',lty=i,col=1,lwd = 3)
      }
      for ( i in y) {
       for ( j in seq(0,.9, by = .1)) {
        panel.xyplot(c(min(x)+ j*(max(x)-min(x)),min(x)+ (j+.1)*(max(x)-min(x))),c(i,i),type='l',lty=1,col=i, lwd = 3)
       }
      }
     }))

     # print(z$x, z$y, ylim=c(0,7))
     spl <-trellis.par.get('superpose.line')
     z <- expand.grid( y = 1:length(spl$lty), x = 0:2)
     print(xyplot( y ~ x , z, ylim =c(0,length(spl$lty)),groups = y, type='b',
            main="superpose.line and .symbol"))

     y <- 10*(0:25)
     x <- 0:9
     print(xyplot( y ~ x, type = 'n', main = 'pch',
        xlab = expression( ~ alpha + beta + gamma + delta[epsilon] + zeta^eta + theta + iota+kappa),
        ylab = expression( ~ lambda + mu + nu + xi + omicron + pi + rho + sigma + tau + upsilon + phi + chi +psi + omega),
      panel = function(x,y,...) {
      for ( i in x) {
       for ( j in y ) {
        panel.xyplot(i,j,pch=i+j,cex = 2)
       }
      }
     }))

     invisible(0)
}

pal <- function(col=c('blue','pink'), border = "light gray", ...) {
     n <- length(col)
     plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1), axes = FALSE,
             xlab = "", ylab = "", ...)
     rect(0, 0:(n - 1)/n, .6, 1:n/n,  col = col, border = border)
     ret <- col2rgb(col)
     dimnames(ret)[[2]] <- col
     ret <- t(ret)
     txt <- paste( as.character(col), "(",
        apply( ret, 1, paste, collapse=" "), ")")
     text( rep(.6, n), (0:(n-1)+.5)/n, txt, adj = 0)
     ret <- col2rgb(col)
     dimnames(ret)[[2]] <- col
     t(ret)
}


pals <- function(pp=30){
    n <- length(cc <- colors())
    ii <- 1
    while( ii < n ){

        pal(cc[ii:min(ii+pp,n)], ask = TRUE)
        ii <- ii + pp + 1
    }
}


if ( FALSE ){
    td( col = 1:36, pch = 0:35)
    sampler()
    td(new=T)
    pals()
    pal(colors()[1:30])
}


## brace() now moved to brace.R




change <- function(x,ll) {
    #
    # Modifies elements in a list
    # Ideal for changing ... arguments in calls to panel.groups etc.
    #
    nams <- names(ll)
    for ( ii in  seq_along(ll) ) {
        x[[nams[ii]]] <- ll[[ii]]
    }
    x
}


panel.subgroups <- function( x, y, subscripts,
     subgroups, subgroups.type = c('p','l'),type,
     panel.subgroups = panel.xyplot, ...) {
help = "Use help: ?panel.subgroups"
         subgroups <- as.factor(subgroups)
         levs <- levels(subgroups)
         subgroups.type <- rep( subgroups.type, length.out = length(levs))

         subgroups = subgroups[subscripts]
         for ( i in seq_along( levs) ) {
             sel <- subgroups == levs[i]
             if ( any( sel )) {
                panel.subgroups( x[sel], y[sel], type = subgroups.type[i],
                     subscripts = subscripts[sel], subgroup.number = i,
                     subgroup = levs[i], ...)
             }
         }
}





"</pre>
== Contingency tables ==
<pre>"




### atotal

atotal <- function(arr, FUN = sum, label = 'Total', ...) {
help <- "
atotal                coursefun.R     for PSYC 6140/MATH 6630 05/06

Adds border of sums to an array

Description:

     'atotal' adds by default a border of sums to an array.
     The function FUN may be used instead of 'sum'. Additional
     arguments to FUN can also be given.

Usage:

     atotal( arr , FUN = sum, label = 'Total', ...)

Arguments:

     arr: array, matrix or vector

     FUN: function to be applied to cross sections of arr

     ...: additional arguments to FUN

Details:

Value:

     An array with dimension dim(arr) + 1

References:

Contributed by:  G. Monette  2005-10-10

Modifications:
     2007-12-17: Fixed bug so dimnames is preserved for one-dimensional tables

"
	d <- dim(arr)
	cls <- class(arr)
	dim1 <- FALSE  # to handle 1-dimensional tables
	if (length(d) == 1) {
    dim1 <- TRUE
    dn <- dimnames(arr)
		arr <- c(arr)
		d <- dim(arr)
	}
	if(is.character(FUN))
		FUN <- get(FUN, mode = "function")
	else if(mode(FUN) != "function") {
		farg <- substitute(FUN)
		if(mode(farg) == "name")
			FUN <- get(farg, mode = "function")
		else stop(paste("\"", farg, "\" is not a function", sep = ""))
	}

	if (is.null(d)) {
        ret <- structure(c(arr,FUN(arr,...)), names = c(names(arr), label), class = cls)
        if (dim1) {
            dn[[1]] <- c(dn[[1]],label)
            ret <- structure( ret, dim = length(ret), dimnames = dn)
        }
        return (ret)
    }
	n <- length(d)
	ret <- arr
	ind <- 1:n
	for ( i in n:1) {
		new <- apply(ret,ind[-i],FUN,...)
		ret <- abind( ret, new, i, label)
	}
	class( ret ) <- cls
    ret
}


# TEST:
#  atotal(tab( Status = c(NA,NA,1,1,2)))




###
###  abind
###


abind <- function(arr1,arr2,d,facename="") {

help <- "
abind                coursefun.R     for PSYC 6140/MATH 6630 05/06

Binds two conformable arrays along a dimension

Description:

     'abind' binds two conformable arrays, arr1 and arr2 along
     the dimension d.

Usage:

     abind( arr1, arr2 , d ,  facename )

Arguments:

     arr1, arr2: arrays, matrices or vectors

       d: dimension along which arr1 and arr2 are joined

     ...: Name of extended index if arr2 has one fewer dimensions than arr1

Details:

      dim( arr1 ) and dim( arr2 ) must be equal except in the dth dimension.
      If the length of dim( arr2 ) is 1 less than that of dim( arr1 ), then
      'arr2' is treated as if it had dropped the dth dimension with size 1.

Value:

     The returned value, ret, is an array with dimension dim(arr1)
     except for the dth dimension where dim(ret)[d] == dim(arr1)[d] + dim(arr2)[d].
     If length(dim(arr2)) == length(dim(arr1)) - 1, then arr2 is treated as if
     it dropped a dimension of size 1 in the dth dimension. 'facename' is used
     as the name of the dimnames list for this dimension.

Notes:

     abind is used by atotal

References:

Contributed by:  G. Monette  2005-10-10

Modifications:

"

	d1 <- dim(arr1)
	n1 <- length(d1)
	dn1 <- dimnames(arr1)
	d2 <- dim(arr2)
	n2 <- length(d2)
	dn2 <- dimnames(arr2)

	arenull <- is.null(dn1) & is.null(dn2)
	if (is.null(dn1)){
		dn1 <- lapply( as.list(d1), function(x) seq(1,x))
		dimnames(arr1) <- dn1
	}

	if ( n1 != n2 ) {
		d2 <- d1
		d2[d] <- 1
		dn2 <- dn1
		dn2[[d]] <- facename
		dimnames(arr2) <- NULL
		dim(arr2) <- d2
		dimnames(arr2) <- dn2
		n2 <- n1
	}
	if (is.null(dn2)){
		dn2 <- lapply( as.list(d2), function(x) seq(1,x))
		dimnames(arr2) <- dn2
	}

	# check input for consistency
	# ... later
	#

	perm <- 1:n1
	perm[c(d,n1)] <- c(n1,d)	# perm is an involution

	#
	# permute arr1
	#

	arr.p1 <- aperm(arr1,perm)

	#
	# permute arr2 if length of dimension same as arr1
	#

	arr.p2 <- aperm(arr2,perm)
	dret <- d1[perm]
	dret[n1] <- dret[n1] + d2[d]
	dnret <- dn1
	dnret[[d]] <- c(dnret[[d]],dn2[[d]])

	ret <- c(arr.p1, arr.p2)
	dim(ret) <-  dret

	#
	# permute response back
	#

	ret <- aperm(ret, perm)

	dimnames(ret) <- dnret
	ret
}





###
###   tab
###


otab <- function(...) {
help <- "
abind                fun.R     for PSYC 6140/MATH 6630 05/06

Cross Tabulation and Table Creation Including Missing Values

Description:

     'tab' does the same thing as 'table' except that it includes
     missing values for factors.  The argument 'exclude = NULL' to 'table'
     results in the inclusion of missing values for numeric variable but
     excludes missing values for factors. 'tab' is intended to remedy this
     deficiency of 'table'.

Usage:

     tab(...)

Arguments:
     ...: objects which can be interpreted as factors (including
          character strings), or a list (or data frame) whose
          components can be so interpreted.

Details:

Value:

Notes:

References:

Bugs:
      Does not use argument name as a dimension name, in contrast with 'table'.

Contributed by:  G. Monette  2005-10-10

Modifications:

"
  args <- list(...)
  if( is.list(args[[1]])) args <- args[[1]]
  # for ( ii in 1:length(a)) if ( is.factor( a[[ii]])) a[[ii]] <- factor(a[[ii]],exclude = NULL)
  for ( ii in 1:length(args))  args[[ii]] <- factor(args[[ii]],exclude = NULL)
  do.call("table", args)
}

"</pre>

== Ellipses: data and confidence ==
<pre>"

ellplus <- function ( center = rep(0,2), shape = diag(2), radius = 1, n = 100,
               angles = (0:n)*2*pi/n,
               fac = chol ,
               ellipse = all,
               diameters = all,
               box = all,
               all = FALSE) {
        help <- "
        ellplus can produce, in addition to the points of an ellipse, the
        conjugate axes corresponding to a chol or other decomposition
        and the surrounding parallelogram.
        "
        rbindna <- function(x,...) {
            if ( nargs() == 0) return(NULL)
            if ( nargs() == 1) return(x)
            rbind( x, NA, rbindna(...))
        }
        if( missing(ellipse) && missing(diameters) && missing(box)) all <- TRUE
        circle <- function( angle) cbind( cos(angle), sin(angle))
        Tr <- fac(shape)
        ret <- list (
            t( c(center) + t( radius * circle( angles) %*% Tr)),
            t( c(center) + t( radius * circle( c(0,pi)) %*% Tr)),
            t( c(center) + t( radius * circle( c(pi/2,3*pi/2)) %*% Tr)),
            t( c(center) + t( radius * rbind( c(1,1), c(-1,1),c(-1,-1), c(1,-1),c(1,1)) %*% Tr)))
        do.call( 'rbindna', ret[c(ellipse, diameters, diameters, box)])
    }


dellplus <- function( x, y,  ...) {
    if ( (is.matrix(x) && (ncol(x) > 1))|| is.data.frame(x)) mat <- as.matrix(x[,1:2])
    else if (is.list(x)) mat <- cbind(x$x, x$y)
    else mat <- cbind( x,y)
    ellplus( apply(mat,2,mean), var(mat), ...)

}

ell <- function(center = rep(0,2) , shape = diag(2) , radius = 1, n = 100,
        angles = (0:n)*2*pi/n) {
       circle <- radius * cbind( cos(angles), sin(angles))
       t( c(center) + t( circle %*% fac(shape)))
}


old.cell <-
function (model, which.coef, levels = 0.95, Scheffe = FALSE, dfn = 2,
    center.pch = 19, center.cex = 1.5, segments = 51, xlab, ylab,
    las = par("las"), col = palette()[2], lwd = 2, lty = 1,
    add = FALSE, ...)
{
help <- "
See help for car::confidence.ellipse.lm
except that 'cell' returns the points to form the ellipse
which must be plotted with plot(...,type='l') or lines(...)
-- Use dfn to determine Sheffe dimension, i.e. dfn = 1 to generate ordinary CIs, dfn = 2 for 2-dim CE, etc.
"
    require(car)
    which.coef <- if (length(coefficients(model)) == 2)
        c(1, 2)
    else {
        if (missing(which.coef)) {
            if (has.intercept(model))
                c(2, 3)
            else c(1, 2)
        }
        else which.coef
    }
    coef <- coefficients(model)[which.coef]
    xlab <- if (missing(xlab))
        paste(names(coef)[1], "coefficient")
    ylab <- if (missing(ylab))
        paste(names(coef)[2], "coefficient")
    if(missing(dfn)) {
        if (Scheffe) dfn <- sum(df.terms(model))
        else 2
    }
    dfd <- df.residual(model)
    shape <- vcov(model)[which.coef, which.coef]
    ret <- numeric(0)
    for (level in rev(sort(levels))) {
        radius <- sqrt(dfn * qf(level, dfn, dfd))
        ret <- rbind(ret, c(NA,NA), ell( coef, shape, radius) )
    }
    colnames(ret) <- c(xlab, ylab)
    ret
}

# from Plot3d.R



    dell <- function( x, y, radius = 1, ...) {
        if ( (is.matrix(x) && (ncol(x) > 1))|| is.data.frame(x)) mat <- as.matrix(x[,1:2])
        else if (is.list(x)) mat <- cbind(x$x, x$y)
        else mat <- cbind( x,y)
        ell( apply(mat,2,mean), var(mat), radius = radius, ...)

    }



    ell <- function( center, shape , radius  = 1, n =100) {
          fac <- function( x )  {
              # fac(M) is a 'right factor' of a positive semidefinite M
              # i.e. M = t( fac(M) ) %*% fac(M)
              # similar to chol(M) but does not require M to be PD.
              xx <- svd(x,nu=0)
          t(xx$v) * sqrt(pmax( xx$d,0))
          }
         angles = (0:n) * 2 * pi/n
         if ( length(radius) > 1) {
            ret <- lapply( radius, function(r) rbind(r*cbind( cos(angles), sin(angles)),NA))
            circle <- do.call( rbind, ret)
         }
         else circle = radius * cbind( cos(angles), sin(angles))
         ret <- t( c(center) + t( circle %*% fac(shape)))
         attr(ret,"parms") <- list ( center = rbind( center), shape = shape, radius = radius)
         class(ret) <- "ell"
         ret
    }


    center <- function( obj, ... ) UseMethod("center")

    center.ell <- function( obj, ...) attr(obj, 'parms') $ center



cell <- function( ... )  {
        UseMethod("cell")
        help <- "
        See help for car::confidence.ellipse.lm
        except that 'cell' returns the points to form the ellipse
        which must be plotted with plot(...,type='l') or lines(...)
        -- Use dfn to determine Sheffe dimension, i.e. dfn = 1 to generate ordinary CIs, dfn = 2 for 2-dim CE, etc.
        -- TODO: extend to 3 dimensions if which.coef has length 3
        "
}


cell.wald <-
 function (obj, which.coef = 1:2, levels = 0.95, Scheffe = FALSE, dfn = 2,
    center.pch = 19, center.cex = 1.5, segments = 51, xlab, ylab,
    las = par("las"), col = palette()[2], lwd = 2, lty = 1,
    add = FALSE, ...)
{

# BUGS: works only on first element of glh list
# glh should be restructured to have two classes: waldList and wald

    obj <- obj[[1]]
    coef <- obj$coef[which.coef]
    xlab <- if (missing(xlab))
        paste(names(coef)[1], "coefficient")
    ylab <- if (missing(ylab))
        paste(names(coef)[2], "coefficient")

    dfd <- obj$anova$denDF
    shape <- obj$vcov[which.coef, which.coef]
    ret <- ell( coef, shape , sqrt( dfn * qf( levels, dfn, dfd)))
    colnames(ret) <- c(xlab, ylab)
    ret
}


cell.default <-
function (model, which.coef, levels = 0.95, Scheffe = FALSE, dfn = 2,
    center.pch = 19, center.cex = 1.5, segments = 51, xlab, ylab,
    las = par("las"), col = palette()[2], lwd = 2, lty = 1,
    add = FALSE, ...)
{

    #require(car)
    which.coef <- if (length(coefficients(model)) == 2)
        c(1, 2)
    else {
        if (missing(which.coef)) {
            if (any(names(coefficients(model)) == "(Intercept)"))
                c(2, 3)
            else c(1, 2)
        }
        else which.coef
    }
    coef <- coefficients(model)[which.coef]
    xlab <- if (missing(xlab))
        paste(names(coef)[1], "coefficient")
    ylab <- if (missing(ylab))
        paste(names(coef)[2], "coefficient")
    if(missing(dfn)) {
        if (Scheffe) dfn <- sum(df.terms(model))
        else 2
    }
    dfd <- df.residual(model)
    shape <- vcov(model)[which.coef, which.coef]
    ret <- numeric(0)

    ret <- ell( coef, shape,sqrt(dfn * qf(levels, dfn, dfd)))
    colnames(ret) <- c(xlab, ylab)
    ret
}













"</pre>

== Diagnostics ==
<pre>"


diags <- function(x, ...) UseMethod("diags")


diags.lm <- function(x, y, ..., ask, labels = names(residuals(x)), showlabs = text)
    {
    # diags.lm
    # graphical diagnostics for lm, locally first-order for glm
    # enlarged version of plot.lm with emphasis on diagnostics
    # G. Monette, Dec. 94
    # modified Nov. 97, May 98
    # Slight modification to pairs adding labels, Jan 03
	if(!missing(ask)) {
		op <- par(ask = ask)
		on.exit(par(op))
	}
	form <- formula(x)
	f <- predict(x)
	r <- residuals(x)
	nams <- names(r)
	if(!missing(labels)) {
		nams <- names(residuals(x))	#
    # if labels not same length as residuals assume it's a vector
    # of len == original data and select elements included in residuals
		if(length(nams) != length(labels))
			labels <- labels[nams]
	}
	ret <- NULL
	if(missing(y)) {
		y <- f + r
		yname <- deparse(form[[2]])
	}
	else yname <- deparse(substitute(y))
	fname <- paste("Fitted:", deparse(form[[3]]), collapse = " ")
	plot(f, y, xlab = fname, ylab = yname, main = "Dependent var. vs. Predicted",
		...)
	abline(0, 1, lty = 1)
	lines(supsmu(f,y))
	showlabs(f, y, labels,...)
    #
    # get influence diags and model matrix while looking at first plot
    #
	lmi <- lm.influence(x)
	hat <- lmi$hat
	sigma <- lmi$sigma	# drop 1 sigma
	mm <- scale(model.matrix(x), scale = F)	# centres each column
	mp <- predict(x, type = "terms")
	comp.res <- mp + r	# effect + residual
    #
    # Absolute residual vs. predicted
    #
	plot(f, abs(r), xlab = fname, ylab = deparse(substitute(abs(resid(x)))),
		main = "Absolute Residual vs. Predicted", ...)
	showlabs(f, abs(r), labels, ...)	#
    #
    # Normal quantile plot
    #
	zq <- qqnorm(r, main = "Normal Quantile Plot", ylab = "Residual", sub
		 = fname)
	qqline(r)
	showlabs(zq, labels,...)	#
    #
    # Symmetry plot of residuals (Lawrence C. Hamilton, Regression with
    #       Graphics, Duxbury, 1992)
	n <- length(r)
	r.o <- sort(r)
	half <- (n + 1)/2
	if(n %% 2 == 1) {
    # n is odd
		med <- r.o[half]
		below <- med - r.o[half:1]
		above <- r.o[half:n] - med
	}
	else {
    # n is even
		med <- sum(r.o[c(half, half + 1)])/2
		below <- med - r.o[(n/2):1]
		above <- r.o[(n/2 + 1):n] - med
	}
	opt <- par(pty = "s")
	ran <- range(c(below, above))
	plot(below, above, main = "Symmetry plot of residuals", xlab =
		"Distance below median", ylab = "Distance above median", xlim
		 = ran, ylim = ran)
	abline(0, 1, lty = 2)
	par(opt)	#
	
    #
    # Studentized residual vs. leverage
    #
    
	std.r <- r/(sigma * sqrt(1 - hat))
	plot(hat, std.r, xlab = "Leverage (hat)", ylab = yname, sub = fname,
		main = "Studentized residual vs. Leverage", ...)
	showlabs(hat, std.r, labels,...)	#	plot(lmi$sig, std.r)	#

    #
    # effect of dropping one observation DFBETA
    #
    
	nams <- dimnames(lmi$coefficients)[[1]]
	pairs(lmi$coefficients)
	pairs(lmi$coefficients, panel = function(x,y,nams){
		points(x,y)
		text(x,y,nams)
	}, nams = nams)

	# main = "Effect of dropping one case", sub = fname)
	invisible(0)
}

model.matrix.lme <- function( fit , data = fit$data, ...)
                 naresid( fit$na.action, model.matrix( formula(fit)[-2], data = data ,...))

diags.lme <- function( ... ) cat("Being implemented")



"
</pre>
== vplot -- a plot function for matrix algebra ==
<pre>
"

##
##  vplot  plots columns of a 2 x n matrix
##  Transferred to coursfun: Nov. 15, 2005

vplot <- function( mat , type = 'p', new = F,  pch = 16, pop = 0, ...) {
help <- "
vplot    - plots the columns of a 2 x n matrix or a vector of length 2
         - vplot adds to the current plot resizing it to include all plotted
           objects in a 'euclidean' frame
         - to start a new plot, use 'new = TRUE'
         - to remove the last element added use 'vplot(pop=1)'
         Associated functions:
         - vell( mean, var) generates an ellipse, default = unit circle
         - vbox() generates a box
         - vobj() generates a circle in a box
         - orthog(theta) generates an orthog matrix rotating through angle theta
         - orthog.proj generates the matrix of an orthog. projection into span (x)
         - vmat( .... ) generates a 2 by n matrix
         Examples:
           vplot( new = TRUE )
           vplot( vell(), 'l' )
           vplot( cbind(c(3,1),c(1,4)) %*% vell())
           vplot( pop = 1)
           vplot( cbind(c(3,1),c(1,4)) %*% vell(), type = 'l', col = 'red')
"
     if (  new || !exists(".vplot")) assign(".vplot", list(list(x=0,y=0,type='n')),pos=1)
     a <- .vplot
     if ( ! missing(mat) ) {
        mat <- cbind(mat)
        if ( type == 'v' ) {
           zz <- rbind( 0*mat, mat, mat/0)
           mat <- matrix( c(zz), nrow = 2)
           type = 'b'
        }
        d <- dim(mat)
        if ( d[1] != 2 && d[2] == 2){
           mat <- t(mat)
           warning("mat is n x 2 and has been transposed")
        }
        a <- c(a,list( list(x=mat[1,],y = mat[2,],type=type, pch = pch, ...)))
     }
     dat <- NULL
     for ( i in seq( along = a )) {
         dat <- c( dat, a[[i]]$x, a[[i]]$y)
     }
     # print(a)
     par ( pty = 's')
     plot( range(na.omit(dat)), range(na.omit(dat)), type = 'n', xlab = '', ylab ='')
     if ( pop > 0 ) {
        keep <- 1:max(1,(length(a)-(pop+1)))
        a <- a[keep]
     }
     abline( h = 0, v = 0)
     for ( i in seq( along = a)) do.call('points', a[[i]])
     assign(".vplot", a, pos = 1)
     invisible(a)
}

vell <- function(...) t( ell(...))
vbox <- function(...) cbind( c(-1,-1), c(-1,1), c(1,1), c(1,-1), c(-1,-1))
vobj <- function(...) {
     cbind( vell(), NA, vbox(), NA, c(0,-1),c(0,1), NA, c(-1,0), c(1,0))
}
vsquare <- function(...) vmat( 0,0,0,1,1,1,1,0,0,0)

vmat <- function(...) {
help <- "
vmat creates a matrix entering data column by column
"
     aa <- list(...)
     aa <- do.call('c', aa)
     matrix(aa, nrow = 2)
}

orthog <- function( theta ) cbind( c( cos(theta), sin(theta)), c( - sin(theta), cos(theta)))
orthog.proj <- function ( x ) {
       x <- cbind(x)
       x %*% solve(t(x) %*% x , x)
}


"
</pre>
== Read.spss and Read.dta ==
<pre>
"

###
###  trim
###

  trim <- function(x) {
       help <- "
trim in fun.R
  removes trailing blanks from character variables or from
  factor levels
  Use mainly to trim .dta files produced with SPSS
" 
       UseMethod("trim")
}

  trim.data.frame <- function(x) {
      for ( nn  in names(x)) x[[nn]] <- trim(x[[nn]])
      x
  }
  trim.factor <- function( x ) {
      levels(x) <- trim(levels(x))
      x
  }
  trim.character <- function( x ) {
      x[] <- sub(" +$", "", x )
      x[] <- sub("^ +", "", x )
      x
  }
  trim.default <- function(x) x
  
  Read.spss <- function( ... ) {
            require("Hmisc")
            dd <- spss.get ( ... )
            trim( dd )
  }
  
  Read.dta <- function ( ... ) {
help <- "
  Read.dta reads Stata files using 'read.dta' in 'library(foreign)'
  This appears to be an ideal way of importing spss files in order
  to keep full variable names. Direct use of 'read.spss' on a SPSS
  '.sav' file abbreviates variable names to 8 characters.
  Note: missing.type = TRUE produces warnings.
"
           require("foreign")
           #  dd <- read.dta(... , missing.type = TRUE)  # Note: missing.type = T produces warnings.
           dd <- read.dta(...)
           cls <- sapply(dd,class)
           ch.nams <- names(dd) [ cls == "character" ]
           for ( nn in ch.nams ) dd[[nn]] <- factor(trim(dd[[nn]]) )
           dd
  }
  
  
  Write.dta <- function( ... ) {
       require("foreign")
       write.dta( ..., version = 7)
  }
  
  Write.spss <- function( dataframe, file, ... ) {
        require(foreign)
        dname <- deparse(substitute(dataframe))
        disp(dname)
        cls <- sapply( dataframe, class )
        for ( nn in names(dataframe)[cls == "Date"] ){
            dataframe[[nn]] <- as.character( dataframe[[nn]], "%Y/%m/%d")
        }
        if ( any ( cls == "Date")) {
            cat("\nOpen .dta file in SPSS and convert following variables to dates\nwith yyyy/mm/dd format:\n")
            for ( nn in names(dataframe) [ cls == "Date" ] ) cat("       ", nn, "\n")
        }
        for ( nn in names(dataframe)[cls == "factor"]) {
            dataframe[[nn]] <- as.character( dataframe[[nn]])
        }
        if ( missing(file) )  file <- paste(dname,".dta", sep ="")
        else file <- sub("\\.dta$|\\.DTA$|$",".dta",file)

        cat(paste("\nData saved in", file,"\n"))
        write.dta( dataframe, file, version = 7, ...)
  }
  
   # zd <- data.frame( x=1:10, a = LETTERS[1:10], d=seq(as.Date("2000-01-01"), by ="4 months", length.out = 10))
   # zd
   # Write.spss(zd)





"
</pre>
== RDC functions -- needs to be organized ==
<pre>
"

capply <- function ( x ,... ) UseMethod("capply")
   
capply.default <- function ( x, by, FUN , ...) {
            if (inherits(by,'formula')) by <- model.frame( by , x , na.action = na.include)
            ret <- unsplit ( lapply ( split ( x , by ), FUN, ...), by )
            if ( !is.null( dim(ret)) && length(dim(ret)) ==1) ret <- c(ret)
            ret
}



###
###   pchisq.mix
###

   pchisq.mix <- function( q, df , mix = rep(1,length(df))/length(df) ) {
         # returns cdf for mixture of chi-squares. Usefule for testing
         # random effects in mixed models
           pc <- function( df ) if ( df == 0 ) 1*(q >= 0) else pchisq( q, df )
           sum ( sapply( df, pc) * mix)
   }


   ch <- as.character

   select.first <- function( x , ... ) {
       ret <- unique( x [ !is.na(x) ])
       if ( length(ret) > 1 && is.character( ret ) && ( "" %in% ret)) ret <- ret [ ret != ""]
       if ( length(ret) > 1 ) {
        cat("\n\n============= Multiple values in select.first ==========\n")
        for ( i in 1: length(ret)) cat("\nValue",i,": <<",ret[i],">>\n")
       }
       ret [1]
   }

   fill <- function( x, ... ) UseMethod("fill")

   fill.factor <- function( x, by, FUN = select.first, ...) {
      levs <- levels(x)
      ret <- capply( ch( x), by, FUN , ... )
      factor( ret, levels = levs)
   }

   fill.Date <- function( x, by, FUN = select.first, ...) {
      as.Date( capply(ch(x), by, FUN, ...))
   }

   fill.default <- function( x , by , FUN = select.first, ...) {
      capply( x, by, FUN, ...)
   }

   fill.data.frame <- function ( x, by ,... ) {
          ret <- list()
          for ( nn in names(x)) {
              ret[[nn]] <- fill( x[[nn]], by)
          }
          as.data.frame(ret)
   }


    ##
    ##   cvar: V0.1 August 15, 2006
    ##   Creating contextual variables for categorical variables
    ##   
    ##   cvar is designed to create contextual variables 
    ##   for factors, as well as for numerical variables.
    ##   If a factor has g levels, convar will create a
    ##   matrix with g-1 columns each of which is the within group
    ##   mean value of the correponding column of the "contrast"
    ##   matrix for the factor.
    ##
    

    cvar <- function( x, id ,... ) {
        help = "
        cvar: creates contextual group mean variables within levels of 'id'.\n
              If 'x' is a factor, 'cvar' returns a matrix labelled so that it\n
              is consistent with labels generated for coding variables for 'x'.\n
              Example:\n
               \n
                dd <- data.frame(x= 1:100, id = rep( LETTERS[1:10], each = 10))\n
                dd$a <- factor(sample( c('a','b','c'), 100, replace = T))\n
                dd$y <- dd$x + rep(rnorm(10), each = 10) + rnorm(100) + as.numeric(dd$a)\n
                library(nlme)\n
                fit <- lme( y ~ x + cvar(x,id), dd, random = ~ 1 + dvar(x,id) | id)\n
                anova( fit , type = 'm')\n
                                        \n
              The output of 'anova' can be used to test whether a contextual effect\n
              needs to be included in the model.\n
                                                \n
              See also: dvar for group-mean centering: x - cvar(x, id)\n
        "
        UseMethod("cvar")
    }
    
    cvar.factor <- function( x, id, ... ) {
        mat <- contrasts( x) [ x,]
        ret <- cvar(mat, id, ...)
        colnames(ret) <- colnames(mat)
        ret
    }
    
    cvar.default <- function( x, id, ... ) {
        if ( is.matrix (x) ) {
            if ( dim(x)[2] == 1) return( cvar( x[,], id, ...))
            else {
                ret <-  cbind( cvar(x[,1], id, ...), cvar(x[,-1],id,...))
                colnames(ret) <- colnames(x)
                return( ret )
            }
        } else {
            capply( x, id, mean, na.rm = T)
        }
    }
    

    dvar <- function( x, id ,... ) {
        help = "
        dvar: produces group mean centering: x - cvar(x, id)
        See 'cvar'
        "
        UseMethod("dvar")
    }
    
    dvar.factor <- function( x, id, ... ) {
        mat <- contrasts( x) [ x,]
        ret <- mat - cvar(mat, id, ...)
        colnames(ret) <- colnames(mat)
        ret
    }
    
    dvar.default <- function( x, id, ... ) {
        if ( is.matrix (x) ) {
            if ( dim(x)[2] == 1) return( dvar( x[,], id, ...))
            else {
                ret <-  cbind( dvar(x[,1], id, ...), dvar(x[,-1],id,...))
                colnames(ret) <- colnames(x)
                return( ret )
            }
        } else {
            x - capply( x, id, mean, na.rm = T)
        }
    }

##
##  sum
##
    
cvars <- function(  x, by, ...) {
      if ( length(x) == 1 && x == 1) {
            n <- nrow(as.data.frame(by))
            capply( rep(1,n), by, sum)
      } else {
            capply( x, by, sum, ...)
      }
}



na20 <- function(x) {
     x[is.na(x)] <- 0
     x
}

describe.vector <-
function (x, descript, exclude.missing = TRUE, digits = 4, weights = NULL,
    normwt = FALSE, ...)
{
    # GM: modified by rounding n and missing
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
    else n <- round(sum(weights[present]),2)
    if (exclude.missing && n == 0)
        return(structure(NULL, class = "describe"))
    missing <- round(sum(weights[!present], na.rm = TRUE),2)
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
        dimnames(counts)[[1]] <- c("Frequency", "%")
    }
    z$values <- counts
    structure(z, class = "describe")
}


na.include  <- function (obj)
{
     # from library(Hmisc)
    if (inherits(obj, "data.frame"))
        for (i in seq(along = obj)) obj[[i]] <- na.include(obj[[i]])
    else {
        if (length(levels(obj)) && any(is.na(obj)))
            obj <- factor(obj, exclude = NULL)
    }
    obj
}


summ <- function(x,...) UseMethod("summ")

summ.lmer <- function(x, ...) {
               ret <- c(AIC = AIC(x@logLik), BIC= BIC(x@logLik), logLik=x@logLik)
               ret
}

pr <- function(x,...) {
   # print to cut and paste as input
   UseMethod("pr")
}
pr.default <- function(x,pre="\t",...) {
          # cat('\nc(')
          for ( xx in x) cat(pre,'"',xx,'",\n',sep='')
          invisible(x)
}

##  reorder.factor has been removed because the new method in
##  in the base package does the same thing
##reorder.factor <- function (x, v, FUN = mean, ...) {
##        warning("gm version produces ordinary factor -- not ordered factor")
##               factor(x, levels = levels(x)[order(tapply(v, x, FUN, ...))])
##}
## NOTE: reorder.factor in <base> is hidden. You must use is as reorder(x)

reorder.default <- function(x,...) reorder( factor(x),...)  # so reorder works on non-factors


Q <- function(x, verbose = 0) {
    # find the Q matrix of a qr decomposition with possible NAs
    miss <- apply(x, 1, function(xx) any(is.na(xx)))
    xc <- x[!miss,]
    qf <- qr.Q(qqr<-qr(xc))
    
    if(verbose > 0) {
               cat("xc:", dim(xc),'\n')
               cat("qf:", dim(qf), '\n')

               print(qf)
    }
    if( ncol(xc) > ncol(qf)) xc <- xc[,1:ncol(qf)]
    ip <- sign(apply(qf*xc,2,sum))   # sign of reln between Q and X
    qf <- qf * rep(ip, rep( nrow(qf), length(ip)))   # change sign of each row
    ret <- array(NA, dim = c(nrow(x),ncol(qf)))
    rownames(ret) <- rownames(x)
    colnames(ret) <- colnames(xc)
    ret[!miss,] <- qf
    attr(ret,'rank') <- qqr$rank
    attr(ret,'miss') <- miss
    ret
}

Contrasts <- function(x) UseMethod("Contrasts")

Contrasts.default <- function(x) contrasts(x)

Contrasts.lmer <- function(x) {
       dd <- x@frame[1,]
       ret <- list()
       for ( nn in names(dd)) {
           vv <- dd[[nn]]
           if (is.factor(vv)) ret[[nn]] <- contrasts(vv)
       }
       vv
}


comp.old <- function(fit, form, varname, varpattern = vname, data = fit@frame, ...) {
     ## Computing regression components
     ## currently works form lmer only
     ## idea is to return a data frame with value a variable and corresponding fitted values
     ## for the 'component' of that variable or variables
     ##
     ## This can be used to fit a model to a prediction data.frame with only
     ## the necessary variables defined!!
     ## However BEWARE not to use 'unsafe' transformations (e.g. Q, poly)
     ##
     model.mat <- model.matrix(form, data,...)
     #print(dim(model.mat))
     ret <- data[rownames(model.mat),varname,drop = F]
     fe <- fixef(fit)
     pos.fe <- grep( varpattern, names(fe))
     pos.mat <- grep( varpattern, colnames(model.mat))
     ret$comp <- c( model.mat[,pos.mat] %*% fe[pos.fe] )
     attr(ret,"predictors") <-  names(fe)[pos.fe]
     ret
}



comp.bak <- function(fit, form, varname, varpattern = vname, data = fit@frame, ...) {
     ## Computing regression components
     ## currently works form lmer only
     ## idea is to return a data frame with value a variable and corresponding fitted values
     ## for the 'component' of that variable or variables
     ##
     ## This can be used to fit a model to a prediction data.frame with only
     ## the necessary variables defined!!
     ## However BEWARE not to use 'unsafe' transformations (e.g. Q, poly)
     ##
     ## Arguments:
     ##
     ## fit : model for which component to be estimated
     ##
     ## form : portion of model formula to generate model variables needed for fitting
     ##
     ## varname: variable names to be included in output data frame
     ##
     ## varpattern:  regular expression so select effects in component
     ##
     ## data:  the 'parent' data.frame -- often a prediction data frame.
     ##        When using the original data frame, it is often necessary to
     ##        'na.action = na.omit' for '...'
     ##
     model.mat <- model.matrix(form, data,...)
     #print(dim(model.mat))
     ret <- data[rownames(model.mat),varname,drop = F]
     fe <- fixef(fit)
     effect.names <- grep( varpattern, names(fe), value = T)
     ret$comp <- c( model.mat[,effect.names] %*% fe[effect.names] )
     attr(ret,"predictors") <-  effect.names
     ret
}


comp <- function(fit, form = terms(fit), varname, varpattern = "", data = getData(fit), ...) {
# this is the version that worked for RDC code
# the version only returns the adjusted predicted value
# The reason for returning a portion of the original data frame was to ensure
# correspondence between components and data values.
# With 'na.comp' the components are padded to match the original data.
# and can do so in a form that conforms to the original data mat
## for lmer: data = fit@frame
     ## Computing regression components
     ## currently works for lmer only
     ## idea is to return a data frame with value a variable and corresponding fitted values
     ## for the 'component' of that variable or variables
     ##
     ## This can be used to fit a model to a prediction data.frame with only
     ## the necessary variables defined!!
     ## However BEWARE not to use 'unsafe' transformations (e.g. Q, poly)
     ##
     ## Arguments:
     ##
     ## fit : model for which component to be estimated
     ##
     ## form : portion of model formula to generate model variables needed for fitting
     ##
     ## varname: variable names to be included in output data frame
     ##
     ## varpattern:  regular expression to select effects in component
     ##
     ## data:  the 'parent' data.frame -- often a prediction data frame.
     ##        When using the original data frame, it is often necessary to
     ##        'na.action = na.omit' for '...'
     ##
     
     getData <- function(x) UseMethod("getData")
     getData.lme <- function(x) x$data
     getData.lmer <- function(x) x@frame
     disp(form)
     disp( dim(data))
     model.mat <- model.matrix(form, data,...)
     disp(dim(model.mat))
     ret <- data[rownames(model.mat),,drop = F]
     fe <- fixef(fit)
     effect.names <- grep( varpattern, names(fe), value = TRUE)
     ret$comp <- c( model.mat[,effect.names] %*% fe[effect.names] )
     attr(ret,"predictors") <-  effect.names
     ret
}



    com <- function(fit, varpattern = "", form = terms(fit), data = getData(fit), ...) {
    # this is a new version of 'comp' that
    # only returns the adjusted predicted value
    # The reason for returning a portion of the original data frame was to ensure
    # correspondence between components and data values.
    # With 'na.com' the components are padded to match the original data.
    # and can do so in a form that conforms to the original data mat
    ## for lmer: data = fit@frame
         ## Computing regression components
         ## currently works for lmer only
         ## idea is to return a data frame with value a variable and corresponding fitted values
         ## for the 'component' of that variable or variables
         ##
         ## This can be used to fit a model to a prediction data.frame with only
         ## the necessary variables defined!!
         ## However BEWARE not to use 'unsafe' transformations (e.g. Q, poly)
         ##
         ## Arguments:
         ##
         ## fit : model for which component to be estimated
         ##
         ## form : portion of model formula to generate model variables needed for fitting
         ##
         ## varname: variable names to be included in output data frame
         ##
         ## varpattern:  regular expression to select effects in component
         ##
         ## data:  the 'parent' data.frame -- often a prediction data frame.
         ##        When using the original data frame, it is often necessary to
         ##        'na.action = na.omit' for '...'
         ##

         getData <- function(x) UseMethod("getData")
         getData.lme <- function(x) x$data
         getData.lmer <- function(x) x@frame
         #disp(form)
         #disp( dim(data))
         model.mat <- model.matrix(form, data,...)
         #disp(dim(model.mat))
         # ret <- data[rownames(model.mat),,drop = F]
         fe <- fixef(fit)
         effect.names <- grep( varpattern, names(fe), value = TRUE)
         ret <- c( model.mat[,effect.names] %*% fe[effect.names] )
         names(ret) <- rownames(model.mat)
         # return only components corresponding to fit$resid
         retnames <- rownames(fit$resid)   # this main only work for 'lme'
         ret <- ret[retnames]
         attr(ret,"predictors") <-  effect.names
         ret
    }

    na.com <- function( fit, ...) {
        ret <- com(fit,...)
        prs <- attributes(ret)$predictors
        ret <- na.pad( fit, ret)
        attr(ret,"predictors") <- prs
        ret
    }

    na.comres <- function( fit, varp = "", level = 0, ...) na.com( fit, varp = varp, ...) + na.resid(fit, level = level)
    na.getResponse <- function( fit) na.pad( fit, getResponse(fit))

    #length(names(na.comp(fit.moda.age,varp='Month')))
    #dim(dd)


#  glh( fit, Lall(fit, 'Wave'))

cond <- function(x) {
     # reporting on conditioning of matrix
     Svd <- svd(x)
     ret <- list(svd = Svd, dim = dim(x), rank = qr(x)$rank, condition = Svd$d[1]/Svd$d[length(Svd$d)])
     if(nrow(x) == ncol(x)) ret <- c(ret, cor= list(cor(x)))
     ret
}



Vcov.rdc <- function( fit, L  = Lmat(fit,"") ) {
     # variance of etahat = L beta.hat
     vc <- getFix(fit)$vcov
     L %*% vc %*% t(L)
}

Vcor.rdc <- function(fit, L = Lmat(fit,"")) {
     ret <- cov2cor(vc <- Vcov(fit, L))
     sv <- svd(vc)$d

     attr(ret,'condition') <- sv[1]/sv[length(sv)]
     attr(ret,'class') <- "correl"
     ret
}


print.correl <- function(x) {
      correl <- format(round(x, 3), nsmall = 3,
                  digits = 4)
      correl[!lower.tri(correl)] <- ""
      if (!is.null(cond <- attr(correl,"condition"))) {
         attr(correl,"condition") <- NULL
      }
      print(correl, quote = FALSE)
      if(!is.null(cond)) cat("Condition:",cond,"\n")
      invisible(x)
}


glh.rdc <- function(fit, Llist, help = FALSE, clevel = 0.95, debug = FALSE) {
    if(help) {
         cat("help!!!")
         return(0)
    }
    if( !is.list(Llist) ) Llist <- list(Llist)
    ret <- list()
    fix <- getFix(fit)
    beta <- fix$fixed
    vc <- fix$vcov
    dfs <- fix$df
    for ( ii in 1:length(Llist)) {
        ret[[ii]] <- list()

        L <- rbind(zz <-  Llist[[ii]])
        heading <- attr(zz, "heading")
        nam <- names(Llist)[ii]
        
        ## Anova

        qqr <- qr(t(L))
       #if(debug) print(qqr)
        L.rank <- qqr$rank
        L.full <- t(qr.Q(qqr)) [ 1:L.rank,,drop=F]
        if(F) {
                  cat("\n+++++++++++++++++++\n")
                  print(nam)
                  print(L)
                  print(svd(L)$d)
                  print(L.full)
                  print(svd(L.full)$d)
         }
         if (debug) {
            print(L)
            print(dim(L.full))
            print(dim(vc))
         }
        vv <- L.full %*% vc %*% t(L.full)
        eta.hat <- L.full %*% beta
        Fstat <- (t(eta.hat) %*% qr.solve(vv, eta.hat))/L.rank
        Fstat2 <- (t(eta.hat) %*% solve(vv) %*% eta.hat)/L.rank
        included.effects <- apply(L, 2, function(x) sum(abs(x))) != 0
        denDF <- min(dfs[included.effects])
        numDF <- L.rank

        ret.anova <- rbind(c(numDF, denDF, Fstat,Fstat2, pf(Fstat, numDF, denDF, lower.tail = FALSE)))
        colnames(ret.anova) <- c("numDF","denDF","F value","F2","Pr(>F)")
        rownames(ret.anova) <-  nam
        ret[[ii]]$anova <- ret.anova
        
        ## Estimate
        
        etahat <- L %*% beta
        etavar <- L %*% vc %*% t(L)
        etasd <- sqrt(diag(etavar))
        denDF <- apply( L , 1, function(x,dfs) min(dfs[x!=0]), dfs = dfs)
        aod <- cbind(
                 c(etahat),
                 etasd,
                 denDF,
                 c(etahat/etasd),
                 2*pt(-abs(etahat/etasd), denDF))
        colnames(aod) <- c("Estimate","Std.Error",'DF','t value','Pr(>|t|)')
        if( !is.null(clevel) ) {
            hw <- qt( 1-(1-clevel)/2, denDF) * etasd
            aod <- cbind(aod, LL = etahat - hw, UL = etahat + hw)
            labs <- paste( c("Lower","Upper"), format(clevel))
            colnames(aod)[ ncol(aod) + c(-1,0)] <- labs

        }
        #aod <- as.data.frame(aod)
        rownames(aod) <- rownames(L)
        ret[[ii]]$estimate <- aod
        
        ## Vcov
        
        ret[[ii]]$vcov <- Vcov( fit, L)
        ret[[ii]]$vcor <- Vcor(fit,L)
        ret[[ii]]$L <- L
        ret[[ii]]$L.full <- L.full
        if ( !is.null(heading)) attr(ret[[ii]],"heading") <- heading
    }
    names(ret) <- names(Llist)
    attr(ret,"class") <- "glh"
    ret
}





formatCoefmat <- function(x ,digits = 6, pdigits = digits-1 ) {
     pformat <- function(x, digits) {
             x <- format(xx <- round(x,digits))
             x[ as.double(xx) == 0 ] <- paste(c("<.",
                              rep('0',digits-1),"1"), collapse = "")
             x
     }
     xx <- array("",dim=dim(x), dimnames = dimnames(x))
     for ( i in 1:ncol(xx)) {
         xx[,i] <- format(round(x[,i],digits),digits = digits)
     }
     if ( length( isp <- grep("^Pr\\(",colnames(x))) > 0) {
         xx[,isp[1]] <- pformat( x[,isp[1]], digits = pdigits)
     }
     xx
}

print.glh.rdc <- function(x, round = 6, pround = round - 1, L  = TRUE, cov = TRUE, ...) {
     # should round by SD, i.e. keep 3 sig digits for sd and others rounded accordingly


     rnd <- function(x, digits) {
             if ( is.numeric( x)) x <- round(x, digits = digits)
             format(x)
     }
     for ( ii in 1:length(x)) {
         nn <- names(x)[[ii]]
         tt <- x[[ii]]
         ta <- tt$anova
         tap <- array("", dim = dim(ta), dimnames = dimnames(ta))

         # ta[["p-value"]] <- pformat( ta[["p-value"]], digits = pround)
         ## print(as.data.frame(ta,row.names = nn))
         cat("\n",nn,"\n",sep='')

         print(formatCoefmat(ta,digits=round,pdigits=pround),quote=F,right=T)
         cat("\n")
         te <- tt$estimate
         #tret <- te
         #mode(tret) <- 'character'
         #tret[,'p-value'] <- pformat( te[,'p-value'], digits = pround)
         #if( !is.null(round)) {
         #    for (i in 1:ncol(te) ) {
         #        tret[,i] <- rnd(te[,i], digits = round)
         #    }
         #print(tret,quote=F)
         if(!is.null(zhead <- attr(tt,'heading'))) cat(zhead,"\n")
         print(formatCoefmat(te,digits=round,pdigits=pround),quote=F,right=T)
         if (L == TRUE ) {
            cat("\nL:\n")
            print(tt$L)
            if( dim(tt$L.full)[1] < dim(tt$L) [1]) {
             cat("\nL (full rank):\n")
             print(tt$L.full)
            }
         }
         if ( cov == TRUE ) {
            cat("\nVar-Cov of estimates:\n")
            print(tt$vcov)
            cat("\nCorrelations:\n")
            print(tt$vcor)
         }

     }
     invisible(x)
}


##z <- glh(fit, list('language'=Lmat(fit, "language")))
##z

##z <- glh(fit, list('Language | Year=1998'=Ldiff(fit, "language",ref="Qc French Multi")))
##z


####################################################
####################################################  ABOVE is OKAY
xanova <- function(object,...) UseMethod("xanova")

#xanova.lmer <- function(fit, Llist ) {
#     if ( is.list(Llist) )
#}


xanova.lmer <- function( fit, Llist , df = NULL, clevel = .95) {
       # performs a Wald test on an object that has a fixef and a vcov methods
       warning("xanova.lmer uses Chi-Square tests")
       ret <- list()
       for ( ii in 1:length(Llist) ) {
           L <- rbind(Llist[[ii]])

           # anova step

           QR <- qr(L)
           R <- qr.R(QR)
           dfH <- QR$rank
           eta <- R %*% fixef(fit)
           vv <- R %*% vcov(fit) %*% t(R)
           chisq <- t(eta) %*% qr.solve(vv, eta)
           test <- list(ChiSquare = chisq, DF = dfH, "p-value" = 1-pchisq(chisq,dfH))
           ret[[ii]]$anova <- test
           
           # estimation
           
           eta <- L %*% fixef(fit)
           vv <- diag( L %*% vcov(fit) %*% t(L))
           etasd <- sqrt(vv)
           zval <- c(eta/etasd)
           aod <- cbind(Estimate=c(eta), Std.Error = etasd,
               "z-value" = zval, "p-value" = 2*pnorm(-abs(zval)))
           if( !is.null(clevel) ) {
               hw <- qnorm(1-(1-clevel)/2) * etasd
               aod <- cbind( aod, LL = eta - hw, UL = eta + hw)
               labs <- paste(c("Lower","Upper",format(clevel)))
               colnames(aod) [ ncol(aod) + c(-1,0)] <- labs
           }
           aod <- as.data.frame(aod)
           class(aod) <- c('estimate.lme','data.frame')
           
           ret[[ii]]$estimate <- aod
        }
}

####################################   NEED TO TEST ABOVE

            
           
           





enc <- function(x) {
		## this function will use the coding table in Stats Can documentation
		## to create a factor with the right names
		#
		tonum <- function(x) as.numeric(gsub(",","",as.character(x)))
		ff <- function(x) format(x, big.mark=',')
		tran.table <- scan(what=list('character','integer','character','character'),
			flush = TRUE)
		# Brief report
		sam <- tonum(tran.table[[3]])
		pop <- tonum(tran.table[[4]])
		samp.pct <- 100 * (sam / pop) / ( sum(sam,na.rm=T)/sum(pop,na.rm=T)) 
		print( data.frame(
			Code = tran.table[[2]], Content = tran.table[[1]], Sample = ff(sam), Popn = ff(pop),
				Sampling.Pct = round(samp.pct,1))) 
		tran( tran.table[[2]], tran.table[[1]], x , tofactor = TRUE)
	}

apct <- function(x,MARGIN=1) {
		if( length(dim(x)) == 1) {
			ret <- cbind( N = x, pct = 100*x/sum(x,na.rm=T))
			ret <- rbind( ret, Total = apply(ret,2,sum,na.rm=T))
			print( round(ret,1))
			return(invisible(ret))
		}	
		# report a table
		ret <- list( N=atotal(x), pct = 100 * acond(x,MARGIN) )		
		cat("\nN:\n")
		print( ret[[1]])
		cat("\nPercentage:\n")
		print( round( ret[[2]],1))
		invisible(ret)
	}

  arep <- apct   # old names

## From library gm

write.sas <- function( df , datafile="G:/SAS/R2SAS.txt", codefile = "G:/SAS/R2SAS.sas"){
	debug <- F
	pr <- function(x) {
		cat(deparse(substitute(x)),"\n")
		print(x)
		cat("==========\n")		
		cat(x)
		cat("\n=============\n")
	}
	
	lrecl <- 256
	if(!debug) { 
		write.table(df, file = datafile, row = FALSE, col = FALSE, 
        		sep = ";", na = ".")
		# compute lrecl
		lines <- scan(file = datafile, what = list("character"), sep = "\n")
		lrecl <- max(nchar(lines))
	}
    	nms <- names(df)
	nms.sas <- gsub("\\.","_",nms)
	## Check for duplicate SAS names

	if ( length( unique(toupper(nms.sas))) != length(nms.sas)) {
		ind <- duplicated( toupper(nms.sas))
		ind.rev <- duplicated( rev(toupper(nms.sas)))
		cat("Warning:\n")
		cat("The following R names may yield duplicate SAS names",
			"\n", paste(nms[ind | rev(ind.rev)],collapse=" "),"\n")
		warning("Possible duplicate SAS names")
	}
	
	factors <- sapply(df, is.factor) | sapply(df, is.character)
	## check for odd types
	classes <- sapply(df, class)
	odd.classes <- setdiff( sapply(df,class), c('numeric','factor','character') )
	if ( length(odd.classes) > 0 ) {
		cat("Warning:\n")
		cat("The following variables have classes that might not be handled properly by SAS\n")
		print( classes[ grep( odd.classes, classes)])
		cat("\n")
	}

	factor.names <- nms[factors]
	factor.names.sas <- nms.sas[factors]	
	dollarsign <- ifelse( factors, "$","")
	factor.lengths <- sapply( df[factor.names], function(x) {
		if(is.factor(x)) max(nchar(levels(x) ) ) else 
			max(nchar(x))
		})
	length.stmt <- paste(paste( "   ", factor.names.sas, "$", factor.lengths,"\n"),collapse = "")
	length.stmt <- paste( "LENGTH\n", length.stmt, ";\n")
	if (debug) pr(length.stmt)
	input.stmt <- paste(paste("    ", nms.sas, dollarsign,"\n"), collapse = "")
	input.stmt <- paste( "INPUT\n", input.stmt, ";\n")
	if (debug) pr(input.stmt)
	
	code <- paste("filename r2sas \'",datafile,"\';\n",   # might have to convert to backslashes
		"libname to \'G:/SAS\';\n",
		"data to.r2sas;\n",
		"infile r2sas delimiter=\';\' dsd LRECL =", lrecl+100, ";\n", sep = "")
	code <- paste(code,length.stmt, input.stmt, "\nrun;\n")
	if (debug) pr(code)
	if(!debug) cat(code, file = codefile)
	invisible(0)
}


## date()
## write.sas(dd[dd$wave > 1,setdiff(sort(names(dd)),c('qday','bday')) ]) # approx 1 min.
## date()
	 
if(F) { # current version in /R/coursefun.R
td <- function( basecol = NULL, col = c(3,5,4,6,7,8,2), lty = 1:7,
	lwd = 1, pch = 1:7, cex = 0.8, font = 1, len = 7, long = FALSE,
	new = FALSE, record = TRUE,  theme = col.whitebg,
	col.symbol = col, col.line = col,
	colsets = c( "plot.symbol","plot.line", "dot.symbol","dot.line","cloud.3d","box.dot"),...) {
	require(lattice)
	if (long) {
		col <- c(3,5,4,6,8,2)
		len <- 42
	}
	if (new) trellis.device( theme = theme, record = record, new = new, ...)
	len <- max( len, length(col.symbol), length(col.line), length( lty), length(lwd), length(pch),
		length(cex), length(font))
	spl <- trellis.par.get("superpose.line")
	spl$lty <- if ( is.null(lty)) spl$lty else lty
	spl$col <- if ( is.null(col.line)) spl$col else col.line
	spl$lwd <- if ( is.null(lwd)) spl$lwd else lwd
	trellis.par.set("superpose.line", Rows(spl, 1:len))
	sps <- trellis.par.get("superpose.symbol")
	sps$pch <- if ( is.null(pch)) sps$pch else pch
	sps$col <- if ( is.null(col.symbol)) sps$col else col.symbol
	sps$cex <- if ( is.null(cex)) sps$cex else cex
	sps$font <- if ( is.null(font)) sps$font else font
	trellis.par.set("superpose.symbol", Rows(sps, 1:len))
	if( !is.null(basecol) ) {
		for(ii in colsets) {
			tt <- trellis.par.get(ii)
			tt$col <- basecol
			trellis.par.set(ii,tt)
		}
	}
	invisible( attr( .Device, "trellis.settings"))
}
} # end of F


cap1 <- function(x, tofactor = is.factor(x),
		stop=c(" The"," Of"," By"," To"," And"),
		blanks=c(" ","(","\"","/","+")) {
	# capitalizes first letters 
	under2blank <- T
	if ( is.factor(x)) {
		ret <- cap1(levels(x))
		if ( length(unique(ret)) != length(ret)) warning("factor levels have been shortened")
		levels(x) <- ret
		return(x)
	}
	ret <- as.character(x)
	for ( ii in 1:length(ret)) {
		z <- ret[[ii]]
		if(under2blank) z <- gsub("_"," ",z)
		n <- nchar(z)
		z <- substring( z, 1:n, 1:n)
		zu <- toupper(z)
		zl <- tolower(z)
		zb <- c(" ",zu[-n])
		z <- 	paste( ifelse(zb %in% blanks, zu, zl), collapse ="")	
		for ( ss in stop) z <- gsub(ss,tolower(ss), z)
		ret[[ii]] <- z
	}
	ret
}

tran <- function( from, to, x, tofactor = is.factor(x)) {
	# This should return something like 'to' whenever all of 'x'
	# is in 'from'
	# Otherwise, if some of x is not modified then it needs to
	# either to change those to NAs or to leave them as is.
	# For a strict function, we can simply index to which will
	# leave it unchanged
	# which conflicts with my versions of tr, ergo this is called tran
	if(is.factor(from)) from <- as.character(from)
	if(is.factor(to)) to <- as.character(to)
	to <- rep(to, length=length(from))
	ret <- x
	if( is.factor(x) ) {
		ret <- as.character(ret)
		levs <- levels(x)
	}
	ret <- c(to,unique(ret)) [ match( ret, c(from, unique(ret)))]
	# DEBUG: print(ret)
	if (tofactor) {
		if(is.factor(x)) tolevs <- tran(from,to,levs)
		else tolevs <- to
		tolevs <- c(tolevs,unique(ret))
		tolevs <- intersect( tolevs, unique( ret ))
		ret <- factor(ret,levels = tolevs)
	}
    ret
}

tr <- function( x, from , to ) tran( from, to, x)

map <- function( from, to, x ) to[ match( x, from) ]

if ( FALSE){
    zf <- c('A',"FromNA",'Z','B',"N")
    zt <- factor(c('a',NA,'z','b',NA),levels=c('z','a','b',NA))
    zx <- c("Z","B",NA,"N","M")
    
    map( zf, zt, zx)
}



abind.rdc <- function( arr1, arr2, d, facename = "") {
	# glue arr1 to arr2 along dimension d (i.e. face 'not d')
	# copied from library gm 05 05 03
	d1 <- dim(arr1)
	n1 <- length(d1)	
	d2 <- dim(arr2)
	n2 <- length(d2)
	dn1 <- dimnames( arr1 )
	dn2 <- dimnames( arr2 )
	
	arenull <- is.null(dn1) & is.null(dn2)
	if ( is.null(dn1)) {
		dn1 <- lapply( as.list(d1), function(x) seq(1,x))
		dimnames(arr1) <- dn1
	}
	if ( n1 != n2 ) {
		d2 <- d1
		d2[d] <- 1
		dn2 <- dn1
		dn2[[d]] <- facename
		dimnames(arr2) <- NULL
		dim(arr2) <- d2
		dimnames(arr2) <- dn2
		n2 <- n1
	}
	if ( is.null(dn2)) {
		dn2 <- lapply( as.list(d2), function(x) seq(1,x))
		dimnames(arr2) <- dn2
	}
	perm <- 1:n1
	perm[c(d,n1)] <- c(n1,d)  # perm is an involution
	
	arr.p1 <- aperm(  arr1, perm )

	arr.p2 <- aperm(  arr2, perm )
	dret <- d1[perm]
	dret[n1] <- dret[n1] + d2[d]
	dnret <- dn1
	dnret[[d]] <- c(dnret[[d]],dn2[[d]])
	ret <- c(arr.p1, arr.p2)
	dim(ret) <- dret

	ret <- aperm( ret, perm )
	dimnames( ret ) <- dnret
	ret
}



# acond( with(dds, table( lang, prov)),2)

atotal.rdc <- function( arr, FUN = sum, name = "Total",...) {
	# copied from library gm  05 05 03
	# 05 05 03: added option for name
	d <- dim(arr)
	if ( length(d) == 1) {
		arr <- c(arr)
		d <- dim(arr)
	}
	if ( is.character( FUN ) ) FUN <- get(FUN,mode='function')
	else if( mode(FUN) != "function") {
		farg <- substitute(FUN)
		if (mode(farg) == 'name' ) FUN <- get(farg,mode='function')
		else stop(paste("\"", farg, "\" is not a function", sep=""))
	}
	if ( is.null(d) ) {
		ret <- c(arr, FUN(arr))
		names(ret)[length(ret)] = name
		return(ret)
	}
	n <- length (d)
	name <- rep(name, length= n)
	ret <- arr
	ind <- 1:n
	for ( i in n:1 ) {
		new <- apply(ret, ind[ -i ], FUN,...)
		ret <- abind( ret, new, i, name[i])
	}
	ret
}

# aprop and acond moved to tab.R

###
###  Utility functions for regression
###

    ##
    ##  na.fitted na.resid   pads with NA to fit original data frame
    ##
    ##  extend to new classes by writing a 'na.pad' methoc


    na.fitted <- function(fit,...) na.pad( fit, fitted(fit,...))
    na.resid <- function(fit,...) na.pad( fit, resid(fit,...))
    na.residuals <- function(fit,...) na.pad( fit, residuals(fit,...))
    # na.predict <- function(fit,...) na.pad( fit, predict(x,...))

    na.pad <- function(fit,...) UseMethod('na.pad')
    na.pad.lme <- function( x, obj) {
        ind <- rep(NA, nrow(x$data))
        names(ind) <- rn <- rownames(x$data)

        ind.part <- 1:nrow(x$resid)
        names(ind.part) <- rownames(x$resid)
        ind[ names(ind.part) ] <- ind.part
        # disp(ind)
        if( !is.null( dim (obj)) ) ret <- obj[ ind,]
        else ret <- obj[ind]
        names(ret) <- rn
        ret
    }

    na.pad.default <- function(x,obj) {
        stop(paste("You need to write a method for na.pad for objects of class",class(x),"\n",
            "see na.pad.lme for an example"))
    }





# Lagging
#    - interesting if 'index' is continuous and/or if there
#      are missing waves.
#    - need some kind of intrapolation to compute a derivative
#      at the time of observation

# The first function here only works well with complete data
# and an integer index ranging from 1 to n. (n can vary from 
# subject to subject)


Lag.0 <- function(x,id,idx,lag = 1,at= NULL, check=T) {
	# computes Lagged values but without intrapolation
	# adds 'at'

	if (check) {
		comb <- paste(id,idx)
		if (any(duplicated(comb))) stop("id not unique in each level of idx")
	}
	if(any(is.na(idx))) stop("NAs in idx")
	if(any( round(idx) != idx)) stop("Non integers in idx")
	ret <- x
	id <- as.character(id)
	names(x) <- id
	for ( i in max(idx):min(idx)) {
		to.pos <- idx == i
		if( is.null(at) ) from <- x[idx == (i-lag)]
		else from <- x[idx == at ]
		ids <- names( x[to.pos])
		ret[to.pos] <- from[ids]
	}
	ret
}


cLag <- function(x, id = rep(1,length(x)), time = 1:length(x), lag=1, at = NULL, check = T, idx) {
   # renamed cLag to avoid conflict with Hmisc::Lag
   ## Takes approx. 30% of time of version Lag.0 on a small problem
  if (check) {
		comb <- paste(id,time)
		if (any(duplicated(comb))) stop("id not unique in each level of time")
	}
	if( !missing(idx)) time = idx    # for compatibility with previous argument names
  ret <- x
  names(ret) <- paste(id,time,sep='::')
  retnn <- if(is.null(at)) paste(id,time - lag,sep='::')  else paste(id,at,sep="::")
  ret [ retnn ]
}

Lag <- cLag # historical name that conflicts with Hmisc

 if(F) {
## small test of Lag
zx <- c(2,3,2,4,2,4, 5,3,4,5,7,8,9)
zid<- c(1,1,2,2,3,3 ,3,4,4,4,4,4,4)
cbind( zid, zx, Lag(zx,zid,zx))
}


cLagI <- function(x,id,time,lag=1,delta=.01,check=T) {
    # renamed from LagI to avoid conflict with Hmisc
	# lags by intrapolating 
	# with complete data at each value of index, this does the same thing
	# as Lag
	# If values of Lag are skipped then we get linear intrapolations.
	# Note that 'delta' should be small enough so that values of x are
	# at least delta apart. However, too small a value for delta introduces
	# numerical error
	#  	
	if (check) {
		comb <- paste(id,time)
		if (any(duplicated(comb))) stop("id not unique in each level of idx")
	}
	ret <- x
	id <- as.character(id)
	names(x) <- id
	names(time) <- id
	for (nn in unique(id)){
		pos <- id == nn
		xx <- x[pos]
		tt <- time[pos]
		topred <- tt-delta
		drop <- is.na(xx)|is.na(tt)
		xxc <- xx[!drop]
		ttc <- tt[!drop]
		nl <- length(xxc)
		if ( nl > 0) {
			if ( nl > 1 ) xx.p <- approx(ttc,xxc,topred)$y
			else xx.p <- NA 
			xx.lag <- xx - lag*(xx - xx.p)/delta
			ret[pos] <- xx.lag
		}	
	}
	ret
}

LagI <- cLagI

cDiffI <- function(xx,...) {
	xx - LagI(xx,...)
}

DiffI <- cDiffI

# tt <- c(1,2,3,1,2,3,1,2,3)
# xx <- c(1,2,3,1,2,3,1,2,3)
# id <- rep(letters[1:3],rep(3,3))
# xx-LagI(xx,id,tt)
# DiffI(xx,id,tt,delta=.1)

###
###  Splines and polynomial: quadratic and cubic
###


qs <- function(x, knots=quantile(x,pc), exclude = 0, pc = c(.25,.75)) {
	# simple quadratic spline
	ret <- cbind(x,x^2)
	nam <- deparse(substitute(x))
    if ( missing(knots) ) warning("May be unsafe for prediction with newdata")
	for ( kk in knots ) {
		z <- x - kk
		z [z<0] <- 0
		ret <- cbind(ret, z^2)
	}

	if ( is.null(knots) ) dimnames(ret)[[2]] <- paste(nam,c('','^2'), sep='')
	else dimnames(ret)[[2]] <- paste(nam,c('','^2',paste('.',knots,sep='')),sep='')
	if (exclude > 0) ret <- ret[, -( 1:exclude)]
	ret

	ret
}



lsp <- function(x, knots=quantile(x,pc), exclude = 0, pc = c(.25,.75)) {
	# linear spline
	ret <- cbind(x)
	nam <- deparse(substitute(x))	
    if ( missing(knots) ) warning("May be unsafe for prediction with newdata")
	for ( kk in knots ) {
		z <- x - kk
		z [z<0] <- 0
		ret <- cbind(ret, z)
	}
	if ( is.null(knots) ) dimnames(ret)[[2]] <- paste(nam,c(''), sep='')
	else dimnames(ret)[[2]] <- paste(nam,c('',paste('.',knots,sep='')),sep='')
	if (exclude > 0) ret <- ret[, -( 1:exclude)]
	ret

	ret
}

cs <- function(x, knots=quantile(x,pc), exclude = 0, pc = c(.25,.75))  {
	# simple cubic spline
    if ( missing(knots) ) warning("May be unsafe for prediction with newdata")
	ret <- cbind(x,(x)^2,(x)^3)
	nam <- deparse(substitute(x))
	for ( kk in knots ) {
		z <- (x) - kk
		z [z<0] <- 0
		ret <- cbind(ret, z^3)
	}
	if ( is.null(knots) ) dimnames(ret)[[2]] <- paste(nam,c('','^2','^3'), sep='')
	else dimnames(ret)[[2]] <- paste(nam,c('','^2','^3',paste('.',knots,sep='')),sep='')
	if (exclude > 0) ret <- ret[, -( 1:exclude)]
	ret
}

Poly <- function(x, order=2, exclude = 0)  {
	# polynomial
	# exclude: NULL to include intercept, otherwise exclude order up to 
	# including exclude
	ret <- cbind(rep(1,length(x)))
	for ( i in 1:order) ret <- cbind(ret, x^i)
	nam <- deparse(substitute(x))
	powers <- paste(nam,"^",0:order,sep="")
	powers[1] <- "Intercept"
	powers[2] <- nam
	dimnames(ret)[[2]] <- powers
	if (!is.null(exclude)) ret <- ret[, -( 1:(exclude+1))]
	ret
}



if (F) {
	## shows that columns of cs span same space as columns of bs	

	zz <- 0:20
	cs(zz, c(5,15))

	#search()
	#  detach(9)
	### comparison with ns

	zd <- data.frame( x = 0:20, y = (-10:10)^3 + rnorm(21))

	bss <- bs(zd$x,knots = c(5,15))
	css <- cs(zd$x,knots = c(5,15))

	fit <- lm(bss ~ css - 1)
	round(coef(fit),7)
	summary(fit)
	
	fit <- lm(css ~ bss - 1)
	summary(fit)
	round(coef(fit),7)

	fit <- lm(y ~ 1+cs(x,c(5,15)), zd)
	summary(fit)
	anova(fit)
	
	fit <- lm(y ~ x + I(x^2) + cs(x,c(5,15),2), zd)
	summary(fit)
	anova(fit)

	fit <- lm(y ~ x + I(x^2) + I(x^3) + cs(x,c(5,15),3), zd)
	summary(fit)
	anova(fit)
	
	fit <- lm(y ~ x + I(x^2) + cs(x,pc=c(.05,.95)), zd, singular.ok=T)
	summary(fit)
	anova(fit)
	
	qqr <- qr.Q(qr(cs(0:20,c(5,15))))
	fit <- lm( css ~ qqr-1)
	summary(fit)
	round(coef(fit),5) == round(qr.R(qr(cs(0:20, c(5,15)))),5)

}


## much better approach to stack:

mergec <- function( ... ,vname = '.type',oname = ".order") {
	
	## stacks data frames adding a new variable .type to identify each source data frame
	## good way to combine data and prediction data frames to show data and fitted
	## values in a common plot
	
	z <- list(...)
	for ( ii in 1:length(z)) {
		z[[ii]][,vname] <- rep(ii, nr <- nrow( z[[ii]]))
        z[[ii]][,oname] <- 1:nr
	}
	ret <- z[[1]]
	for ( ii in 2:length(z)) ret <- merge(ret,z[[ii]], all = TRUE)
	ret [ order( ret[[vname]], ret[[oname]]),]
}

if(F) {
	## test mergec
	zd1 <- data.frame(x=1:4, y = letters[1:4], z = 1:4,d=LETTERS[1:4])
	zd2 <- data.frame(x=5:7, y = letters[5:7])
	zd3 <- data.frame(x=8:9,  w=1:2, d=1:2)
	sapply(zz <- mergec(zd1,zd2,zd3),class)
	zz
	levels(zz$d)
}

###
###  extended merge for typical merging of data frames by 'ID' with possible
###  other common variables
###
###

xmerge <- function(x, y, by , all = TRUE, dropdots = FALSE , verbose = FALSE, debug = TRUE, from = FALSE, ... ) {
    help <-
"This is a modification of merge that combines consistent variables
even if not specified in 'by' to keep a common name.
-- Some errors fixed Apr 24, 2007"

    xm <- function( a, b, tofac = is.factor(a)||is.factor(b)) {
          # replace values of a with non-missing values of b
          if ( tofac ) {
              levs <- union( levels(b), levels(a))
              a <- as.character(a)
              b <- as.character(b)
          }
          b [ is.na(b) ] <- a[is.na(b)]

          if ( tofac ) {
				levs <- union( levs, unique(b))
				b <- factor(b,levels = levs)
			}
          b
    }
    na2f <- function(x) {
        x[is.na(x)] <- F
        x
    }
    consistent <- function(a,b) {
            # which values are consistent if neither is NA
            if( is.factor(a)) a <- as.character(a)
            if( is.factor(b)) b <- as.character(b)
            !na2f(a != b)
    }
    if ( from ) {
        xname <- deparse(substitute(x))
        yname <- deparse(substitute(y))
        # x[[paste("F.",xname,sep="")]] <- rep(T, nrow(x))
        # y[[paste("F.",yname,sep="")]] <- rep(T, nrow(y))
        x[[".F"]] <- rep("x", nrow(x))
        y[[".F"]] <- rep("y", nrow(y))
    }
	## ids in each file
	xby <- x[,by,drop=F]
	yby <- y[,by,drop=F]
	xby$.file <- rep('x',nrow(xby))
	yby$.file <- rep('y',nrow(yby))
	by2 <- rbind( xby, yby)
	if ( verbose ) cat("\nby in x and y:\n")
	if ( verbose ) print(atotal(do.call("tab",by2),sum,"Total"))
    nams <- union( names(x), names(y))
   	if ( verbose ) print( c( DimX = dim(x), DimY = dim(y)))
    if ( verbose ) cat("\nVariables in both:\n")
    if ( verbose ) print(intersect( names(x), names(y)))
    if ( verbose ) cat("\nVariables in X only:\n")
    if ( verbose ) print( setdiff( names(x), names(y)))
    if ( verbose ) cat("\nVariables in Y only:\n")
    if ( verbose ) print( setdiff( names(y), names(x)))
    # compare two data frames, e.g. for merging
    x$FromX <-  1:nrow(x)
    y$FromY <- 1:nrow(y)
    mm <- merge( x, y, by, all = TRUE, ...)
    # names in both data frames
    newroots <- setdiff( intersect(names(x),names(y)), by)

    if ( verbose ) cat("\nDimension of merged data frames:\n")
    if ( verbose ) print(c("DimMerge"=dim(mm)))
    if ( verbose ) cat("\nNames of variables in merged data frame:\n")
    if ( verbose ) print(names(mm))

    # Are similarly named variables consistent

    if(F){
        dotx <- grep("\\.x", names(mm), value =T)
        if ( verbose )print( c(dotx = dotx))
        doty <- grep("\\.y", names(mm),value = TRUE)
        if ( verbose )print( c(doty = doty))
        rootx <- substring( dotx, 1, nchar(dotx) - 2 )
        rooty <- substring( doty, 1, nchar(doty) - 2 )
        newroots <- intersect( rootx, rooty )
    }
    FromBoth <- !is.na(mm$FromX) & !is.na(mm$FromY)
    Xonly <-  !is.na(mm$FromX) & is.na(mm$FromY)
    Yonly <-  is.na(mm$FromX) & !is.na(mm$FromY)
    if ( verbose )cat("\nRows in:\n")
    if ( verbose )print( c( Both = sum(FromBoth), Xonly = sum(Xonly), Yonly = sum(Yonly)))
    if ( verbose )cat("\nThe following variables occur in both data frames:\n")
    if ( verbose )print( newroots )
    drop.list <- character(0)
    for ( nn in newroots) {
        nn.x <- paste(nn,".x", sep = '')
        nn.y <- paste(nn,".y", sep = '')
        mm[[nn]] <- xm( mm[[nn.x]], mm[[nn.y]])
        # Note that consistent NAs should be okay here
        if( all( same <- consistent(mm[[nn.x]],mm[[nn.y]] )) ) {
            if ( verbose )cat ("Variable ",nn," is consistent\n")
            drop.list <- c(drop.list, nn)
        } else {
            if ( verbose )cat("Variable ",nn," is inconsistent in the following rows:\n")
            if ( verbose )print(mm[ same, c(by, nn.x, nn.y, nn)])
        }
    }
	if( dropdots ) drop.list <- newroots
    drop <- if( length( drop.list)>0) {
          c(paste(drop.list, "x", sep = "."), paste(drop.list, "y", sep = ".") )
          } else character(0)
	if ( verbose )cat("\nDrop list:\n")
	if ( verbose )print( drop)
    if( length(drop) > 0) {
            if ( verbose )print( c(drop=drop))
            mm <- mm[, - match( drop, names(mm))]
    }
    # nice ordering
    onams <- 1:length(nams)
    #print(onams)
    onams <- c(onams, onams+.1, onams +.2)
    #print(onams)
    names(onams) <- c(nams, paste(nams,".x",sep=''), paste(nams,".y",sep =''))
    keep <- intersect( names(sort(onams)), names(mm))
    mm[, keep]
}


if ( F ) {  # test xmerge
    d1 <- data.frame(id = 1:5, x = 1:5, d = 1:5, v1 = 1:5)
    d2 <- data.frame(id = 3:7, x = 3:7, d = 1:5, v2 = 3:7)

    # need to specify 'id'

    xmerge(d1,d2,by='id',verbose=T) # right row but inconsistent variable 'd' is created with 2nd df replacing first -- no warning
    xmerge(d1,d2,by='id',all=T)   # ditto
}


part <- function( A, B) {
   # partition two sets
	if ( is.factor(A)) A <- as.character(A)
	if (is.factor(B)) B <- as.character(B)
	ret <- list( "A and B" = intersect(A,B),
        "A - B" = setdiff( A, B),
        "B - A" = setdiff(B, A))
	attr(ret,"N") <- sapply(ret, length)
    ret
}


###
###  Rbind to stack data and prediction frames
###

Rbind <- function(x, ...) UseMethod("Rbind")

Rbind.data.frame <-
function( ... ,vname = '.which',oname = ".order") {

        ## stacks data frames adding a new variable .which to identify each source data frame
        ## good way to combine data and prediction data frames to show data and fitted
        ## values in a common plot
        Unique <- function(x, ...) UseMethod("Unique")
        Unique.factor <- function( x, ...) levels(x)
        Unique.default <- function( x, ...) unique(x)
        Maxp1 <- function( x ) {
              mm <- max(c(0, as.numeric( as.character(x))), na.rm = TRUE)
              ceiling( mm + 1)
        }
        z <- list(...)
        cum.which <- numeric(0)
        for ( ii in 1:length(z)) {
                ddi <- as.data.frame(z[[ii]],stringsAsFactors = FALSE )
                if ( vname %in% names(ddi)) cum.which <- c(cum.which, Unique(ddi[[vname]]))
                else {
                     ddi[[vname]] <- rep(vv <- Maxp1(cum.which), nrow( ddi))
                     cum.which <- c(cum.which, vv)
                }
                ddi[,oname] <- 1:nrow(ddi)
                z[[ii]] <- ddi
        }
        ret <- z[[1]]
        if( length(z) > 1 )for ( ii in 2:length(z)) ret <- merge(ret,z[[ii]], all = TRUE, sort = FALSE)
        ret [ order( ret[[vname]], ret[[oname]]),]
}

Rbind.list <- function(x,...) {
    ret <- Rbind.data.frame(x, ...)
    as.list(ret)
}




if (FALSE) {
    df1 <- data.frame( a = c('z','y','b'), x = 1:3, y1 = 1:3)
    df2 <- data.frame( a = c('a','d','b'), x = c("A","A","B"), y2 = 4:6)

    mat1 <- cbind( x=1:3, y1 = 1:3)

    Rbind( mat1, df2)

    df1
    df2
    z <- Rbind( df1, df2)
    z
    levels(z$a)
    levels(z$x)

    list1 <- list( a = 1:4, b = letters[1:4])
    list2 <- list( a = 5:6, b = letters[5:6], c = factor(c("A","B")))
    Rbind(list1, list2)

}






###
###   Reshaping data: function to call 'reshape'
###

long <- function ( data , varying=NULL , sep.varying = "\\.", v.names = names(Varying),
                 timevar = 'time', idvar = 'id', ids = 1:NROW(data),
                 times = seq(length=length(Varying[[1]])),
                 drop = NULL, new.row.names = NULL,
                 split = list(regexp="\\.", include = FALSE),
                 debug = FALSE){
    help <- "
    Example:   long(dwl, varying = c('wmrss', 'lm1tot','lm2tot'), sep='')
    "
           nn <- names(data)
           if ( debug ) disp(nn)
           Varying <- list()
           for ( ii in 1:length(varying)){
               if(length( prefix <- varying[[ii]]) == 1) {
                          ns <- grep(paste("^",prefix,sep.varying,sep=""), nn,value = T)
                          if (debug) disp(ns)
                          Varying[[ii]] <- ns
                          names(Varying)[ii] <- prefix
               }  else {
                    Varying[[ii]] <- varying[[ii]]
                    names(Varying)[ii] <- names(varying)[ii]
               }

           } 
           #print( varying) 
           #print( v.names )
           if ( debug ) disp(Varying)
           if ( debug ) disp(times)
         ret <- stats::reshape( data, Varying, v.names , timevar, idvar, ids,
                times , drop, direction = 'long', new.row.names,
                split)
         ret [ order( ret[[idvar]], ret[[timevar]]),]
         
} 

#zd <- data.frame( sub = c('a','b','c'), x.1 = 10+1:3, x.2 = 20+1:3, y.1 = c('a','b','c'), y.2 = factor(2:4))
#long( zd, varying = list("y","x"))
#long(zd, varying = list( y = c('y.1','y.2'), x = c("x.1","x.2")))




##
##  constant: to check whether something is constant
##

constant <- function(x,...) UseMethod("constant")

constant.default <- function(x, na.rm = FALSE,...) {
    if (na.rm) x <- na.omit(x)
    length(unique(x)) <= 1
}


constant.data.frame <- function( x, id  , all = FALSE , ...) {
   ## Description:    (G. Monette, June 13, 2005)
   ## check if values of variables are constant within levels of id
   ## if no id, check each variable
   ## if all == TRUE, report overall instead of within levels of id
   ## Possible improvements:
   ##    allow nested formulas: ~id1/id2 for id and report level of
   ##    each variable
   ##    [see varLevel()]
   
   ## note that the following code allows id to be given as a name or as
   ## a formula
   
   if (missing(id)) ret <- sapply(x, constant,...)
   else {
        id <- eval(substitute(id), x, parent.frame())
        if( inherits(id,"formula") ) id <- c( model.frame(id,x) )
        ret <- sapply(x, function(xx) tapply(xx, id, constant, ...))
        if ( all ) ret <- apply( ret, 2, all)
   }
   ret
}

varLevel <- function(x, form, ...) {
         ## Description:    (G. Monette, June 13, 2005)
         ## shows levels of each variable with respect to grouping formula
         ## of form ~id or nested ids ~id1/id2
         ## Level 0 is a constant for the whole data frame
         ## Level <= 1 implies the variable is constant within levels of id1 
         ## Level <= 2 implies the variable is constant within levels of id2 
         ##    ... etc.
         ## NOTE: NA counts as a distinct value
        sel <- model.frame( form, x )
        z <- list()
        idx <- rep('',nrow(x))
        z[[1]] <- constant(x,...)
        for ( ii in 1:ncol(sel)) {
            idx <- paste(idx, as.character(sel[[ii]]), sep = ";")
            z[[ii+1]] <- constant( x, idx, all = TRUE,...)
        }
        # print(z)
        ret <- do.call("rbind", z)
        # print(ret)
        ret <- length(z) - apply( ret*1, 2 , sum)
        ret
}

if (FALSE) {
        up <- function( dd, form , all = FALSE, keep = ncol(sel) ) {
               ##
               ## Replaced Nov 11, 2007  See below

              ## Description:    (G. Monette, June 10, 2006)
              ## Creates a higher level data set by selecting first
              ## after ordering and keeping only invariant variables
              ## To have summary variables they need to have been created
              ## with, e.g., capply
              sel <- model.frame( form , dd , na.action = na.include )
              if( !all) vl <- varLevel( dd, form ) else vl <- rep(0,ncol(dd))
              dd [ na.omit(sapply( split( 1:nrow(dd), sel ), function(x) x[1])), vl < keep + 1, drop = FALSE]
        }


        up <- function( dd, form , all = FALSE, keep = ncol(sel) ) {
              ## Description:    (G. Monette, June 10, 2006)
              ## Creates a higher level data set by selecting first
              ## after ordering and keeping only invariant variables
              ## To have summary variables they need to have been created
              ## with, e.g., capply
              #sel <- model.frame( form , dd , na.action = na.include )
              #if( !all) vl <- varLevel( dd, form ) else vl <- rep(0,ncol(dd))
              #dd [ na.omit(sapply( split( 1:nrow(dd), sel ), function(x) x[1])), vl < keep + 1, drop = FALSE]
              require(nlme)

              gsummary( dd, form = form, invariantsOnly = ! all )
        }
} # end of FALSE


up <-
function ( object, form = formula(object),
         all = FALSE,
         FUN = function(x) mean(x, na.rm = TRUE),
         omitGroupingFactor = FALSE,
         groups, invariantsOnly = !all , ...)
{
help <- "

     object     data frame
     form       one-sided formula defining groups, e.g. ~school/student
     all        if FALSE keep only variables that are invariant within
                clusters. If TRUE, summarize all variables, using mean or Mode
                
     Examples:

          up( dd, ~ id, all = TRUE)
          
     To plot summary panel:

          ddu <- up( dd, ~ id, all = TRUE)
          ddu$id <- 'summary'
          xyplot( Y ~ X | id, rbind(dd,ddu),
                  panel = function( x, y, ...) {
                        panel.xyplot( x, y, ...)
                        panel.lmline( x, y, ...)
                        
                  }
          )
                
     Created from nlme::gsummary, modified for two reasons
      1) so gsummary would be available when using lme4
      2) modified to handle multiple levels which you definitely
         need when making successive summaries

      Possible changes:
          create attribute that identifies extent to which a variable is 'varying'

"


    if (!inherits(object, "data.frame")) {
        stop("Object must inherit from data.frame")
    }
    sel.mf <- model.frame( form , object , na.action = na.include )
    if ( ncol(sel.mf) > 1) {
       sel <- apply( sel.mf, 1 , paste, collapse = "/")
       groups <- as.factor(sel)
    } else {
       groups <- as.factor(sel.mf[[1]])
    }
    gunique <- unique(groups)
    firstInGroup <- match(gunique, groups)
    asFirst <- firstInGroup[match(groups, gunique)]
    value <- as.data.frame(object[firstInGroup, , drop = FALSE])
    row.names(value) <- as.character(gunique)
    value <- value[as.character(sort(gunique)), , drop = FALSE]
    varying <- unlist(lapply(object, function(column, frst) {
        aux <- column
        if( is.matrix( aux))aux[] <- as.character( aux )
        else aux <- as.character(aux)
        if ( is.matrix( aux )) any(!identical( aux, aux[frst,]))
        else any(!identical(aux, aux[frst]))
    }, frst = asFirst))
    if (any(varying) && (!invariantsOnly)) {
        Mode <- function(x) {
            aux <- table(x)
            names(aux)[match(max(aux), aux)]
        }
        if (data.class(FUN) == "function") {
            FUN <- list(numeric = FUN, ordered = Mode, factor = Mode)
        }
        else {
            if (!(is.list(FUN) && all(sapply(FUN, data.class) ==
                "function"))) {
                stop("FUN can only be a function or a list of functions")
            }
            auxFUN <- list(numeric = mean, ordered = Mode, factor = Mode)
            aux <- names(auxFUN)[is.na(match(names(auxFUN), names(FUN)))]
            if (length(aux) > 0)
                FUN[aux] <- auxFUN[aux]
        }
        for (nm in names(object)[varying]) {
            dClass <- if (is.ordered(object[[nm]]))
                "ordered"
            else if (is.factor(object[[nm]]))
                "factor"
            else mode(object[[nm]])
            if (dClass == "numeric") {
               if( is.matrix ( object[[nm]])){
                 zmat <- object[[nm]]
                 ret <- list()
                 for ( jj in seq_len(ncol(zmat))) {
                     ret[[jj]] <- as.vector( tapply( zmat[,jj],
                               groups, FUN[['numeric']],...))
                 }
                 value[[nm]] <- do.call(cbind, ret)

               } else {
                value[[nm]] <- as.vector(tapply(object[[nm]],
                  groups, FUN[["numeric"]], ...))
                }
            }
            else {
                value[[nm]] <- as.vector(tapply(as.character(object[[nm]]),
                  groups, FUN[[dClass]]))
                if (inherits(object[, nm], "ordered")) {
                  value[[nm]] <- ordered(value[, nm], levels = levels(object[,
                    nm]))[drop = TRUE]
                }
                else {
                  value[[nm]] <- factor(value[, nm], levels = levels(object[,
                    nm]))[drop = TRUE]
                }
            }
        }
    }
    else {
        value <- value[, !varying, drop = FALSE]
    }
    if (omitGroupingFactor) {
        if (is.null(form)) {
            stop("Cannot omit grouping factor without \"form\"")
        }
        grpForm <- getGroupsFormula(form, asList = TRUE)
        if (missing(level))
            level <- length(grpForm)
        grpNames <- names(grpForm)[level]
        whichKeep <- is.na(match(names(value), grpNames))
        if (any(whichKeep)) {
            value <- value[, whichKeep, drop = FALSE]
        }
        else {
            return(NULL)
        }
    }
    value
}




###
###  Row.names heading on data.frame
###


print.data.frame.lab <-    
    function (x, ..., digits = NULL, quote = FALSE, right = TRUE) 
    {
        labs <- attributes(x)$labs
        if (length(x) == 0) {
            cat("NULL data frame with", length(row.names(x)), "rows\n")
        }
        else if (length(row.names(x)) == 0) {
            print.default(names(x), quote = FALSE)
            cat("<0 rows> (or 0-length row.names)\n")
        }
        else {
            mat <- as.matrix(format.data.frame(x, digits = digits, 
                na.encode = FALSE))
            labs <- c(labs,"","")
            labs <- labs[1:2]
            names(dimnames(mat)) <- labs    
            print(mat , ..., quote = quote, right = right)
        }
        invisible(x)
    }

"[.data.frame.lab" <- function(x, ...){
    lab <- labs(x)
    ret <- get("[.data.frame")(x,...)
    if( inherits(ret, "data.frame")) labs(ret) <- lab
    ret
} 



"labs<-" <- function(x,...) UseMethod("labs<-")

"labs<-.data.frame" <- function( x, value ) {
      value <- c( value, "", "") [ 1:2 ]
      attr(x,"labs") <- value
      if( !inherits(x,"data.frame.lab")) class(x) <- c( "data.frame.lab", class(x))
      x
  } 
"labs<-.default" <- function(x, value) {
      nd <- length(dim(x))
      value <- c( value, rep("",nd))[1:nd]
      names(dimnames(x)) <- value
      x
}

labs <- function(x,...) UseMethod("labs")
labs.data.frame.lab <- function( x ,...) attr(x,"labs")
labs.default <- function(x,...) names(dimnames(x))



###
###  Trellis additions
###


###
###   Miscellaneous utility functions
###


val2na <- function( x, val) {
	## val2na(1:10, c(3,5))
	## val2na(factor(c('a','b','b','c','a')), c('b','a'))
	x[match(x,val,0)>0] <- NA
	x
}

# zf <- factor(c('a','a','b',NA,'c','a',NA))


grepv <- function(...) grep( ..., value = TRUE)

# grepl <- function(pattern,x,...) match( 1:length(x) , grep(pattern,x,...), 0) >0

p <- function(...) paste(..., sep ="")

ch <- function(x) as.character(x)

"%less%" <- function( a, b) setdiff(a,b)
"%or%" <- function( a, b) union(a,b)
"%and%" <- function( a, b) intersect( a, b)


lib <- function(x) {
    xn <- deparse(substitute(x))
    if ( !do.call("require", list(as.name(xn))) ) {
        install.packages(xn)
        do.call("library",list(as.name(xn)))
    }
}

###
###  Interfacing with SAS
###


sasin <- function(file, tfile = tempfile() ) {
# moved to "/R/fun.R"
    help = "
    sasin reads a .csv file created by SAS with
       ODS CSV FILE = 'file';
        < SAS procedure statements >
       ODS CSV CLOSE;
    The tables produced by SAS are elements in the list
    returned by sasin.
    "

    todf <- function(ll) {
        if ( length(ll) < 3) return (character(0))
        if ( length(ll) == 3) return (ll[2])
        cat(ll[2],"\n" , file = tfile)
        for ( ii in 3:(length(ll)-1)) {
            cat(ll[ii], "\n",file = tfile ,append = TRUE)
        }
        df <- read.csv( tfile , header = F)
        if ( !any ( sapply( df, is.numeric ))) df <- read.csv(tfile)
        df
    }
    readin <- scan(file,what="",sep="\n",blank.lines.skip=F)
    blanks <- which(readin == "")
    head.pos <- c(1,1+head(blanks,-1))
    heads  <- gsub('\\"|,',"",readin[head.pos])
    # disp(heads)
    reps   <- diff( c(head.pos, 1+length(readin)))
    # disp(reps)
    heads  <- rep(heads, reps)
    readin <- split( readin, heads)
    readin <- lapply( readin , todf)
    readin
}




#########
#########   SCRAPS
#########

###
### Stack
###

## NOTE THAT R HAS A FUNCTION CALLED stack
if (F) {
stack <- function(...) {
	# simple version of stack, uses names of first df.
      # 05-04-19: THis function might be obsolete d/t mergec below
	ll <- list(...)
	llclass <- sapply(ll,function(x) inherits(x,'data.frame'))
	if(!all(llclass)){
		dfr <- ll[[1]]
		keep <- ll[[2]]
		add <- ll[[3]]
		ll <- list(0)
		for ( i in 1:length(add) ) ll[[i]] <- dfr[,c(keep,add[i])]
	} 
	nam <- names(ll[[1]])
	for ( ii in 2:length(ll)) names(ll[[ii]]) <- nam
	ret <- do.call('rbind',ll)
	nrows <- sapply(ll, function(x) dim(x)[1])
	ret$Type. <- rep(1:length(ll),nrows)
	ret
}
}



#####
#####  anova.lme
#####

xanova.lme <- function (object, ..., test = TRUE, type = c("sequential", "marginal"),
    adjustSigma = TRUE, Terms, L, verbose = FALSE)
{
    warning("This is a modified version of anova.lme that uses min dfs for the denominator")
    Lmiss <- missing(L)
    dots <- list(...)
    if ((rt <- (length(dots) + 1)) == 1) {
        if (!inherits(object, "lme")) {
            stop("Object must inherit from class \"lme\" ")
        }
        vFix <- attr(object$fixDF, "varFixFact")
        if (object$method == "ML" && adjustSigma == TRUE) {
            vFix <- sqrt(object$dims$N/(object$dims$N - ncol(vFix))) *
                vFix
        }
        c0 <- solve(t(vFix), fixef(object))
        assign <- attr(object$fixDF, "assign")
        nTerms <- length(assign)
        if (missing(Terms) && Lmiss) {
            type <- match.arg(type)
            Fval <- Pval <- double(nTerms)
            nDF <- integer(nTerms)
            dDF <- object$fixDF$terms
            for (i in 1:nTerms) {
                nDF[i] <- length(assign[[i]])
                if (type == "sequential") {
                  c0i <- c0[assign[[i]]]
                }
                else {
                  c0i <- c(qr.qty(qr(vFix[, assign[[i]], drop = FALSE]),
                    c0))[1:nDF[i]]
                }
                Fval[i] <- sum(c0i^2)/nDF[i]
                Pval[i] <- 1 - pf(Fval[i], nDF[i], dDF[i])
            }
            aod <- data.frame(nDF, dDF, Fval, Pval)
            dimnames(aod) <- list(names(assign), c("numDF", "denDF",
                "F-value", "p-value"))
            attr(aod, "rt") <- rt
        }
        else {
            nX <- length(unlist(assign))
            if (Lmiss) {
                if (is.numeric(Terms) && all(Terms == as.integer(Terms))) {
                  if (min(Terms) < 1 || max(Terms) > nTerms) {
                    stop(paste("Terms must be between 1 and",
                      nTerms))
                  }
                }
                else {
                  if (is.character(Terms)) {
                    if (any(noMatch <- is.na(match(Terms, names(assign))))) {
                      stop(paste("Term(s)", paste(Terms[noMatch],
                        collapse = ", "), "not matched"))
                    }
                  }
                  else {
                    stop("Terms can only be integers or characters")
                  }
                }
                dDF <- unique(object$fixDF$terms[Terms])
                if (length(dDF) > 1) {
                  # stop("Terms must all have the same denominator DF")
                  warning("Terms do not all have the same denominator DF -- using the minimum")
                  dDF <- min(dDF)
                }
                lab <- paste("F-test for:", paste(names(assign[Terms]),
                  collapse = ", "), "\n")
                L <- diag(nX)[unlist(assign[Terms]), , drop = FALSE]
            }
            else {
                L <- as.matrix(L)
                if (ncol(L) == 1)
                  L <- t(L)
                nrowL <- nrow(L)
                ncolL <- ncol(L)
                if (ncol(L) > nX) {
                  stop(paste("L must have at most", nX, "columns"))
                }
                dmsL1 <- rownames(L)
                L0 <- array(0, c(nrowL, nX), list(NULL, names(object$fixDF$X)))
                if (is.null(dmsL2 <- colnames(L))) {
                  L0[, 1:ncolL] <- L
                }
                else {
                  if (any(noMatch <- is.na(match(dmsL2, colnames(L0))))) {
                    stop(paste("Effects", paste(dmsL2[noMatch],
                      collapse = ", "), "not matched"))
                  }
                  L0[, dmsL2] <- L
                }
                L <- L0[noZeroRowL <- as.logical((L0 != 0) %*%
                  rep(1, nX)), , drop = FALSE]
                nrowL <- nrow(L)
                if (is.null(dmsL1)) {
                  dmsL1 <- 1:nrowL
                }
                else {
                  dmsL1 <- dmsL1[noZeroRowL]
                }
                rownames(L) <- dmsL1
                dDF <- unique(object$fixDF$X[noZeroColL <- as.logical(c(rep(1,
                  nrowL) %*% (L != 0)))])
                if (length(dDF) > 1) {
                  ## stop("L may only involve fixed effects with the same denominator DF")
                  warn <- paste( "L involves fixed effects with the different denominator DF:",
                          paste(dDF, collapse=" "), collapse = " ")
                  warning(warn)
                  dDF <- min(dDF)
                }
                lab <- "F-test for linear combination(s)\n"
            }
            nDF <- sum(svd(L)$d > 0)
            c0 <- c(qr.qty(qr(vFix %*% t(L)), c0))[1:nDF]
            Fval <- sum(c0^2)/nDF
            Pval <- 1 - pf(Fval, nDF, dDF)
            aod <- data.frame(nDF, dDF, Fval, Pval)
            names(aod) <- c("numDF", "denDF", "F-value", "p-value")
            attr(aod, "rt") <- rt
            attr(aod, "label") <- lab
            if (!Lmiss) {
                if (nrow(L) > 1)
                  attr(aod, "L") <- L[, noZeroColL, drop = FALSE]
                else attr(aod, "L") <- L[, noZeroColL]
            }
        }
    }
    else {
        ancall <- sys.call()
        ancall$verbose <- ancall$test <- NULL
        object <- list(object, ...)
        termsClass <- unlist(lapply(object, data.class))
        if (!all(match(termsClass, c("gls", "gnls", "lm", "lmList",
            "lme", "nlme", "nlsList", "nls"), 0))) {
            stop(paste("Objects must inherit from classes \"gls\", \"gnls\"",
                "\"lm\",\"lmList\", \"lme\",\"nlme\",\"nlsList\", or \"nls\""))
        }
        resp <- unlist(lapply(object, function(el) deparse(getResponseFormula(el)[[2]])))
        subs <- as.logical(match(resp, resp[1], FALSE))
        if (!all(subs))
            warning(paste("Some fitted objects deleted because",
                "response differs from the first model"))
        if (sum(subs) == 1)
            stop("First model has a different response from the rest")
        object <- object[subs]
        rt <- length(object)
        termsModel <- lapply(object, function(el) formula(el)[-2])
        estMeth <- unlist(lapply(object, function(el) {
            val <- el[["method"]]
            if (is.null(val))
                val <- NA
            val
        }))
        if (length(uEst <- unique(estMeth[!is.na(estMeth)])) >
            1) {
            stop("All fitted objects must have the same estimation method.")
        }
        estMeth[is.na(estMeth)] <- uEst
        REML <- uEst == "REML"
        if (REML) {
            aux <- unlist(lapply(termsModel, function(el) {
                aux <- terms(el)
                val <- paste(sort(attr(aux, "term.labels")),
                  collapse = "&")
                if (attr(aux, "intercept") == 1) {
                  val <- paste(val, "(Intercept)", sep = "&")
                }
                val
            }))
            if (length(unique(aux)) > 1) {
                warning(paste("Fitted objects with different fixed effects.",
                  "REML comparisons are not meaningful."))
            }
        }
        termsCall <- lapply(object, function(el) {
            if (is.null(val <- el$call)) {
                if (is.null(val <- attr(el, "call"))) {
                  stop("Objects must have a \"call\" component or attribute.")
                }
            }
            val
        })
        termsCall <- unlist(lapply(termsCall, function(el) paste(deparse(el),
            collapse = "")))
        aux <- lapply(object, logLik, REML)
        if (length(unique(unlist(lapply(aux, function(el) attr(el,
            "nall"))))) > 1) {
            stop("All fitted objects must use the same number of observations")
        }
        dfModel <- unlist(lapply(aux, function(el) attr(el, "df")))
        logLik <- unlist(lapply(aux, function(el) c(el)))
        AIC <- unlist(lapply(aux, AIC))
        BIC <- unlist(lapply(aux, BIC))
        aod <- data.frame(call = termsCall, Model = (1:rt), df = dfModel,
            AIC = AIC, BIC = BIC, logLik = logLik, check.names = FALSE)
        if (test) {
            ddf <- diff(dfModel)
            if (sum(abs(ddf)) > 0) {
                effects <- rep("", rt)
                for (i in 2:rt) {
                  if (ddf[i - 1] != 0) {
                    effects[i] <- paste(i - 1, i, sep = " vs ")
                  }
                }
                pval <- rep(NA, rt - 1)
                ldf <- as.logical(ddf)
                lratio <- 2 * abs(diff(logLik))
                lratio[!ldf] <- NA
                pval[ldf] <- 1 - pchisq(lratio[ldf], abs(ddf[ldf]))
                aod <- data.frame(aod, Test = effects, L.Ratio = c(NA,
                  lratio), "p-value" = c(NA, pval), check.names = FALSE)
            }
        }
        row.names(aod) <- unlist(lapply(as.list(ancall[-1]),
            deparse))
        attr(aod, "rt") <- rt
        attr(aod, "verbose") <- verbose
    }
    class(aod) <- c("anova.lme", "data.frame")
    aod
}










if( FALSE ) {


#### NOT RUN: ####


glmmPQL <-
function (fixed, random, family, data, correlation, weights,
    control, niter = 10, verbose = TRUE, ...)
{
    if (!require("nlme"))
        stop("package 'nlme' is essential")
    if (is.character(family))
        family <- get(family)
    if (is.function(family))
        family <- family()
    if (is.null(family$family)) {
        print(family)
        stop("'family' not recognized")
    }
    m <- mcall <- Call <- match.call()
    nm <- names(m)[-1]
    keep <- is.element(nm, c("weights", "data", "subset", "na.action"))
    for (i in nm[!keep]) m[[i]] <- NULL
    allvars <- if (is.list(random))
        allvars <- c(all.vars(fixed), names(random), unlist(lapply(random,
            function(x) all.vars(formula(x)))))
    else c(all.vars(fixed), all.vars(random))
    Terms <- if (missing(data))
        terms(fixed)
    else terms(fixed, data = data)
    off <- attr(Terms, "offset")
    if (length(off <- attr(Terms, "offset")))
        allvars <- c(allvars, as.character(attr(Terms, "variables"))[off +
            1])
    m$formula <- as.formula(paste("~", paste(allvars, collapse = "+")))
    environment(m$formula) <- environment(fixed)
    m$drop.unused.levels <- TRUE
    m[[1]] <- as.name("model.frame")
    mf <- eval.parent(m)
    off <- model.offset(mf)
    if (is.null(off))
        off <- 0
    w <- model.weights(mf)
    if (is.null(w))
        w <- rep(1, nrow(mf))
    mf$wts <- w
    fit0 <- glm(formula = fixed, family = family, data = mf,
        weights = wts, ...)
    w <- fit0$prior.weights
    eta <- fit0$linear.predictor
    zz <- eta + fit0$residuals - off
    wz <- fit0$weights
    fam <- family
    nm <- names(mcall)[-1]
    keep <- is.element(nm, c("fixed", "random", "data", "subset",
        "na.action", "control"))
    for (i in nm[!keep]) mcall[[i]] <- NULL
    fixed[[2]] <- quote(zz)
    mcall[["fixed"]] <- fixed
    mcall[[1]] <- as.name("lme")
    mcall$random <- random
    mcall$method <- "ML"
    if (!missing(correlation))
        mcall$correlation <- correlation
    mcall$weights <- quote(varFixed(~invwt))
    mf$zz <- zz
    mf$invwt <- 1/wz
    mcall$data <- mf
    for (i in 1:niter) {
        if (verbose) {
            cat("iteration", i, "\n")
            print(names(mcall))
        }
        fit <- eval(mcall)
        etaold <- eta
        eta <- fitted(fit) + off
        if (sum((eta - etaold)^2) < 1e-06 * sum(eta^2))
            break
        mu <- fam$linkinv(eta)
        mu.eta.val <- fam$mu.eta(eta)
        mf$zz <- eta + (fit0$y - mu)/mu.eta.val - off
        wz <- w * mu.eta.val^2/fam$variance(mu)
        mf$invwt <- 1/wz
        mcall$data <- mf
    }
    attributes(fit$logLik) <- NULL
    fit$call <- Call
    fit$family <- family
    oldClass(fit) <- c("glmmPQL", oldClass(fit))
    fit
}

#log <-function (x, base = exp(1)) {
#    x <- pmax(x,.00000001)
#    if (missing(base)) .Internal(log(x)) else .Internal(log(x, base))
#}




logLik.reStruct  <-
function (object, conLin, ...)
{
    if (any(!is.finite(conLin$Xy)))
        return(-1000000)
    .C("mixed_loglik", as.double(conLin$Xy), as.integer(unlist(conLin$dims)),
        as.double(pdFactor(object)), as.integer(attr(object,
            "settings")), loglik = double(1), double(1), PACKAGE = "nlme")$loglik
}

}


##
##
##  General polynomial splines
##
##
##


##
##
##  fun-new.R   August 8, 2008
##
##




oldtab = function(...) {
    UseMethod("tab")
}

oldtab.default = function( ..., total.margins = TRUE ) {
    # This might -- finally -- be equivalent to table(..., exclude = NULL)
        aa <- list(...)
        if ( length(aa) == 1 && is.list ( aa[[1]]) ) {
             return( do.call("tab", aa[[1]]))
    }
        #disp( names(aa))
        #disp( is.null( names(aa)))
        if ( is.null(names(aa))) {
            nns = names( match.call())
            #disp( nns)
            names(aa) = nns[ 2:(1+length(aa)) ]
        }
        for ( ii in 1:length(aa) ) aa[[ii]] <- factor(aa[[ii]], exclude = NULL)
        #disp( "aa" )
        #disp( aa)
        ret <- do.call("table",aa)
        if ( total.margins) ret = atotal(ret)
        ret
}

oldtab.data.frame = function( dd, fmla , total.margins = TRUE) {
     if ( missing( fmla )) return ( do.call('tab', as.list( dd , total.margins = total.margins)))
     xx = model.frame( fmla , dd , na.action = na.include)
     xx = c(xx, total.margins = total.margins)
     do.call( 'tab', xx)
}

oldtab.formula <- function( fmla, dd, total.margins = TRUE, ...){
               tab( dd, fmla, total.margins = total.margins, ...)
}

if (FALSE) {

    tt = function ( ... ) {
          disp( list(...))
          disp( nargs())
          disp( names(match.call()))
          disp( deparse( match.call()) )
    }

    tt( 1:3, b = 3:4)
    tt( a = 1:3, b = 3:4)


    zdd = data.frame( x = c(1,2,3,2,1,2,NA,1,2,NaN,Inf), A = c(NA,rep( c('a','b'), 5)), B = c(rep( c('x','z'),each = 5),NA))
    tab( zdd, ~ x + A + B)
    with( rima, tab( Category, Year))
    tab( rima[c('Category','Year')])
}


if (FALSE) {
# moved to fun.R and rewritten so it does not need nlme, otherwise
# it can't be used with lme4
up <- function( dd, form , all = F) {


      ## Description:    (G. Monette, June 10, 2006)
      ## Creates a higher level data set by selecting first
      ## after ordering and keeping only invariant variables
      ## To have summary variables they need to have been created
      ## with, e.g., capply
      #sel <- model.frame( form , dd , na.action = na.include )
      #if( !all) vl <- varLevel( dd, form ) else vl <- rep(0,ncol(dd))
      #dd [ na.omit(sapply( split( 1:nrow(dd), sel ), function(x) x[1])), vl < keep + 1, drop = F]
      require(nlme)
      if ( all == F ) gsummary( dd, form = form, invariantsOnly = TRUE )  else gsummary( dd, form = form)
}
}

###
###
###   ith derivative of an expression
###
###
###
###

d <- function( ex, xv,i ) if ( i == 0) ex else D( d(ex, xv, i-1), xv)

# Example
if (FALSE) {
   d( quote( exp( x^2 /2 )), 'x', 3)      # admittedly not efficient code
   d( quote( 3*x^3 + a*x^2), 'x', 2)
}







"
</pre>
"



reorder.factor <- function (x, v, FUN = mean, ...) {

      # Hmisc returns an ordered factor,
      # This returns the same as its input
      if ( inherits(x,'ordered')) ordered(x, levels(x)[order(tapply(v, x, FUN, ...))])
      else factor(x, levels(x)[order(tapply(v, x, FUN, ...))])
}








