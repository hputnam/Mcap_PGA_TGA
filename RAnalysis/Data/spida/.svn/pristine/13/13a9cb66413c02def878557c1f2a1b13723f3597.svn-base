
## wald.R
## This a collection of functions designed to facilitate testing hypotheses
## with Wald tests.
## The methods are readily extended to any fitting method that returns a vector
## of estimated parameters and an estimated variance matrix.
## Extensions are implemented through the 'getFix' generic function.


glh <- function( ...) wald( ...)    # previous name for 'wald' function



wald <- function(fit, Llist = "",clevel=0.95, data = NULL, debug = FALSE , maxrows = 25,
     full = FALSE, fixed = FALSE, invert = FALSE, method = 'svd', help = FALSE) {
if(help) {
cat("
wald:  General Linear Hypothesis with Wald test
      for lm, glm, lme, nlme and lmer objects.

      Can be extended to other objects (e.g.) 'glm' by writing 'getFix.glm'

      Usage:
         wald ( fit, L ) where L is a hypothesis matrix
         wald ( fit, 'pat' ) where 'pat' a regular expression (see ?regex)
           used to match names of coefficients of fixed effects.
           e.g. wald( fit, ':.*:') tests all 2nd and higher order interactions.
         wald ( fit, c(2,5,6)) to test 2nd, 5th and 6th coefficients.
         wald ( fit, list( hyp1= c(2,5,6), H2 = 'pat')) for more than one hypothesis
           matrix

      There are number of functions to help construct hypothesis matrices:

         Lform ( fit, list( arg1, arg2, arg3), data = dframe)
           creates an L matrix by evaluating arg1, arg2, arg3 in the dframe
           environment to generate columns of the L matrix. 'dframe' is
           model.frame(fit) by default. For example, the derivative with respect
           to ses could be evaluated at each point in the following:

               hs <- read.csv('http://www.math.yorku.ca/~georges/Data/hs.csv')
               fit <- lme( mathach ~ (ses + I(ses^2)) * Sex, hs, random = ~ 1 + ses| school)
               wald( fit, 'Sex')  # sig. overall effect of Sex
               wald( fit, ':Sex') # but no evidence of interaction with ses
               wald( fit, '\\^2') # nor of curvature

               # but we continue for the sake of illustration

               L <- Lform( fit, list( 0, 1, 2*ses, 0, Sex == 'Male', (Sex == 'Male')*2*ses), hs)

               ww <- wald ( fit, L )
               dd <- as.data.frame(ww, se = 2)
               head( dd)
               library(lattice)
               xyplot( coef + coefp + coefm ~ ses | Sex, dd, main=
                       'Increase in predicted mathach per unit increase in ses')

" )
 return(invisible(NULL))
 }
 
        if (full ) return( wald ( fit, model.matrix(fit)))
          dataf <- function(x,...) {
           x <- cbind(x)
           rn <- rownames(x)
           if( length( unique(rn)) < length(rn)) rownames(x) <- NULL
           data.frame(x,...)
         }
     as.dataf <- function(x,...) {
           x <- cbind(x)
           rn <- rownames(x)
           if( length( unique(rn)) < length(rn)) rownames(x) <- NULL
           as.data.frame(x,...)
         }

           unique.rownames <- function(x) {
                ret <- c(tapply(1:length(x), x , function(xx) {
                    if ( length(xx) == 1) ""
                    else 1:length(xx)
                })) [tapply(1:length(x),x)]
                ret <- paste(x,ret,sep="")
                ret
           }
           if(debug) disp( Llist)
           if( is.character(Llist) ) Llist <- structure(list(Llist),names=Llist)
           if(!is.list(Llist)) Llist <- list(Llist)

           ret <- list()
           fix <- getFix(fit)
           if(debug) disp(fix)
           beta <- fix$fixed
           vc <- fix$vcov

           dfs <- fix$df
           if(debug) disp(Llist)
           for (ii in 1:length(Llist)) {
               ret[[ii]] <- list()
               Larg <- Llist[[ii]]
               if(debug) {
                    disp(ii)
                    disp(Larg)
               }
               L <- NULL
               if ( is.character(Larg)) {
                  L <- Lmat(fit,Larg, fixed = fixed, invert = invert)
               } else {
                  if ( is.numeric(Larg)) {
                     if ( is.null(dim(Larg))) {
                        if(debug) disp(dim(Larg))
                        if ( (length(Larg) < length(beta)) && (all(Larg>0)||all(Larg<0)) ) {
                              L <- diag(length(beta))[Larg,]
                              dimnames(L) <- list( names(beta)[Larg], names(beta))
                        } else L <- rbind( Larg )
                     }
                     else L <- Larg
                  }
               }
               if (debug) {
                  disp(Larg)
                  disp(L)
               }
               ## Delete coefficients that are NA
               Ldata <- attr( L , 'data')
               L <- L[, !is.na(beta),drop = FALSE]
               attr(L,'data') <- Ldata
               beta <- beta[ !is.na(beta) ]
               
               ## Anova
               if( method == 'qr' ) {
                   qqr <- qr(na.omit(t(L)))
                   #Qqr <- Q(t(L))
                   L.rank <- qqr$rank
                   #L.rank <- attr(Qqr,'rank')
                   #L.miss <- attr(Qqr,'miss')
                   if(debug)disp( t( qr.Q(qqr)))
                   L.full <- t(qr.Q(qqr))[ 1:L.rank,,drop=FALSE]
                   #L.full <- t(Qqr[!L.miss,])[ 1:L.rank,,drop=F]
               } else if ( method == 'svd' ) {
                   if(debug)disp( t(na.omit(t(L))))
                   sv <- svd( t(na.omit(t(L))) , nu = 0 )
                   if(debug)disp( sv )
                   tol.fac <- max( dim(L) ) * max( sv$d )
                   if(debug)disp( tol.fac )
                   if ( tol.fac > 1e6 ) warning( "Poorly conditioned L matrix, calculated numDF may be incorrect")
                   tol <- tol.fac * .Machine$double.eps
                   if(debug)disp( tol )
                   L.rank <- sum( sv$d > tol )
                   if(debug)disp( L.rank )
                   if(debug)disp( t(sv$v))
                   L.full <- t(sv$v)[seq_len(L.rank),,drop = FALSE]
               } else stop("method not implemented: choose 'svd' or 'qr'")
               
               # from package(corpcor)
               # Note that the definition tol= max(dim(m))*max(D)*.Machine$double.eps
               # is exactly compatible with the conventions used in "Octave" or "Matlab".


               if (debug && method == "qr") {
                  disp(qqr)
                  disp(dim(L.full))
                  disp(dim(vc))
                  disp(vc)
               }

               if (debug) disp(L.full)
               if (debug) disp( vc )
               
               vv <-  L.full %*% vc %*% t(L.full)
               eta.hat <- L.full %*% beta
               Fstat <- (t(eta.hat) %*% qr.solve(vv,eta.hat,tol=1e-10)) / L.rank
               included.effects <- apply(L,2,function(x) sum(abs(x),na.rm=TRUE)) != 0
               denDF <- min( dfs[included.effects])
               numDF <- L.rank
               ret[[ii]]$anova <- list(numDF = numDF, denDF = denDF,
                               "F-value" = Fstat,
                               "p-value" = pf(Fstat, numDF, denDF, lower.tail = FALSE))
               ## Estimate
               etahat <- L %*% beta
               if( nrow(L) <= maxrows ) {
                    etavar <- L %*% vc %*% t(L)
                    etasd <- sqrt( diag( etavar ))
               } else {
                    etavar <- NULL
                    etasd <- sqrt( apply( L * (L%*%vc), 1, sum))
               }

               denDF <- apply( L , 1 , function(x,dfs) min( dfs[x!=0]), dfs = dfs)
               
               aod <- cbind( Estimate=c(etahat),
                   Std.Error = etasd,
                   DF = denDF,
                   "t-value" = c(etahat/etasd),
                   "p-value" = 2*pt(abs(etahat/etasd), denDF, lower.tail =FALSE))
                colnames(aod)[ncol(aod)] <- 'p-value'
             if (debug ) disp(aod)
             if ( !is.null(clevel) ) {
                #print(aod)
                #print(aod[,'DF'])
                #print(aod[,'etasd'])
                 hw <- qt(1 - (1-clevel)/2, aod[,'DF']) * aod[,'Std.Error']
                 #print(hw)
                 aod <- cbind( aod, LL = aod[,"Estimate"] - hw, UL = aod[,"Estimate"] + hw)
                 #print(aod)
                 if (debug ) disp(colnames(aod))
                 labs <- paste(c("Lower","Upper"), format(clevel))
                 colnames(aod) [ ncol(aod) + c(-1,0)] <- labs
             }
             if (debug ) disp(rownames(aod))
             aod <- as.dataf(aod)

             rownames(aod) <- rownames(as.dataf(L))
             labs(aod) <- names(dimnames(L))[1]
             ret[[ii]]$estimate <- aod
             ret[[ii]]$coef <- c(etahat)
             ret[[ii]]$vcov <- etavar
             ret[[ii]]$L <- L
             ret[[ii]]$se <- etasd
             ret[[ii]]$L.full <- L.full
             ret[[ii]]$L.rank <- L.rank
             if( debug ) disp(attr(L,'data'))
             ret[[ii]]$data <- attr(L,'data')
        }
        names(ret) <- names(Llist)
        attr(ret,"class") <- "wald"
        ret
}

print.wald <- function(x,round = 6, pround = 5,...) {
    pformat <- function(x, digits = pround) {
        x <- format(xx <- round(x,digits))
        x[ as.double(xx) == 0 ] <- paste(c("<.",rep('0',digits-1),'1'),collapse="")
        x
    }
    rnd <- function(x,digits) {
        if (is.numeric(x)) x <- round(x,digits=digits)
        format(x)
    }
             for( ii in 1:length(x)) {
                  nn <- names(x)[ii]
                  tt <- x[[ii]]
                  ta <- tt$anova

                  ta[["p-value"]] <- pformat(ta[["p-value"]])
                  print(as.data.frame(ta,row.names=nn))
                  te <- tt$estimate
                  rowlab <- attr(te,"labs")

                  te[,'p-value'] <- pformat( te[,'p-value'])
                  if ( !is.null(round)) {
                     for ( ii in 1:length(te)) {
                         te[[ii]] <- rnd(te[[ii]],digits=round)
                     }
                  }
                  labs(te) <- rowlab
                  print(te,digits=round,...)
                  cat("\n")
             }
        invisible(x)
}

as.data.frame.wald <- function( obj, se , digits = 3 , sep = "" , which = 1 ) {
# modified by GM 2010_09_20 to avoid problems with coefs with duplicate rownames
         dataf <- function(x,...) {
           x <- cbind(x)
           rn <- rownames(x)
           if( length( unique(rn)) < length(rn)) rownames(x) <- NULL
           data.frame(x,...)
         }
         obj = obj [which]
         if ( length(obj) == 1) {

            cf <- obj[[1]]$coef
            ret <- if ( missing(se)) data.frame( coef = cf, se = obj[[1]]$se)
            else {
               if ( is.null( names(se))) names(se) <-
                             sapply(se,function(x) as.character(round(x,digits)))
               SE <- obj[[1]]$se
               SEmat <- cbind(SE) %*% rbind(se)
               cplus <- cf + SEmat
               cminus <- cf - SEmat
               colnames(cplus) <- paste( "U",colnames(cplus),sep=sep)
               colnames(cminus) <- paste( "L",colnames(cminus),sep=sep)
               ret <- data.frame( coef = cf)
               ret <- cbind( ret, cplus, cminus)
            }
            if( is.null(dd <- obj[[1]]$data)) return( ret)
            else return( cbind(ret, dd))
         }
         else ret <- lapply( obj, as.data.frame.wald)
         ret
}

coef.wald <- function( obj , se = FALSE ) {
         if ( length(obj) == 1) {
            ret <-
            ret <- obj[[1]]$coef
            if ( is.logical(se) && (se == TRUE) ) {
               ret <- cbind( coef = ret, se = obj[[1]]$se)

            } else if ( se > 0 ){
               ret <- cbind( coef = ret, coefp = ret+se*obj[[1]]$se,
                      coefm = ret - se*obj[[1]]$se)
               attr(ret,'factor') <- se
            }
         }
         else ret <- sapply( obj, coef.wald )
         ret
}


###  getFix: function designed to be used internally to get coef, var(coef) and df.resid

getFix.rdc <- function(fit, ...) UseMethod("getFix")

getFix.default <- function(fit, ...) stop("not yet implemented")

getFix.rdc.lmer <- function(fit, ...) {
            ret <- list()
            ret$fixed <- fixef(fit)
            ret$vcov <- vcov(fit)
            ret$df <- Matrix:::getFixDF(fit)
            ret
}

getFix.rdc.lm <- function(fit, ...) {
          ret <- list()
          ret$fixed <- coef(fit)
          ret$vcov <- vcov(fit)
          ret$df <- fit$df.residuals
          ret
}


Lform <- function( fit, form, data = model.frame(fit)) {
help <- "
      Creates an L matrix using formulas evaluated in 'data' for
      each column of the L matrix
      Example:

      library(car)
      fit <- lm( income ~ (education + I(education^2) )* type, Prestige)
      summary(fit)

      . . .
      Coefficients:
                               Estimate Std. Error t value Pr(>|t|)
       (Intercept)                891.3    23889.1   0.037  0.97032
       education                  210.0     5638.8   0.037  0.97037
       I(education^2)              38.3      328.3   0.117  0.90740
       typeprof                191523.2    63022.0   3.039  0.00312 **
       typewc                   25692.3    73888.0   0.348  0.72887
       education:typeprof      -28133.0    10236.0  -2.748  0.00725 **
       education:typewc         -4485.4    14007.9  -0.320  0.74956
       I(education^2):typeprof   1017.5      451.8   2.252  0.02679 *
       I(education^2):typewc      170.9      671.6   0.255  0.79967
      . . .

      # estimate the marginal value of occupation for each occupation in the data set

      L <- list( 'marginal value of education' =Lform( fit,
             form = list( 0,1,2*education, 0,0, type == 'prof', type == 'wc',
                2*education *(type=='prof'), 2* education * (type == 'wc')),
             data = Prestige))
      wald( fit, L )
      chat <- coef( wald( fit, L ), se = 2)
      xyplot( coef +coefp+coefm ~ education | type, cbind(Prestige,chat)[order(Prestige$education),],
              type = 'l')
      xyplot( chat~ education | type, Prestige)

"
      gg <- getFix(fit)
      #print(gg)
      Lsub <- do.call(cbind,eval( substitute( form ), data))
      if( (nrow(Lsub) != nrow( data))) {
          if ((nrow(Lsub)==1)) Lsub <- Lsub[rep(1,nrow(data)),]
          else stop('nrow(Lsub) != nrow(data)')
      }
      if( is.null( colnames(Lsub))) colnames(Lsub) <- rep('',ncol(Lsub))
      # disp(Lsub)
      L <- matrix( 0, nrow = nrow(Lsub), ncol = length( gg$fixed))
      rownames(L) <- rownames(data)
      colnames(L) <- names( gg$fixed)
      Lpos <- Lsub[, colnames(Lsub) == '', drop = FALSE]
      # disp(Lpos)
      Lnamed <- Lsub[ , colnames(Lsub) !='', drop  = FALSE]
      # disp(Lnamed)
      for ( ip in seq_len( ncol( Lpos ))) L[,ip] <- Lpos[,ip]
      if ( ncol( Lnamed ) > 0 ) {
         if ( length( unknown <- setdiff( colnames(Lnamed) , colnames(L)))) {
              stop( paste("Unknown effect(s):" , unknown, collapse = " "))
         }
         for ( nn in colnames(Lnamed)) L[,nn] <- Lnamed[,nn]
      }
      attr(L,"data") <- data
      L
}

"
 L <- Lform( fit,
             form = list( 0,1,2*education, 0,0, type == 'prof', type == 'wc',
                2*education *(type=='prof'), 2* education * (type == 'wc')),
             data = Prestige)

Lform( fit, list(age, Sex))
model.frame(fit )
do.call( 'cbind',
Lform(fit, list( 0, 1, a, a*c, a*(f=='a')), list( a = 1:4, c = 11:14, f = factor(c('a','a','b','c'))))
)

Lform(fit, list( 'a:b'= 0, b=1, a, a*c, a*(f=='a')), list( a = 1:4, c = 11:14, f = factor(c('a','a','b','c'))))
"




Lmat <- function(fit, pattern, fixed = FALSE, invert = FALSE, debug = FALSE) {
     # pattern can be a character used as a regular expression in grep
     # or a list with each component generating  a row of the matrix
     umatch <- function( pat, x ) {
            ret <- rep(0,length(pat))
            for ( ii in 1:length(pat)) {
                imatch <- grep(pat[ii], x, fixed= fixed, invert = invert)
                if ( length(imatch) != 1) {
                   cat("\nBad match of:\n")
                   print(pat)
                   cat("in:\n")
                   print(x)
                   stop("Bad match")
                }
                ret[ii] <- imatch
            }
            ret
     }
     if ( is.character(fit)) {
        x <- pattern
        pattern <- fit
        fit <- x
     }
     fe <- getFix(fit)$fixed
     ne <- names(fe)
     if (is.character(pattern)) {
        L.indices <- grep(pattern,names(fe))
        ret <- diag( length(fe)) [L.indices,,drop = FALSE]
        if (debug) disp(ret)
        rownames(ret) <- names(fe) [L.indices]
        labs(ret) <- "Coefficients"
     } else if (is.list(pattern)){
        ret <- matrix(0, nrow = length(pattern), ncol = length(fe))
        colnames(ret) <- ne
        for ( ii in 1:length(pattern)) {
            Lcoefs <- pattern[[ii]]
            pos <- umatch(names(Lcoefs), ne)
            if ( any( is.na(pos))) stop("Names of L coefs not matched in fixed effects")
            ret[ii, pos] <- Lcoefs
        }
        rownames(ret) <- names(pattern)
      }
      labs(ret) <- "Coefficients"
      ret
}

Ldiff.old <- function(fit, pat, levnames = c(reflevel,substring(rownames(L),cut+1)),
                       reflevel = "<ref>", cut = nchar(pat)) {
         L <- Lmat(fit, pat)
         nam <- rownames(L)
         n <- nrow(L)
         if(n < 2) return(L)
         plus <- unlist( apply( rbind( 2:n), 2, seq, n))
         minus <- rep(1:(n-1), (n-1):1)
         Lp <- L[ plus, ]
         Lm <- L[ minus, ]
         Lret <- rbind( L, Lp - Lm)
         rn <- paste( levnames[ c(1:n,plus) + 1], levnames[ c(rep(0,n), minus)+1], sep = " - ")
         rownames(Lret) <- rn
         Lret
}

Ldiff.rdc <- function( fit, nam , ref = "no longer used") {
      # based on Lmu
      # Tests all pairwise difference in factor with model with Intercept term
      Lm <- Lmu(fit, nam)
      levs <- rownames(Lm)
      n <- nrow(Lm)
      if (n < 2) return (Lm)
      plus <- unlist( apply ( rbind(2:n), 2, seq, n))
      minus <- rep(1:(n-1), (n-1):1)
      Lret <- Lm[plus,] - Lm[minus,]
      rn <- paste( levs [plus], levs[minus] , sep = " - ")
      rownames(Lret) <- rn
      Lret
}



Ldiff <- function( fit, pat, levnames = c(reflevel,substring(rownames(L),cut+1)),
         reflevel ="<ref>", cut=nchar(pat),verbose=F) {
      L <- Lmat(fit, paste("^",pat,sep=""))
      nam <- rownames(L)
      n <- nrow(L)
      zm <- matrix(1:n,nrow=n,ncol=n)
      plus <- zm[ col(zm) < row(zm)]
      minus <- rep(1:(n-1), (n-1):1)
      Lp <- L[plus,]
      Lm <- L[minus,]
      Lret <- rbind( L, Lp - Lm)
         pnames <- levnames [ c(1:n, plus) +1]
      mnames <- levnames [ c(rep(0,n), minus) + 1]
      if (verbose) {
        print(levnames)
        print(plus)
        print(minus)
        print(Lp)
        print(Lm)
        print(L)
        print(Lret)
        print(pnames)
        print(mnames)
      }
      rn <- paste( levnames[ c(1:n,plus)+1], levnames[ c(rep(0,n),minus) + 1], sep = " - ")
      rownames(Lret) <- rn
      Lret
}
 


Lmu <- function(fit, nam, verbose = 0) {
       ## "Works only if 'nam' is a factor and a main effect and model has Intercept")
       if ( class(fit) != 'lmer' ) stop( "only implemented for lmer")
       v <- fit@frame[[nam]]
       if( !is.factor(v)) stop ("nam needs to specify the name of a factor")
       levs <- levels(v)
       if( verbose > 0) print(levs)
       cmat <- contrasts(v)
       if( verbose > 0) print(cmat)
       #  print(cmat)
       fe <- getFix(fit)$fixed
       if( verbose > 0) print(fe)
       if ( substring(nam,1,1) != '^') nam <- paste("^",nam,sep="")
       L.indices <- grep(nam,names(fe))
       if( verbose > 0) print(L.indices)
       L <- matrix(0,nrow=length(levs), ncol = length(fe))

       colnames(L) <- names(fe)
       if( verbose > 0) print(L)
          rownames(L) <- levs
       L[,L.indices] <- cmat
       if('(Intercept)' %in% colnames(L)) L[,'(Intercept)'] <- 1
       L
}



Lc <- function(fit, nam, ref = 1, verbose = 0) {
       ## Comparisons with one level
       ## Use Lmu
       ## "Works only if 'nam' is a factor and a main effect and model has Intercept?")
       if ( class(fit) != 'lmer' ) stop( "only implemented for lmer")
       L <- Lmu( fit, nam)
       Lref <- L[ ref,,drop = FALSE]
       index <- 1:nrow(L)
       names(index) <- rownames(L)
       refind <- index[ref]
       if (length(refind) != 1) stop( paste( ref, "does not refer to a single level"))
       Lret <- L[-refind,]
       Lret <- Lret - cbind( rep(1,nrow(Lret))) %*% Lref
       attr(Lret,"heading") <- paste("Comparisons with reference level:", rownames(L)[refind])
       Lret
}

Lrm <- function(fit, nam, vals = 1:nrow(L.mu)) {
    ## Repeated measures polynomial contrasts
    ## Uses Lmu
    ## "Works only if 'nam' is a factor and a main effect and model has Intercept?")
    ##
    L.mu <- Lmu(fit, nam)
    # print(L.mu)
    pp <- cbind( 1, Poly(vals, nrow(L.mu) -1))
    # print(pp)
    ortho <- Q(pp)[,-1] # (linear, quad, etc.)
    # print(ortho)
    ortho <- ortho[,-1]
    maxp <- max( 5, nrow(L.mu))
    colnames(ortho) <- c('linear','quadratic','cubic', paste("^",4:maxp,sep=''))[1:ncol(ortho)]
    L <- t(ortho) %*% L.mu
    L
}
# Lrm(fit, "SRH94")




Lequal <- function(fit, pat) {
       # Test for equality within levels of pat using all differences
         L <- Lmat(fit, pat)
         nam <- rownames(L)
         n <- nrow(L)
         if(n < 2) return(L)
         plus <- unlist( apply( rbind( 2:n), 2, seq, n))
         minus <- rep(1:(n-1), (n-1):1)
         Lp <- L[ plus, ]
         Lm <- L[ minus, ]
         Lret <- rbind( Lp - Lm)
         rn <- paste( nam[plus], nam[minus], sep = " - ")
         rownames(Lret) <- rn
         Lret
}




#Lc <- function(fit, vec ){
#   fe <- getFix(fit)$fixed
#   ret <- 0 * fe
#   if ( is.null(names(vec))) ret[] <-
#}
# Lmu(fit,"SRH")

Lall <- function( fit , nam ) {
        if ( class(fit) != 'lmer' ) stop( "only implemented for lmer")
        v <- fit@frame[[nam]]
        if( !is.factor(v)) stop ("nam needs to specify the name of a factor")
        lev0 <- levels(v)[1]
        ret <-list()
        namf <- nam
        if ( substring(namf,1,1) != "^") namf <- paste("^", namf, sep ="")
        ret[[ nam ]] <- Lmat( fit, namf)
        ret[[ paste(nam,"mu",sep = '.') ]] <- Lmu(fit,nam)
        ret[[ paste(nam,"diff",sep = '.') ]] <- Ldiff( fit , nam)
        ret

}

