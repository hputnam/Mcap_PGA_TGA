\name{anova.lme}
\Rdversion{1.1}
\alias{anova.lme}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
anova for lme
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
anova.lme(object, ..., test = TRUE, type = c("sequential", "marginal"), adjustSigma = TRUE, Terms, L, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
  \item{test}{
%%     ~~Describe \code{test} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{adjustSigma}{
%%     ~~Describe \code{adjustSigma} here~~
}
  \item{Terms}{
%%     ~~Describe \code{Terms} here~~
}
  \item{L}{
%%     ~~Describe \code{L} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
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
function (object, ..., test = TRUE, type = c("sequential", "marginal"), 
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
                L <- L0[noZeroRowL <- as.logical((L0 != 0) \%*\% 
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
                  nrowL) \%*\% (L != 0)))])
                if (length(dDF) > 1) {
                  warn <- paste("L involves fixed effects with the different denominator DF:", 
                    paste(dDF, collapse = " "), collapse = " ")
                  warning(warn)
                  dDF <- min(dDF)
                }
                lab <- "F-test for linear combination(s)\n"
            }
            nDF <- sum(svd(L)$d > 0)
            c0 <- c(qr.qty(qr(vFix \%*\% t(L)), c0))[1:nDF]
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
                  lratio), `p-value` = c(NA, pval), check.names = FALSE)
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
