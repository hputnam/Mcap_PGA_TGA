\name{gsp}
\Rdversion{1.1}
\alias{gsp}
\alias{sc}
\alias{smsp}
\alias{Xf}
\alias{Xmat}
\alias{qs}
\alias{lsp}
\alias{Cmat}
\concept{spline}
\concept{smoothing spline}
\title{
General splines with variable degrees and smoothnes, smoothing splines
}
\description{
These functions implement a general spline with possibly different degrees in each interval and smoothness
at each knot. The function \code{sc} helps in the construction of linear hypothesis
matrices to estimate and test levels and derivatives of splines at arbitrary points and the saltus
of derivatives that have discontinuities at knots.
}
\usage{
gsp(x, knots, degree = 3, smooth = pmax(pmin(degree[-1], degree[-length(degree)]) - 1, 0), intercept = 0, signif = 3)

lsp(x, knots = quantile(x, pc), exclude = 0, pc = c(0.25, 0.75))

qs(x, knots = quantile(x, pc), exclude = 0, pc = c(0.25, 0.75))

cs(x, knots = quantile(x, pc), exclude = 0, pc = c(0.25, 0.75))

sc(sp, x, D = 0, type = 1)

smsp(x, knots)

The following functions are intended for internal use:

Xf(x, knots, degree = 3, D = 0, right = TRUE, signif = 3)

Xmat(x, degree, D = 0, signif = 3)

Cmat(knots, degree, smooth, intercept = 0, signif = 3)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
value(s) where spline is evaluated
}
  \item{knots}{
vector of knots
}
  \item{degree}{
vector giving the degree of the spline in each interval.
Note the number of intervals is equal to the number of knots + 1. A value of 0
corresponds to a constant in the interval. If the spline should evaluate to 0
in the interval, use the \code{intercept} argument to specify some value in the
interval at which the spline must evaluate to 0.
}
  \item{smooth}{
vector with the degree of smoothness at each knot (0 = continuity, 1 = smooth with continuous first derivative, 2 = continuous second derivative, etc.
}
  \item{intercept}{
value(s) of x at which the spline has value 0, i.e. the value(s) of x for which
yhat is estimated by the intercept term in the model. The default is 0. If NULL,
the spline is not constrained to evaluate to 0 for any x.
}
  \item{signif}{
number of significant digits used to label coefficients
}
  \item{exclude}{number of leading columns to drop from spline matrix: 0: excludes the intercept column,
    1: excludes the linear term as well.
    Terms that are excluded from the spline matrix can be modeled explicitly.
    }
  \item{sp}{
a spline function defined by \code{gsp}. See the examples below.
}
  \item{D}{the degree of a derivative: 0: value of the function, 1: first derivative, 2:
  second derivative, etc.
}
  \item{type}{how a derivative or value of a function is measured at a possible discontinuity
    at a knot: 0: limit from the left, 1: limit from the right, 2: saltus (limit from the right minus
    the limit from the left)
    }
}
\details{
 \code{gsp} generates a matrix of regressors for a spline with knots, degree of polynomials
 in each interval and the degree of smoothness at each knot. Typically, \code{gsp} is used to
 define a function that is then used in a model equation. See the examples below.
 
A function to fit a cubic spline with knots at 5 and 10 is generated with:

\preformatted{
      sp <- function( x ) gsp( x, c(5,10), c(3,3,3), c(2,2))
}

      indicating that a cubic polynomial is used in each of the three
      intervals and that the second derivative is continuous at each knot.

      A 'natural cubic spline' with linear components in each unbounded interval
      would have the form:

\preformatted{
      sp <- function( x ) gsp( x, c(0,5,10,15), c(1,3,3,3,1), c(2,2,2,2))
}

      Quadratic and linear splines, respectively:
\preformatted{
         sp.quad <- function( x ) gsp( x, c(5,10), c(2,2,2), c(1,1))
         sp.lin  <- function( x ) gsp( x, c(5,10), c(1,1,1), c(0,0))
}
      Where the same degree is used for all intervals and knots, it suffices
      to give it once:
\preformatted{
         sp.quad <- function( x ) gsp( x, c(5,10), 2, 1)
         sp.lin  <- function( x ) gsp( x, c(5,10), 1, 0)
}
      An easy way to specify a model in which a knot is dropped is to force
      a degree of continuity equal to the degree of adjoining polynomials, e.g.
      to drop the knot at 10, use:
\preformatted{
         sp.1 <- function( x ) gsp( x, c(5,10), c(3,3,3), c(2,3))
}
      This is sometimes easier than laboriously rewriting the spline function
      for each null hypothesis.

      Depending on the maximal degree of the spline, the range of x
      should not be excessive to avoid numerical problems. The spline
      matrix generated is 'raw' and values of max(abs(x))^max(degree)
      may appear in the matrix.  For example, for a cubic spline, it
      might be desirable to rescale x and/or recenter x so abs(x) < 100
      if that is not already the case. Note that the knots need to be
      correspondingly transformed.

      The naming of coefficients should allow them to be easily interpreted.
      For example:
\preformatted{
        >  zapsmall( gsp ( 0:10, c(3, 7) , c(2,3,2), c(1,1)) )

      D1(0) D2(0) C(3).2   C(3).3 C(7).2
f(0)      0   0.0    0.0  0.00000    0.0
f(1)      1   0.5    0.0  0.00000    0.0
f(2)      2   2.0    0.0  0.00000    0.0
f(3)      3   4.5    0.0  0.00000    0.0
f(4)      4   8.0    0.5  0.16667    0.0
f(5)      5  12.5    2.0  1.33333    0.0
f(6)      6  18.0    4.5  4.50000    0.0
f(7)      7  24.5    8.0 10.66667    0.0
f(8)      8  32.0   12.5 20.66667    0.5
f(9)      9  40.5   18.0 34.66667    2.0
f(10)    10  50.0   24.5 52.66667    4.5
}

      The coefficient for the first regressor is the first derivative
      at x = 0; for the second regressor, the second derivative at 0;
      the third, the saltus (change) in the second derivative at x = 3,
      the fourth, the saltus in the third derivative at x = 3 and, finally,
      the saltus in the second derivative at x = 7.

      Example:
\preformatted{
       >  sp <- function(x) gsp ( x, c(3, 7) , c(2,3,2), c(1,1))
       >  zd <- data.frame( x = seq(0,10, .5), y = seq(0,10,.5)^2 + rnorm( 21))
       >  fit <- lm( y ~ sp( x ), zd)
       >  summary(fit)
       >         Ls <-cbind( 0, sc( sp, c(1,2,3,3,3,5,7,7,7,8), D=3,
       +                        type = c(0,0,0,1,2,0,0,1,2,0)))
       >  zapsmall( Ls )

                D1(0) D2(0) C(3).2 C(3).3 C(7).2
D3(1)         0     0     0      0      0      0
D3(2)         0     0     0      0      0      0
D3(3-)        0     0     0      0      0      0
D3(3+)        0     0     0      0      1      0
D3(3+)-D3(3-) 0     0     0      0      1      0
D3(5)         0     0     0      0      1      0
D3(7-)        0     0     0      0      1      0
D3(7+)        0     0     0      0      0      0
D3(7+)-D3(7-) 0     0     0      0     -1      0
D3(8)         0     0     0      0      0      0

        >  wald( fit, list( 'third derivatives' = Ls))

                  numDF denDF  F.value p.value
third derivatives     1    15 2.013582 0.17634

                 Estimate Std.Error  DF   t-value p-value Lower 0.95 Upper 0.95
  D3(1)          0.000000  0.000000  15 -1.123777 0.27877   0.000000   0.000000
  D3(2)          0.000000  0.000000  15 -1.123777 0.27877   0.000000   0.000000
  D3(3-)         0.000000  0.000000  15 -1.123777 0.27877   0.000000   0.000000
  D3(3+)         0.927625  0.653714  15  1.419008 0.17634  -0.465734   2.320984
  D3(3+)-D3(3-)  0.927625  0.653714  15  1.419008 0.17634  -0.465734   2.320984
  D3(5)          0.927625  0.653714  15  1.419008 0.17634  -0.465734   2.320984
  D3(7-)         0.927625  0.653714  15  1.419008 0.17634  -0.465734   2.320984
  D3(7+)         0.000000  0.000000 Inf       NaN     NaN   0.000000   0.000000
  D3(7+)-D3(7-) -0.927625  0.653714  15 -1.419008 0.17634  -2.320984   0.465734
  D3(8)          0.000000  0.000000 Inf       NaN     NaN   0.000000   0.000000

Warning messages:
1: In min(dfs[x != 0]) : no non-missing arguments to min; returning Inf
2: In min(dfs[x != 0]) : no non-missing arguments to min; returning Inf

}
          Note that some coefficients that are 0 by design may lead to invalid
          DRs and t-values.


      \code{sc} generates a portion of a hypothesis matrix for the coefficients
     of a general spline constructed with \code{gsp}
With:
\preformatted{
   sc( sp, x, D, type ):

     sp    is the spline function for which coefficients are required
     x     values at which spline is evaluated
     D     order of derivative: 0 = value of spline, 1 = first derivative, etc.
     type  at knots: 0 limit from the left, 1 from the right, 2 is saltus (i.e. jump
           from left to right)
}

Warning: \code{sc} will not work correctly if the function defining the spline
              transforms the variable, e.g. sp <- function(x) gsp( x/100, knot=2 )

Example:

\preformatted{
     simd <- data.frame( age = rep(1:50, 2), y = sin(2*pi*(1:100)/5)+ rnorm(100),
             G = rep( c("male","female"), c(50,50)))

     sp <- function(x) gsp( x, knots = c(10,25,40), degree = c(1,2,2,1),
                       smooth = c(1,1,1))

     fit <- lm( y ~ sp(age)*G, simd)
     xyplot( predict(fit) ~ age , simd, groups = G,type = "l")
     summary(fit)  # convenient display

# output:

Call:
lm(formula = y ~ sp(age) * G, data = simd)

Residuals:
    Min      1Q  Median      3Q     Max
-2.5249 -0.7765 -0.0760  0.7882  2.6265

Coefficients:
                      Estimate Std. Error t value Pr(>|t|)
(Intercept)           0.733267   0.605086   1.212    0.229
sp(age)D1(0)         -0.084219   0.055163  -1.527    0.130
sp(age)C(10).2        0.010984   0.006910   1.590    0.115
sp(age)C(25).2       -0.023034   0.012881  -1.788    0.077 .
Gmale                -0.307665   0.855721  -0.360    0.720
sp(age)D1(0):Gmale    0.058384   0.078012   0.748    0.456
sp(age)C(10).2:Gmale -0.010556   0.009773  -1.080    0.283
sp(age)C(25).2:Gmale  0.026410   0.018216   1.450    0.150
---


Residual standard error: 1.224 on 92 degrees of freedom
Multiple R-squared: 0.0814,     Adjusted R-squared: 0.0115
F-statistic: 1.165 on 7 and 92 DF,  p-value: 0.3308
# end of output

     L0 <- list( "hat" =
        rbind( "females at age=20" = c( 1, sc(sp,20), 0, 0* sc(sp,20)),
               "males at age=20" = c( 1, sc(sp,20), 1, 1* sc(sp,20))),
        "male-female" = rbind( "at 20" = c( 0 , 0*sc(sp,20), 1, 1*sc(sp,20))))
     wald( fit, L0 )

     ...

     L1 <- list("D(yhat)/D(age)"=
        rbind( "female at age = 25" = c(0, sc(sp,25,1), 0, 0*sc(sp,25,1)),
               "male at x = 25" = c(0, sc(sp,25,1), 0, 1*sc(sp,25,1))))
     wald( fit, L1)
# output:
               numDF denDF  F.value p.value
D(yhat)/D(age)     2    92 1.057307 0.35157

                      Estimate Std.Error DF   t-value p-value Lower 0.95 Upper 0.95
female at age = 25  0.080544  0.056974 92  1.413694 0.16083  -0.032612   0.193700
male at x = 25     -0.019412  0.056974 92 -0.340712 0.73410  -0.132568   0.093744
}

Overview of utility functions:
\preformatted{
Xmat = function( x, degree, D = 0, signif = 3)
     design/estimation matrix for f[D](x) where f(x) is polynomial of degree degree.

Xf =  function(   x, knots, degree = 3, D = 0, right = TRUE , signif = 3)
     uses Xmat to form 'full' matrix with blocks determined by knots intervals

Cmat = function( knots, degree, smooth, intercept = 0, signif = 3)
     linear constraints

Emat = function(  knots, degree, smooth , intercept = FALSE, signif = 3)
     estimates - not necessarily a basis

basis  = function( X , coef = FALSE )
     selects linear independent subset of columns of X

spline.T = function( knots, degree, smooth, intercept = 0, signif = 3 )
     full transformation of Xf to spline basis and constraints

spline.E = function( knots, degree, smooth, intercept = 0, signif = 3 )
     transformation for spline basis (first r columns of spline.T)

gsp = function( x , knots, degree = 3 , smooth = pmax(pmin( degree[-1],
        degree[ - length(degree)]) - 1,0 ), intercept = 0, signif = 3)

 
}
}
\section{Warning}{
The variables generated by \code{gsp} are designed so the coefficients are interpretable as changes
in derivatives at knots. The resulting matrix is not designed to have optimal numerical properties.

The intermediate matrices generated by \code{gsp} will contain \code{x} raised to the power equal to
the highest degree in \code{degree}. The values of \code{x} should be scaled to avoid
excessively large values in the spline matrix and ensuing numerical problems.
}
\value{
\code{gsp} returns a matrix generating a spline. \code{cs}, \code{qs} and \code{lsp} return matrices
generating cubic,
quadratic and linear splines respectively.

\code{smsp}, whose code is adapted from function in the package \code{lmeSplines},
generates a smoothing spline to be used in the random effects portion of
a call to \code{lme}.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Monette, G. \email{georges@yorku.ca}
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{wald}}
}
\examples{
## Fitting a quadratic spline
simd <- data.frame( age = rep(1:50, 2), y = sin(2*pi*(1:100)/5) + rnorm(100),
          G = rep( c('male','female'), c(50,50)))
# define a function generating the spline
sp <- function(x) gsp( x, knots = c(10,25,40), degree = c(1,2,2,1),
          smooth = c(1,1,1))

fit <- lm( y ~ sp(age)*G, simd)

require(lattice)
xyplot( predict(fit) ~ age , simd, groups = G,type = 'l')
summary(fit)

## Linear hypotheses
L <- list( "Overall test of slopes at 20" = rbind(
      "Female slope at age 20" =  c( F20 <- cbind( 0 , sc( sp, 20, D = 1), 0 , 0 * sc( sp, 20, D = 1))),
      "Male slope at age 20" =  c( M20 <- cbind( 0 , sc( sp, 20, D = 1), 0 , 1 * sc( sp, 20, D = 1))),
      "Difference" = c(M20 - F20))
      )
wald( fit, L)

## Right and left second derivatives at knots and saltus
L <- list( "Second derivatives and saltus for female curve at knot at 25" =
          cbind( 0, sc( sp, c(25,25,25), D = 2, type =c(0,1,2)), 0,0,0,0))
L
wald( fit, L )

## Smoothing splines
      library(nlme)
      data(Spruce)
      Spruce$all <- 1
      range( Spruce$days)
      sp <- function(x) smsp ( x, seq( 150, 675, 5))
      spruce.fit1 <- lme(logSize ~ days, data=Spruce,
          random=list(all= pdIdent(~sp(days) -1),
          plot=~1, Tree=~1))
      summary(spruce.fit1)
      pred <- expand.grid( days = seq( 160, 670, 10), all = 1)
      pred$logSize <- predict( spruce.fit1, newdata = pred, level = 1)
      require( lattice )
      xyplot( logSize ~ days, pred, type = 'l')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
