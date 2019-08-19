\name{wald}
\Rdversion{1.1}
\alias{wald}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
General Linear Hypothesis for the 'fixed' portion of a model with Wald test
}
\description{
Tests a general linear hypothesis for the linear fixed portion of a model. The hypothesis
can be specified in a variety of ways such as a hypothesis matrix or a pattern that is used as
a regular expression to be matched with the names of coefficients of the model.
A number of tools are available to facilitate the generation of hypothesis matrices.
}
\usage{
wald(fit, Llist = "", clevel = 0.95, data = NULL, debug = FALSE, maxrows = 25,
full = FALSE, fixed = FALSE, invert = FALSE, method = "svd", help = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{ a model for which a \code{getFix} method exists.
%%     ~~Describe \code{fit} here~~
}
  \item{Llist}{ a hypothesis matrix or a pattern to be matched or a list of these
%%     ~~Describe \code{Llist} here~~
}
  \item{clevel}{ level for confidence intervals
%%     ~~Describe \code{clevel} here~~
}
  \item{data}{ used for 'data' attribute of value returned
%%     ~~Describe \code{data} here~~
}
  \item{debug}{
%%     ~~Describe \code{debug} here~~
}
  \item{maxrows}{ maximum number of rows of hypothesis matrix for which a full variance-covariance matrix is returned
%%     ~~Describe \code{maxrows} here~~
}
  \item{full}{ if TRUE, the hypothesis matrix is the model matrix for \code{fit} such that
  the estimated coefficients are the predicted values for the fixed portion of the model. This is designed
  to allow the calculation of standard errors for models for which the \code{predict} method does not provide them.
%%     ~~Describe \code{full} here~~
}
  \item{fixed}{ if \code{Llist} is a character to be used a regular expression, if \code{fixed} is TRUE
  \code{Llist} is interpreted literally, i.e. characters that have a special meaning in regular expressions
  are interpreted literally.
%%     ~~Describe \code{fixed} here~~
}
  \item{invert}{ if \code{Llist} is a character to be used a regular expression, \code{invert == TRUE}
  causes the matches to be inverted so that coefficients that do not match will be selected.
%%     ~~Describe \code{invert} here~~
}
  \item{method}{
     'svd' (current default) or 'qr' is the method used to find the full rank version of the hypothesis matrix.
     'svd' has correctly identified the rank of a large hypothesis matrix where 'qr' has failed.
}
  \item{help}{  obsolete
%%     ~~Describe \code{help} here~~
}
}
\details{
General Linear Hypothesis with Wald test
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
           model.frame(fit) by default. See the example below to estimate a derivative.



}
\value{
An object of class \code{wald}, with the following components
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
\code{\link{Lall}},\code{\link{Lc}},\code{\link{Lequal}},\code{\link{Lform}},
\code{\link{Ldiff}},\code{\link{Lmu}},\code{\link{Lmat}},\code{\link{Lrm}},
\code{\link{as.data.frame.wald}}. To extend to
new models see \code{\link{getFix}}. To generate
hypothesis matrices for general splines see \code{\link{gsp}} and \code{\link{sc}}.
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
# The derivative with respect
# to ses could be evaluated at each point in the following:

data(hs)
require( nlme )
fit <- lme( mathach ~ (ses + I(ses^2)) * Sex, hs, random = ~ 1 + ses | school)

wald( fit, 'Sex')  # sig. overall effect of Sex
wald( fit, ':Sex') # but no evidence of interaction with ses
wald( fit, '\\\\^2') # nor of curvature

# but we continue for the sake of illustration

L <- Lform( fit, list( 0, 1, 2*ses, 0, Sex == 'Male', (Sex == 'Male')*2*ses), hs)
L
(ww <- wald ( fit, L ))
wald.dd <- as.data.frame(ww, se = 2)
head( wald.dd )

require(lattice)
xyplot( coef + U2 + L2 ~ ses | Sex, wald.dd,
 main= 'Increase in predicted mathach per unit increase in ses')

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
