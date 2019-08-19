\name{getFix}
\Rdversion{1.1}
\alias{getFix}
\alias{getFix.default}
\alias{getFix.glm}
\alias{getFix.glmer}
\alias{getFix.lm}
\alias{getFix.lme}
\alias{getFix.lmer}
\alias{getFix.mer}
\alias{getFix.rdc}
\alias{getFix.rdc.lm}
\alias{getFix.rdc.lmer}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get information on fixed effects from a model object
}
\description{
\code{getFix} extracts coefficients, covariance matrix and degrees of
freedom from model objects.
}
\usage{

getFix(fit, ...)

\method{getFix}{default}(fit, ...)

\method{getFix}{glm}(fit, ...)

\method{getFix}{glmer}(fit, ...)

\method{getFix}{lm}(fit, ...)

\method{getFix}{lme}(fit, ...)

\method{getFix}{lmer}(fit, ...)

\method{getFix}{mer}(fit, ...)

\method{getFix}{rdc}(fit, ...)

\method{getFix}{rdc.lm}(fit, ...)

\method{getFix}{rdc.lmer}(fit, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{fit}{
A fitted model object
}
  \item{\dots}{
Other arguments [unused]
}


}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
Returns a list with the following components:
%%  If it is a LIST, use
\item{fixed}{Fixed effect parameter estimates}
\item{vcov}{Covariance matrix of the parameters}
\item{df}{degrees of freedom for each effect}
%% ...
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
\code{\link{wald}}
}
\examples{
fit <- lme( mathach ~ (ses + I(ses^2)) * Sex, hs, random = ~ 1 + ses| school)
getFix(fit)

data(Prestige, package="car")
mod.prestige <- lm(prestige ~ education + income, data=Prestige)
getFix(mod.prestige)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
\keyword{models}% __ONLY ONE__ keyword per line
       