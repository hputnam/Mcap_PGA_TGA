\name{mscm}
\Rdversion{1.1}
\alias{mscm}
\concept{spida 2009}
\concept{longitudinal data}
\concept{dichotomous response}
\concept{observational data}
\concept{time-varying causal factor}
\docType{data}
\title{
Mothers' Stress and Children's Morbidity (MSCM) study: Longitudinal dichotomous response
}
\description{
The Mothers' Stress and Children's Morbidity (MSCM) study is a longitudinal observational study of the causal effect of maternal stress on childhood illness (Zeger and Liang 1986, pp. 125-128).
}
\usage{data(mscm)}
\format{
  A data frame with 5010 observations on 167 mother/child dyads each observed on 30 successive days with the following 14 variables.
  \describe{
    \item{\code{id}}{id for dyad}
    \item{\code{day}}{day of observation: 1 to 30 for each dyad}
    \item{\code{stress}}{indicator for maternal stress}
    \item{\code{illness}}{indicator for child illness}
    \item{\code{married}}{a factor with levels \code{married} \code{other}}
    \item{\code{education}}{a factor with levels \code{less than high school} \code{some high school} \code{high school graduate} \code{some college} \code{college graduate}}
    \item{\code{employed}}{indicator for employment}
    \item{\code{chlth}}{ child health status at baseline:  1=very poor
                                                     2=poor
                                                     3=fair
                                                     4=good
                                                     5=very good}
    \item{\code{mhlth}}{mother health status at baseline: 1=very poor
                                                     2=poor
                                                     3=fair
                                                     4=good
                                                     5=very good}
    \item{\code{race}}{a factor with levels \code{white} \code{non-white}}
    \item{\code{csex}}{a factor with levels \code{male} \code{female}}
    \item{\code{housize}}{a factor with levels \code{2-3 people} \code{more than 3 people}}
    \item{\code{bIllness}}{proportion of child days ill in first 7 days}
    \item{\code{bStress}}{proportion of maternal days with Stress in first 7 days}
  }
}
\details{
The Mothers' Stress and Children's Morbidity (MSCM) study is a longitudinal observational study of the
causal effect of maternal stress on childhood illness (Zeger and Liang 1986, pp. 125-128).
In the MSCM data, the daily prevalence of childhood illness was 14%. Consider the following causal
questions. How would the prevalence change if an ongoing, fully effective stress-reduction intervention
program were instituted? How would prevalence change if conditions worsened and all mothers were subjected
to substantial stress on a daily basis? To attempt to answer these questions, we use a formal model for causal
effects in longitudinal studies introduced by Robins (1986, 1987a,b). This model extends Neyman's (1923) counterfactual causal model for "point" treatment studies to longitudinal studies with time-varying treatments and confounders. We show that the methods for causal inference developed by Robins
provide a better justified basis for answering the foregoing causal questions in longitudinal data in general and in the MSCM study in particular than do methods based on generalized estimating equations (GEE's).
}
\source{
Estimation of the Causal Effect of a Time-Varying Exposure on the Marginal Mean of a Repeated Binary Outcome
Journal article by James M. Robins, Sander Greenland, Fu-Chang Hu; Journal of the American Statistical Association, Vol. 94, 1999
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(mscm)
## maybe str(mscm) ; plot(mscm) ...
}
\keyword{datasets}
