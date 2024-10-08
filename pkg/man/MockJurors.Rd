\name{MockJurors}
\alias{MockJurors}

\title{Confidence of Mock Jurors in Their Verdicts}

\description{
Data with responses of naive mock jurors to the conventional
conventional two-option verdict (guilt vs. acquittal) versus a
three-option verdict setup (the third option was the Scottish
'not proven' alternative), in the presence/absence of conflicting
testimonial evidence.
}

\usage{data("MockJurors", package = "betareg")}

\format{
A data frame containing 104 observations on 3 variables.
\describe{
  \item{verdict}{factor indicating whether a two-option or
    three-option verdict is requested. (A sum contrast rather
    than treatment contrast is employed.)}
  \item{conflict}{factor. Is there conflicting testimonial evidence?
    (A sum contrast rather than treatment contrast is employed.)}
  \item{confidence}{jurors degree of confidence in his/her verdict,
    scaled to the open unit interval (see below).}
}
}

\details{
  The data were collected by Daily (2004) among first-year psychology
  students at Australian National University. Smithson and Verkuilen (2006)
  employed the data scaling the original confidence (on a scale 0--100)
  to the open unit interval: \code{((original_confidence/100) * 103 - 0.5) / 104}.

  The original coding of \code{conflict} in the data provided from Smithson's
  homepage is -1/1 which Smithson and Verkuilen (2006) describe to mean
  no/yes. However, all their results (sample statistics, histograms, etc.)
  suggest that it actually means yes/no which was employed in \code{MockJurors}.
}

\source{
Example 1 from Smithson and Verkuilen (2006) supplements.
}

\references{
Deady S (2004).
  The Psychological Third Verdict: 'Not Proven' or 'Not Willing to Make a Decision'?
  \emph{Unpublished honors thesis}, The Australian National University, Canberra.

Smithson M, Verkuilen J (2006).
  A Better Lemon Squeezer? Maximum-Likelihood Regression with
  Beta-Distributed Dependent Variables.
  \emph{Psychological Methods}, \bold{11}(7), 54--71.
}

\seealso{\code{\link{betareg}}, \code{\link{ReadingSkills}}, \code{\link{StressAnxiety}}}

\examples{
data("MockJurors", package = "betareg")
library("lmtest")

## Smithson & Verkuilen (2006, Table 1)
## variable dispersion model
## (NOTE: numerical rather than analytical Hessian is used for replication,
##  Smithson & Verkuilen erroneously compute one-sided p-values)
mj_vd <- betareg(confidence ~ verdict * conflict | verdict * conflict,
  data = MockJurors, hessian = TRUE)
summary(mj_vd)

## model selection for beta regression: null model, fixed dispersion model (p. 61)
mj_null <- betareg(confidence ~ 1 | 1, data = MockJurors)
mj_fd <-   betareg(confidence ~ verdict * conflict | 1, data = MockJurors)
lrtest(mj_null, mj_fd)
lrtest(mj_null, mj_vd)
## McFadden's pseudo-R-squared
1 - as.vector(logLik(mj_null)/logLik(mj_vd))

## visualization
if(require("lattice")) {
  histogram(~ confidence | conflict + verdict, data = MockJurors,
    col = "lightgray", breaks = 0:10/10, type = "density")
}

## see demo("SmithsonVerkuilen2006", package = "betareg") for more details
}

\keyword{datasets}
