\name{LossAversion}
\encoding{latin1}

\alias{LossAversion}

\title{(No) Myopic Loss Aversion in Adolescents}

\description{
Data for assessing the extent of myopic loss aversion among
adolescents (mostly aged 11 to 19).
}

\usage{data("LossAversion")}

\format{
A data frame containing 570 observations on 7 variables.
\describe{
  \item{invest}{numeric. Average proportion of points invested across
    all 9 rounds.}
  \item{gender}{factor. Gender of the player (or team of players).}
  \item{male}{factor. Was (at least one of) the player(s) male (in the team)?}
  \item{age}{numeric. Age in years (averaged for teams).}
  \item{treatment}{factor. Type of treatment: long vs. short.}
  \item{grade}{factor. School grades: 6-8 (11-14 years) vs. 10-12 (15-18 years).}
  \item{arrangement}{factor. Is the player a single player or team of two?}
}
}

\details{
  Myopic loss aversion is a phenomenon in behavioral economics,
  where individuals do not behave economically rationally when making short-term
  decisions under uncertainty. Example: In lotteries with positive expected payouts
  investments are lower than the maximum possible (loss aversion). This effect
  is enhanced for short-term investments (myopia or short-sightedness).

  The data in \code{LossAversion} were collected by Matthias Sutter and
  Daniela Glätzle-Rützler (Universität Innsbruck) in an experiment with
  high-school students in Tyrol, Austria (Schwaz and Innsbruck). The students
  could invest X points (0-100) in each of 9 rounds in a lottery. The payouts
  were 100 + 2.5 * X points with probability 1/3 and 100 - X points with
  probability 2/3. Thus, the expected payouts were 100 + 1/6 * X points.
  Depending on the treatment in the experiment, the investments could either be
  modified in each round (treatment: "short") or only in round 1, 4, 7
  (treatment "long"). Decisions were either made alone or in teams of two. The
  points were converted to monetary payouts using a conversion of
  EUR 0.5 per 100 points for lower grades (Unterstufe, 6-8) or EUR 1.0 per 100
  points for upper grades (Oberstufe, 10-12).

  From the myopic loss aversion literature (on adults) one would expect that the
  investments of the players (either single players or teams of two) would
  depend on all factors: Investments should be \itemize{
    \item lower in the short treatment (which would indicate myopia),
    \item higher for teams (indicating a reduction in loss aversion),
    \item higher for (teams with) male players,
    \item increase with age/grade.
  }
  See Glätzle-Rützler et al. (2015) for more details and references to the
  literature. In their original analysis, the investments are analyzes using
  a panel structure (i.e., 9 separate investments for each team). Here, the
  data are averaged across rounds for each player, leading to qualitatively
  similar results. The full data along with replication materials are available
  in the Harvard Dataverse.
}

\source{
Glätzle-Rützler D, Sutter M, Zeileis A (2020).
  Replication Data for: No Myopic Loss Aversion in Adolescents? - An Experimental Note.
  \emph{Harvard Dataverse}, UNF:6:6hVtbHavJAFYfL7dDl7jqA==.
  \doi{10.7910/DVN/IHFZAK}
}

\references{
Glätzle-Rützler D, Sutter M, Zeileis A (2015).
  No Myopic Loss Aversion in Adolescents? - An Experimental Note.
  \emph{Journal of Economic Behavior & Organization}, \bold{111}, 169-176.
  \doi{10.1016/j.jebo.2014.12.021}
}

\seealso{\code{\link{betareg}}}

\examples{
## data for students in higher grades
data("LossAversion", package = "betareg")
LossAversion <- subset(LossAversion, grade == "10-12")

## ad hoc scaling (a la Smithson & Verkuilen)
LossAversion$invests <- ((LossAversion$invest * (nrow(LossAversion) - 1) + 0.5)/nrow(LossAversion))

## fraction of boundary observations for ad hoc extended-support beta specification (xbeta)
p01 <- mean(LossAversion$invest <= 0 | LossAversion$invest >= 1)

## main effects models: Gaussian, beta, extended-support beta, extended-support beta mixture
la_gr  <- glm(invest      ~ arrangement + male + age,                      data = LossAversion)
la_br  <- betareg(invests ~ arrangement + male + age | arrangement + male, data = LossAversion)
la_xbx <- betareg(invest  ~ arrangement + male + age | arrangement + male, data = LossAversion)
la_xb  <- betareg(invest  ~ arrangement + male + age | arrangement + male, data = LossAversion,
  dist = "xbeta", nu = p01/2)

}

\keyword{datasets}
