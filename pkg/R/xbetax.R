## extended-domain beta mixture distribution (XBetaX)
## based exponential mixture of the censored symmetric four-parameter beta distribution in regression parameterization
## (mean = mu, precision = phi, latent support = (-Inf, Inf) censored to [0, 1])

## auxiliary quadrature function
quadtable <- function(nquad = 20) {
  matrix(unlist(
    statmod::gauss.quad(n = nquad, kind = "laguerre", alpha = 0)
  ), nrow = nquad)
}

dxbetax <- function(x, mu, phi, nu = 0, log = FALSE, quad = 20) {
  stopifnot(
    "parameter 'mu' must always be in [0, 1]" = all(mu >= 0 & mu <= 1),
    "parameter 'phi' must always be non-negative" = all(phi >= 0),
    "parameter 'nu' must always be non-negative" = all(nu >= 0)
  )
  if(isTRUE(all(nu == 0))) return(dbeta(x, shape1 = mu * phi, shape2 = (1 - mu) * phi, log = log))

  ## unify lengths of all variables
  n <- max(length(x), length(mu), length(phi), length(nu))
  x <- rep_len(x, n)
  mu <- rep_len(mu, n)
  phi <- rep_len(phi, n)
  nu <- rep_len(nu, n)

  ## quadrature
  if(length(quad) == 1L) quad <- quadtable(quad)
  out <- apply(quad, 1, function(rule) {
    e <- rule[1] * nu
    rule[2] * dbeta((x + e)/(1 + 2 * e), shape1 = mu * phi, shape2 = (1 - mu) * phi)/(1 + 2 * e)
  })
  out <- if (is.null(dim(out))) sum(out) else rowSums(out)

  ## censoring
  out[x <= 0] <- pxbetax(0, mu = mu[x <= 0], phi = phi[x <= 0], nu = nu[x <= 0])
  out[x >= 1] <- pxbetax(1, mu = mu[x >= 1], phi = phi[x >= 1], nu = nu[x >= 1], lower.tail = FALSE)
  out[x < 0 | x > 1] <- 0

  ## additional arguments
  if(log) out <- log(out)

  return(out)
}

pxbetax <- function(q, mu, phi, nu = 0, lower.tail = TRUE, log.p = FALSE, quad = 20, censored = TRUE) {
  stopifnot(
    "parameter 'mu' must always be in [0, 1]" = all(mu >= 0 & mu <= 1),
    "parameter 'phi' must always be non-negative" = all(phi >= 0),
    "parameter 'nu' must always be non-negative" = all(nu >= 0)
  )
  if(isTRUE(all(nu == 0))) return(pbeta(q, shape1 = mu * phi, shape2 = (1 - mu) * phi, lower.tail = lower.tail, log.p = log.p))

  ## unify lengths of all variables
  n <- max(length(q), length(mu), length(phi), length(nu))
  q <- rep_len(q, n)
  mu <- rep_len(mu, n)
  phi <- rep_len(phi, n)
  nu <- rep_len(nu, n)

  ## quadrature
  if(length(quad) == 1L) quad <- quadtable(quad)
  out <- apply(quad, 1, function(rule) {
    e <- rule[1] * nu
    rule[2] * pbeta((q + e)/(1 + 2 * e), shape1 = mu * phi, shape2 = (1 - mu) * phi)
  })
  out <- if (is.null(dim(out))) sum(out) else rowSums(out)

  ## tail of the distribution
  if(!lower.tail) out <- 1 - out

  ## censoring
  if (censored) {
    if(lower.tail) {
      out[q <  0] <- 0
      out[q >= 1] <- 1
    } else {
      out[q <= 0] <- 1
      out[q >  1] <- 0
    }
  }

  ## additional arguments
  if(log.p) out <- log(out)
  return(out)
}

qxbetax <- function(p, mu, phi, nu = 0, lower.tail = TRUE, log.p = FALSE, quad = 20) {
  stopifnot(
    "parameter 'mu' must always be in [0, 1]" = all(mu >= 0 & mu <= 1),
    "parameter 'phi' must always be non-negative" = all(phi >= 0),
    "parameter 'nu' must always be non-negative" = all(nu >= 0)
  )
  ## unify lengths of all variables
  n <- max(length(p), length(mu), length(phi), length(nu))
  p <- rep_len(p, n)
  mu <- rep_len(mu, n)
  phi <- rep_len(phi, n)
  nu <- rep_len(nu, n)
  if(length(quad) == 1L) quad <- quadtable(quad)
  obj <- function(pn, mu, phi, nu, p)
      p - pxbetax(qnorm(pn), mu, phi, nu, lower.tail, log.p, quad, censored = FALSE)
  q <- vapply(seq_along(p), function(i) uniroot(obj, c(0, 1), mu = mu[i], phi = phi[i], nu = nu[i], p = p[i])$root, 0.0)
  q <- qnorm(q)
  q[q < 0] <- 0
  q[q > 1] <- 1
  return(q)
}

rxbetax <- function(n, mu, phi, nu = 0) {
  stopifnot(
    "parameter 'mu' must always be in [0, 1]" = all(mu >= 0 & mu <= 1),
    "parameter 'phi' must always be non-negative" = all(phi >= 0),
    "parameter 'nu' must always be non-negative" = all(nu >= 0)
  )
  if(isTRUE(all(nu == 0))) {
    rbeta(n, shape1 = mu * phi, shape2 = (1 - mu) * phi)
  } else {
    rxbeta(n, mu = mu, phi = phi, nu = rexp(n, 1/nu))
  }
}


mean_xbetax <- function(mu, phi, nu, quad = 20, ...) {
    pu <- pxbetax(1, mu, phi, nu, lower.tail = FALSE, quad = quad)
    int <- integrate(function(x) x * dxbetax(x, mu, phi, nu, quad = quad), 0, 1)
    int$value + pu
}

var_xbetax <- function(mu, phi, nu, quad = 20, ...) {
    pu <- pxbetax(1, mu, phi, nu, lower.tail = FALSE, quad = quad)
    int <- integrate(function(x) x^2 * dxbetax(x, mu, phi, nu, quad = quad), 0, 1)
    e <- mean_xbetax(mu, phi, nu)
    int$value + pu - e^2
}


## distributions3 interface

XBetaX <- function(mu, phi, nu = 0) {
  n <- c(length(mu), length(phi), length(nu))
  stopifnot("parameter lengths do not match (only scalars are allowed to be recycled)" = all(n %in% c(1L, max(n))))
  stopifnot(
    "parameter 'mu' must always be in [0, 1]" = all(mu >= 0 & mu <= 1),
    "parameter 'phi' must always be non-negative" = all(phi >= 0),
    "parameter 'nu' must always be non-negative" = all(nu >= 0)
  )
  d <- data.frame(mu = mu, phi = phi, nu = nu)
  class(d) <- c("XBetaX", "distribution")
  d
}

mean.XBetaX <- function(x, ...) {
  m <- vapply(seq_along(x), function(i) mean_xbetax(mu = x$mu[i], phi = x$phi[i], nu = x$nu[i], ...), 0.0)
  setNames(m, names(x))
}

variance.XBetaX <- function(x, ...) {
  v <- vapply(seq_along(x), function(i) var_xbetax(mu = x$mu[i], phi = x$phi[i], nu = x$nu[i], ...), 0.0)
  setNames(v, names(x))
}

skewness.XBetaX <- function(x, ...) {
  stop("not yet implemented")
}

kurtosis.XBetaX <- function(x, ...) {
  stop("not yet implemented")
}

random.XBetaX <- function(x, n = 1L, drop = TRUE, ...) {
  stopifnot(requireNamespace("distributions3"))
  n <- distributions3::make_positive_integer(n)
  if (n == 0L) return(numeric(0L))
  FUN <- function(at, d) rxbetax(n = at, mu = d$mu, phi = d$phi, nu = d$nu)
  distributions3::apply_dpqr(d = x, FUN = FUN, at = n, type = "random", drop = drop)
}

pdf.XBetaX <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  FUN <- function(at, d) dxbetax(x = at, mu = d$mu, phi = d$phi, nu = d$nu, ...)
  distributions3::apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise)
}

log_pdf.XBetaX <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  FUN <- function(at, d) dxbetax(x = at, mu = d$mu, phi = d$phi, nu = d$nu, log = TRUE, ...)
  distributions3::apply_dpqr(d = d, FUN = FUN, at = x, type = "logLik", drop = drop, elementwise = elementwise)
}

cdf.XBetaX <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  FUN <- function(at, d) pxbetax(q = at, mu = d$mu, phi = d$phi, nu = d$nu, ...)
  distributions3::apply_dpqr(d = d, FUN = FUN, at = x, type = "probability", drop = drop, elementwise = elementwise)
}

quantile.XBetaX <- function(x, probs, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  FUN <- function(at, d) qxbetax(p = at, mu = d$mu, phi = d$phi, nu = d$nu, ...)
  distributions3::apply_dpqr(d = x, FUN = FUN, at = probs, type = "quantile", drop = drop, elementwise = elementwise)
}

support.XBetaX <- function(d, drop = TRUE, ...) {
  stopifnot(requireNamespace("distributions3"))
  distributions3::make_support(rep.int(0, length(d)), rep.int(1, length(d)), d, drop = drop)
}

is_discrete.XBetaX <- function(d, ...) {
  setNames(rep.int(FALSE, length(d)), names(d))
}

is_continuous.XBetaX <- function(d, ...) {
  setNames(d$nu <= 0, names(d))
}
